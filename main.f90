        ! LiouVEL: code for computing polar wander and viscoelastic deformation of planetary mantles
        ! For implementation documentation see Patocka et al. (2017) and Mgr thesis at http://geo.mff.cuni.cz/users/patocka/

        ! trel is the time for which the body is initially rotated in order to reach lithostatic equilibrium
        ! tmax is the time of the studied process (total time is trel+tmax)
        ! MAXSTPINI: maximum time steps of initial relaxation; 
        ! MAXSTP: maximum timesteps for the studied process.
        ! model: defines the radial structure of the body that we work with
        ! 1-homogeneous, 2-constant gradient of density, 3-PREM (input file needed), 4-asteroid, 5-discrete layers
        ! loadstate: loading initial state of the body from a file - only the TPW process is then computed (from trel to trel+tmax)
        ! savestate: saves final deformation of the body into inistate.dat, to be loaded in a TPW simulation
        ! couple: applies only during the initial relaxation - 3rd component of the Liouville equation can be switched on (1) or off (0)

        ! ilibini: determines the time stepping criterion in the relaxation phase (pr=0) 
        ! 0-constant dt, 1-Runge-Kutta 5th order used for advancing in time, 2-dt computed from change in displacement vector
        ! ilibproc: determines the solver (and time stepping) of the Liouville equation (evolutionary ODE)
        ! 0-Adams-Bashforth 5th order, 1-Runge-Kutta 5th order, ... for other see inside the code, 5 and 6 need MKL library
        ! memschema: 0-explicit euler, 1-implicit crank-nicholson, 2-crank-nicholson computed iteratively, 3-implicit euler  
        ! ilib: solver for the main system AAx=y
      
      MODULE mConst
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: opt=0, check_solver=1*opt, modd=1*opt, N=100, NN=6*N+2, jmin=2, jmax=2
      INTEGER, SAVE :: ilib=4, ilibini=2, ilibproc=0, MAXSTPINI=1000000, MAXSTP=100000
      INTEGER, SAVE :: selfg=1, savestate=0, loadstate=0, tcross=0 
      INTEGER, SAVE :: couple=0, model=1, memschema=1, qfit=5, sgit=5, eta_mode=1
      REAL, SAVE :: sidday=86164.1                    ! rotation period of the studied body [seconds]
      REAL, SAVE :: trel, tmax, tol0=1.e-11, eps0=1.e-9
      INTEGER, PARAMETER :: Nphi=100, Ntheta=50, kl=7, ku=7, nvar=3, Nvinp=6, nintp=3
      INTEGER, SAVE :: kmax=50000, kmax2=100000, fixload=0, load_shape=0
      REAL,SAVE :: rmax, rmin, dr, sclJ, sclOmega, sclOrbit, omg0, absT0, absH0
      REAL,PARAMETER :: Gconst=6.6732e-11, TINY=1.e-20, pi=3.141592653589793, yr=31558149.5
      REAL, SAVE :: rhocore=10750., rcore_nsync=0., rhoc_nsync=10750.
      REAL, SAVE :: etaref=1.e21, homorho=3000., Gref=7.e10, etamax=1.e29, g0const=9.8, dgrain=10.e-3
      REAL, SAVE :: tgrowth=0.0, kick_ang=0., kick_amp=1., Mpoint=0., SPAf=1., hostmassf=1.e10
      REAL, SAVE :: hice=1500., rhoice=931., cap_width=10., cap_lon=75., cap_col=25., tstzoom=0.0, lithick=0.0
      REAL, SAVE :: r_vinp(Nvinp), rho_inp(Nvinp), visko_inp(Nvinp), shear_inp(Nvinp)
      LOGICAL, SAVE :: sgres=.false., order_analyse=.false., isotest=.false., layered=.false., fosslit=.false.
      LOGICAL, SAVE :: bot_noslip=.false., konst_g=.false., ice00=.false., icetest=.false., sinsmooth=.true.
      LOGICAL, SAVE :: homocore=.true., shell_only=.false., extended=.false., readeta=.false.
      LOGICAL, SAVE :: tides=.false., qfit_flexible=.false., locked=.true., doubrel=.true.
      INTEGER,SAVE :: crash, sourcef=0, ilibtest(0:6)=(/0,0,1,2,4,5,6/), Hu_mode=0, wMIA=0
      COMPLEX :: extload(0:jmax)=-1000., polarcap(0:jmax)
      INTEGER, PARAMETER :: kref=3000000, Npicsavmax=30, mp=0
      CHARACTER(20) :: rule='simpson' !'trapezo' 'simpson'
      
      namelist /switches/ loadstate, savestate, model, konst_g, isotest, layered, couple, extended, &
         tides, Hu_mode, wMIA, locked, bot_noslip
      namelist /params/ MAXSTPINI, MAXSTP, sgit, kmax, kmax2, g0const, sidday, hostmassf, fosslit, doubrel, shell_only
      namelist /tpw_load/ trel, tmax, icetest, tgrowth, kick_ang, Mpoint, sinsmooth, kick_amp, cap_lon, cap_col, &
         hice, rhoice, cap_width, fixload, SPAf, ice00, load_shape, tstzoom
      namelist /solvers/ ilib, ilibini, ilibproc, eps0, qfit, sgit, qfit_flexible
      namelist /model_prof/ eta_mode, rhocore, rmin, rmax, r_vinp, rho_inp, shear_inp, visko_inp, &
         readeta, sourcef, lithick, dgrain, etaref, homorho, Gref, etamax, rcore_nsync, rhoc_nsync

      !eta_mode     0... const viscosity;   1... list of viscosities;   2... exponential increase from etaref to etamax
                  
      END MODULE
    
      MODULE mShared
      USE mConst
      IMPLICIT NONE
      
      COMPLEX :: selfg1(2:N,0:jmax),selfg2(2:N,0:jmax),Memlast(3*N,0:jmax),ylast(NN,0:jmax),rotemp
      COMPLEX :: selfg1Tr(0:jmax,2),selfg2Tr(0:jmax,2),Vpres(0:jmax),Vsurf(0:jmax),PhiIn(2:N,0:jmax),PhiEx(2:N,0:jmax)
      REAL :: tlast,redge(N),rcent(N+1),rho0(N+1),rho0r(N),g0(N+1),A(6,6),B(6,6),Venload(3),tidax(3),eas,loadax(3),lamp
      REAL :: sclBC=1.e-7,sclRK=1.e0,sclREO=1.e0,shear(N),visko(N),Tedge(N),k2t(0:jmax),meanqfit(2)=0.,yode(nvar)
      INTEGER :: j,pr,indx(NN),nstp,nrhse,imprnt,npicsav,sjmp(Nvinp),Nm5l,nsav,qfit_used
      INTEGER :: ldab=2*kl+ku+1,ipiv(NN),info,lwork=-1,lworkj=-1,liwork=-1
      REAL :: w(NN),wmin,wmax,rx(NN),ix(NN),al(NN,kl),nrd,Inertia(3,3),J0(3,3),Isph,Iload(3,3),Ifoss(3,3)
      REAL :: d,trialdtmax,rozjezd=0.*yr,Ball(3,3),Smat(3,3),J0nolit(3,3),hmin
      REAL :: tfakt,taasolve,tliou,tnapln,th1,th2,ur20,k2Te,k2Tf,crotfaktor,m3jump,dc33
      REAL,PARAMETER :: a0=sqrt(2./5.), a1=-sqrt(2./15.), a2=sqrt(1./3.), a3=-sqrt(7./30.)
      REAL,PARAMETER :: b0=-sqrt(3./5.), b1=sqrt(3./15.), b3=sqrt(1./35.), b4=-sqrt(4./7.) 
      LOGICAL :: write_totalE=.false., add_Ifb=.false.
      COMPLEX, ALLOCATABLE :: ymidp(:,:)
      
      REAL :: Mbody, Mcentr, adistsq
      REAL :: Edisip,Esgrav,Edef,Eellost,Eel,Erot,Etid,EdrhoPhi,EdrhoV,Esrc(2)
      REAL :: Erot0,Etid0,Esgrav0,Eel0,EnG0,EdrhoPhi0,EdrhoV0,dJdtw0(3)
      REAL :: SEng,SEng_CMB,Ssg,Ssg_CMB,SEng0,SEng_CMB0,Ssg0,Ssg_CMB0
      REAL :: BoFoEn,BoRoEn,BoTiEn,SuTrEn,CMBTopoEn,CMBSgPrEn,CMBRoPrEn,CMBTiPrEn,BougrEn
      COMPLEX :: ur2(N,0:jmax),ur(N+1,0:jmax),ycheck(NN),ut2(N,0:jmax),Vsurf0,Vsurf0nolit
      
      REAL,ALLOCATABLE :: v(:,:),yhelp(:,:,:),AA(:,:),AAcheck(:,:),ab(:,:),abini(:,:),refer(:,:),current(:,:)
      REAL, ALLOCATABLE :: D3vx(:,:,:), D3vy(:,:,:), D3vz(:,:,:), D3vrad(:,:,:)
      REAL,ALLOCATABLE :: work(:),workj(:)
      INTEGER, ALLOCATABLE :: iwork(:)
      
      END MODULE    
    
      INCLUDE 'mProc.f90'
      
      PROGRAM modul_verze
      USE mConst                ! Shared constants
      USE mShared               ! Shared variables that change over time
      USE nr                    ! Contains ludcmp/ludskb, ODE, and other solvers
      USE mProc                 ! Contains derivs() used by ODE solvers, and supplementary procedures
      
      IMPLICIT NONE
      
      INTEGER :: m,k1,k2,k3
      REAL :: t1,t2,t3,t4,csdef,csref,beta,absH,lovangl
      
      REAL, PARAMETER :: up=1.01, down=1.0, uptol=0.04, downtol=0.07 
      INTEGER :: i,nbad,nok,kount,nstpcurr,nstpref,ilibevo,cutstp,solver
      REAL :: x,dxsav,dxsav2,h,hdid,hnext,xsav,xsav2
      REAL :: dydx(nvar),yscal(nvar),yerr(nvar)
      REAL :: urold, urder, maxwellmax, maxwellmin, chyba(2), solchange, pureL(3,3)
      
      REAL :: param(50),eps,tol,truedtmax
      REAL :: r1,r2
      COMPLEX :: ic=(0.,1.)
      REAL, PARAMETER :: epsrange=1.e6,epsjump=1.4,eps0test=1.e-10,tolrange=1.e5,toljump=epsjump,tol0test=1.e-9
      INTEGER :: ntest,nstp_proceed=0,l2imax=3,lnorm=0,dispatch=0,nroz
      LOGICAL :: ref_curve, hotovo=.false., test=.false., use_urvalue=.false.
      
      open(1,file='param.in',status='old')
      read(1,switches); rewind(1);
      read(1,params); rewind(1);
      read(1,tpw_load); rewind(1);
      read(1,solvers); rewind(1);
      read(1,model_prof); rewind(1);
      close(1)
      call LIve()
      sclOmega = 2.*pi/sidday
     
      IF(isotest) THEN
        print *,'Running ISOTEST, setting couple=0 '
        couple=0
        IF(ilibini==1) THEN; ilibini=2; print *,'Overriding ilibini'; ENDIF;
        IF(tmax/=0.) print *,'WARNING: tmax/=0., but ISOTEST does not work in process pr=1'
      ENDIF
      dr = (rmax-rmin)/(N-1)
      trel = trel*yr; tmax = tmax*yr; tgrowth = tgrowth*yr;

      call cpu_time(t1)
      ! formulae for selfgravity potential are j=2 only, k2t in inistate, and possibly other things, I must check:(..
      IF (jmax/=2 .or. jmin/=2) stop 'jmax must be jmin - code is now designed for j==2 only';      
      IF (extended.and.model/=5) stop 'extended formulation works only for model 5'
      IF(order_analyse) print *,'Warning: power of eq. of motion terms is computed only for order: ', mp

      IF (ilib==2) allocate(ab(ldab,NN),abini(ldab,NN)); IF (ilib==4) allocate(ab(NN,kl+ku+1),abini(NN,kl+ku+1));
      IF (ilib==3) allocate(v(NN,NN)); IF(test) allocate(refer(kref,4),current(kref,4));
      IF (modd/=0) THEN
          allocate(D3vx(2:N,0:Ntheta,0:Nphi),D3vy(2:N,0:Ntheta,0:Nphi),D3vz(2:N,0:Ntheta,0:Nphi),D3vrad(2:N,0:Ntheta,0:Nphi))
          D3vx=0. ; D3vy=0. ; D3vz=0. ; D3vrad=0. ;          
      ENDIF
      IF (check_solver/=0) allocate(AAcheck(NN,NN))
      allocate(yhelp(NN,2,0:jmax))
      IF(ilib/=2.or.opt/=0) allocate(AA(NN,NN))
      IF(memschema==4) allocate(ymidp(NN,0:jmax))                

      ! INITIALIZATION of rcent, redge, rho0, shear, Isph
      shear = Gref
      call initialize_redrec
      call init_rhoGcore
      SELECT CASE(model)
      CASE(1) ! Radially uniform model
        rho0 = homorho
        rho0(N+1) = rhocore-rho0(N)
      CASE(2) ! Linear increase from homorho to 1.5*homorho
        rho0(2:N) = homorho*(/(1.+real(k1-2)/(2*N), k1=2,N)/)
        rho0(N+1) = rhocore-rho0(N)
      CASE(4) ! Asteroid - no core (noslip on the inside)  
        rho0 = homorho
        rho0(N+1) = 0. 
        eta_mode = 0
        bot_noslip = .true.
      END SELECT
      ! rho0 on grid edges for non-layered models
      rho0r(:) = 0.5*(rho0(1:N)+rho0(2:)); 
      rho0r(N) = rho0(N);
      ! For model==5 the above initializations are overridden
      IF(model==5) call layered_model()

      call comp_Mbody
      call initialize_g
      call initialize_visc(loadstate)
      call comp_Iload(Iload)
      IF(fixload==1) cap_lon=0.0
      Venload = (/cos(cap_lon*pi/180.)*sin(cap_col*pi/180.), sin(cap_lon*pi/180.)*sin(cap_col*pi/180.), cos(cap_col*pi/180.)/)
      IF(fixload>0) print *,'Load unit vector (using cap_col, cap_lon) ', Venload
      call manexp(4.*pi*Isph,i)
      sclJ = 10.**(i-1)
      call write_rprofs

      ! Initial condition for the rotation vector
      yode = (/0.,0.,1./)
      IF(isotest) yode = (/0.,0.,0./)

      ! TIDES
      eas = 0.; tidax = 0.; sclOrbit = 0.; adistsq = 0.;
      IF (tides) THEN
        Mcentr = hostmassf * Mbody
        ! When the body is tidally locked, its orbital frequency is the same as its rotational frequency
        sclOrbit = sclOmega
        ! adistsq = ( Gconst*(Mbody+Mcentr) / (sclOrbit**2.) )**(2./3.)
        ! Setting the scale for "equivalent angular speed" of the tidal potential
        eas = sqrt(3.0 * Mcentr/(Mcentr+Mbody)) * sclOrbit
        tidax = (/1.0, 0.0, 0.0/)
      ENDIF
      ! I use the transpose of S from Hu et al. (2017), because S^T has a similar meaning and usage to their Q matrix        
      Smat = reshape((/ 0., 0., -1., 0., 1., 0., 1., 0., 0. /), shape(Smat))

      IF(test) THEN 
        open(26,file='test/presnost.dat'); open (33,file='test/crashstat.dat'); 
        current=0.; refer=0.;
      ENDIF
      IF (loadstate/=1) THEN 
        open(52,file='run/tpw_ini.dat');
        open(54,file='run/urad_ini.dat'); 
        open(66,file='run/enip_ini.dat');
      ELSE
        open(52,file='run/tpw_process.dat');
        open(54,file='run/urad_process.dat');
        open(66,file='run/enip_process.dat');
        IF(icetest) open(99,file='run/icetest.dat')
        IF(icetest.or.fixload>0) open(24,file='run/venus.dat')
      ENDIF      

      Memlast=0.; tlast=0.; npicsav=0; dJdtw0=0.;
      selfg1=0.; selfg2=0.; Vpres=0.; pr=0; k2t=0.;
      maxwellmax=maxval(visko/shear); maxwellmin=minval(visko/shear);     

      IF(.not.isotest) extload=0.      
      print '(A,2e20.10)',' max(maxwell)/trel, min(maxwell)/trel: ',maxwellmax/trel, maxwellmin/trel

      tol=tol0; eps=eps0; cutstp=1; 
      ref_curve=.true.
      ilibevo = merge(ilibtest(0),ilibproc,test)
      IF(wMIA>0) ilibevo=0
      
    DO WHILE (.not.hotovo)   

      dxsav = trel/kmax;  xsav = -2.*dxsav;
      dxsav2 = tmax/kmax2;  xsav2 = -2.*dxsav2;     
      call faktorizace(0.,1,jmax)        
      IF (test) print *,'newrun. eps, tol:  ',eps,tol,'ilibevo, cutstp   ',ilibevo,cutstp
      call cpu_time(t2)
      tfakt=0.; taasolve=0.; tliou=0.; tnapln=0.;
      x=0.;  
      omg0 = absv(yode)
      hdid=0.;  nok=0;  nbad=0;  kount=0;  nsav=0;
      crash=0;  truedtmax=0.;  trialdtmax=0.;
      nrhse=0;  param=0.;
      param(4) = MAXSTP*merge(10,1,ref_curve);
      SuTrEn=0.; CMBTopoEn=0.; CMBSgPrEn=0.; CMBRoPrEn=0.; CMBTiPrEn=0.;
      BoFoEn=0.; BoRoEn=0.; BoTiEn=0.; BougrEn=0.; Edisip=0.; Eellost=0.; Esrc=0.; 

      ! Loading deformation of the body from file inistate.dat
      IF (loadstate==1 .or. (test.and.(.not.ref_curve))) THEN
          IF (fosslit.and.doubrel) stop 'inistate cannot be loaded when fossil bulge is computed via double relaxation'
          call readstate
          print '(a,e15.5)',' READING from time [years] ', tlast/yr
          print *,'fluid limit of degree 2 tidal Love number (k2Tf) ', k2Tf
          print *,'elastic limit of degree 2 tidal Love number (k2Te) ', k2Te
          IF(tmax==0.) THEN; 
            print *,'WARNING: running stage pr=0 again, because loadstate=1 and tmax=0' 
          ELSE; 
            pr=1;            
            IF(tlast/=trel) print *,'WARNING: incorrect load time'
          ENDIF 
      ENDIF
      IF (pr==1) THEN
          tlast = 0.0
          x = 0.0
          Vsurf0 = Vsurf(0) ! Vsurf is not stored in inistate.dat, run without loadstate if Vsurf0 is needed
          omg0 = absv(yode); J0 = Inertia; 
          absH0 = absv(matmul(J0,yode)) + omg0*Mbody/(1.+1./hostmassf)*adistsq / sclJ
          absT0 = absv(matmul(J0,tidax))
          ur20 = real(0.5*(sqrt(2./5.)*(ylast(1,0)+ylast(7,0))-sqrt(3./5.)*(ylast(2,0)+ylast(8,0))))
          print *,'Reading degree two deformation in [m]: ', ur20
          print '(a,e20.10)',' Reading state yode(3) = ', yode(3)
          print '(a,e20.10,a,i2,a,e15.5)',' J0(3,3)-J0(1,1): ', J0(3,3)-J0(1,1), ' x 10^(', i-1,'). Ball(1,1):', Ball(1,1)
          print *,'Non-spherical J0 prior to loading:'
          print '(3e17.7)', J0*sclJ - Ball
          print *,'Estimated free wobble period using fluid and elastic Love numbers [yr]: ',&
           2.*pi /( (J0(3,3)-J0(1,1))/J0(1,1) * (k2Tf-k2Te)/k2Tf * sclOmega )/yr
      ENDIF
      print *,'Angular frequency: ',omg0*sclOmega,' [rad/s]'
      
      IF (.not.ref_curve) THEN
        IF (ilibevo==0) THEN; cutstp=cutstp+1;  
        ELSEIF(ilibevo<5) THEN;   eps=eps*epsjump; 
        ELSE;  tol=tol*toljump; cutstp=cutstp+1; ENDIF;
        IF((eps/eps0test>epsrange).or.(tol/tol0test>tolrange).or.(ilibevo==0.and.cutstp==20)) THEN 
            eps=eps0test;   tol=tol0test;   cutstp=1;
            IF (ntest<=size(ilibtest)-1) THEN; ilibevo=ilibtest(ntest); ntest=ntest+1; ELSE; hotovo=.true.; ENDIF;
        ENDIF;        
      ENDIF

      IF (icetest .and. (isotest.or.pr==1)) THEN
        ! Computing the FULL extload size regardless of tgrowth, designed for the below analysis of Inertia
        call load_with_cap(cap_col,cap_lon,hice,rhoice,cap_width,tgrowth,extload)
        print *,'Loading the body with a surface ice cap. Degree two components in [m] ' 
        print '(2f15.7)', extload
      ENDIF
      ! INERTIA TENSOR CONTRIBUTIONS ANALYSIS
      IF (pr==1) THEN
        call print_eige(J0*sclJ,'HYDROSTATIC BULGE',0,x)
        pureL = Inertia_Tensor(4.*pi*Gconst/5.*rho0(1)*rmax*extload) + Iload
        call print_eige(pureL,'PURE LOAD',0,x)
        call print_eige(pureL,'PURE_LOAD  0',2,x)
        imprnt=0;  call derivs(tgrowth,yode,dydx);  imprnt=0;
        ! The part of bulge that does not adjust immediately to new omega, as balanced by FLOAD + F.B. (i.e. t=0 equilibrium)
        call print_eige((k2Tf-k2Te)/k2Tf*(J0*sclJ - Ball) + (Inertia - J0)*sclJ,'EQUIL_MIA  0',2,x)
        call print_eige((k2Tf-k2Te)/k2Tf*(J0*sclJ - Ball) + (Inertia - J0)*sclJ,'EQUIL MIA',0,x)
        IF (add_Ifb) THEN
            call load_with_cap(0.,0.,hice,rhoice,cap_width,tgrowth,polarcap)
            ! WARNING: Standard definition multiplies the prescribed load by 1+k2^L, i.e. takes the t->inf value
            print *,'pure load Q factor (degree 2 order 0: pure load / the remnant figure) ', &
             -real( 4.*pi*Gconst/5.*rho0(1)*rmax*polarcap(0) / (Vsurf0nolit - Vsurf0))
            call print_eige((Inertia-J0nolit)*sclJ,'LOAD + RESPONSE',0,x)  ! FLOAD = Inertia-J0 - (J0nolit-J0)
            call print_eige((J0nolit-J0)*sclJ,'FOSSIL BULGE',0,x)
            call print_eige((Inertia-J0)*sclJ,'FLOAD + F.B.',0,x)
            call print_eige((Inertia-J0)*sclJ,'FLOAD_+_F.B.  0',2,x)
        ELSEIF (fosslit) THEN
            print *,'WARNING: EQUIL_MIA is not evaluated correctly (k2Tf corresponds to the pr=0 model with no lithosphere)'
            print *,'-- use doubrel=.true. if you want fossil bulge to be computed explicitely (and thus eqMIA correctly)'
        ELSE
            call print_eige((Inertia-J0)*sclJ, 'LOAD + RESPONSE',0,x)
            call print_eige(J0*0.0,'FOSSIL BULGE',0,x)
        ENDIF
        !call print_eige((k2Tf-k2Te)/k2Tf*(J0*sclJ - Ball),'FLUIDO BULGE',0,x)
        call print_eige(Inertia*sclJ-Ball,'INERTIA - BALL',0,x)
      ENDIF

      ! Employing an instantaneous change in |w| due to load emplacement
      IF(pr==1 .and. (icetest.or.fixload>0) .and. wMIA==0 .and. tgrowth==0.0) THEN
        IF(ice00.and.icetest) THEN
            DO k3=1,3
                Ball(k3,k3) = Ball(k3,k3) + 8.*pi/3.*legendre(0,hice,rhoice,cap_width)*rmax**4.
            ENDDO
            print *, 'Icecap mass [kg]: ', 4*pi*legendre(0,hice,rhoice,cap_width)*rmax**2.
        ENDIF
        ! Using the above Inertia (computed directly after loading), avoiding the use of any hard-wired Love numbers
        print *,'ur20 (t=0)', real(0.5*(sqrt(2./5.)*(ylast(1,0)+ylast(7,0))-sqrt(3./5.)*(ylast(2,0)+ylast(8,0))))
        yode(3) = yode(3) * absH0 / absv(matmul(Inertia,yode))
        print *,'Employed jump in omega(3) ', (yode(3)/omg0-1.0)*100, ' [percent]'
      ENDIF

      SELECT CASE(pr)
      CASE(0)
        h = trel/MAXSTPINI;
        IF(ilibini==2) h = maxwellmin/1000;
        nstp_proceed = MAXSTPINI+1 
      CASE(1)
        hmin = tmax/(10*MAXSTP)
        h = 2.*hmin
        nstp_proceed = merge(10*MAXSTP, MAXSTP+1+int(10*MAXSTP*rozjezd/tmax), ref_curve.and.test)
        print *,'Euler frequency [yr] ', sidday*sqrt(J0(1,1)*J0(2,2)/( (J0(3,3)-J0(1,1))*(J0(3,3)-J0(2,2)) )) / yr
        IF ( (abs(kick_amp)/=0.0) .or. (abs(kick_ang)/=0.0) ) THEN
            yode(3) = yode(3)*kick_amp
            yode = yode(3)*(/sin(kick_ang)**2.,sin(kick_ang)*cos(kick_ang),cos(kick_ang)/)
            print '(a,f6.1,a,f6.1,a)',' Kick: |w|=|w| x ', kick_amp, ', tilting w by ', kick_ang, ' deg'
            print '(a,e12.5,a)', ' Initial change of LOD by: ', (1.-yode(3)/omg0)*sidday*1000., ' [ms]'
            IF(wMIA>0) THEN
                print '(a,i1)',' Using one of the LE approximations (wMIA, hybrid, LLE, Hu_et_al), option ', wMIA
                IF(wMIA==4) print *,'Hu et al.: solving LE in the bulge-fixed frame. Simplification (LLE, SLOW_ROTATOR, LLE_EXT): ', Hu_mode
                IF(qfit<=1) STOP 'Set number of quasi-fluid iterations to more than 1'
                ! ilibproc is overwritten above, where ilibevo is set (ilibevo is the actual switch used later in the code)
                IF(ilibproc/=0) print *,'WARNING: setting ilibproc=0, because LE approximations are not combined with high-order solvers'
            ENDIF              
        ENDIF
      END SELECT

      ! MAIN TIME-STEPPING LOOP
      DO nstp=1,nstp_proceed
        ! updating extload in case it grows smoothly
        IF(icetest.and.(x<1.5*tgrowth).and.pr==1) call load_with_cap(cap_col,cap_lon,hice,rhoice,cap_width,x,extload)
                
        imprnt=1
        call derivs(x,yode,dydx);  
        imprnt=0
        
        ! setting a measure for computing the accuracy of obtained solution
        ! yscal=abs(yode)+abs(h*dydx)+TINY
        yscal=(/1.,1.,1./)*sqrt(dot_product(yode,yode)+dot_product(h*dydx,h*dydx))+TINY
        
        SELECT CASE(pr)
        ! HYDROSTATIC RELAXATION of the body, rotation vector changes only if couple=1
        CASE(0) 
            IF ((abs(x-xsav).gt.abs(dxsav)).or.nstp>nsav+9) THEN
                IF (kount.lt.kmax-1) THEN                     
                    kount=kount+1;  xsav=x;  nsav=nstp;
                    IF(x==0.) absH0 = absv(matmul(Inertia,yode)) + omg0*Mbody/(1.+1./hostmassf)*adistsq / sclJ
                    IF(x==0.) absT0 = absv(matmul(Inertia,tidax))
                    !IF(x==0.) print *,'absH0: SPIN vs ORBIT', absv(matmul(Inertia,yode)), omg0*Mbody/(1.+1./hostmassf)*adistsq / sclJ
                    call output(1);
                ELSE
                    print *,'maximum number of records, kmax, exceeded'
                ENDIF
            ENDIF

            ! Initial relaxation has reached requested time
            IF (x==trel) THEN
                ! WARNING: transition from pr=0 to pr=1 is not procedurally smooth: dJdt() is restarted
                IF (savestate==1) call output(4)
                IF(ilibevo==0) h=tmax/MAXSTP;
                call output(5)
                IF(tmax==0.) THEN
                    hotovo=.true.
                ELSEIF (.not.(fosslit.and.doubrel)) THEN
                    close(54); close(66);
                    print *,'-----------------------------------------------------------------------------'
                    print *,'Starting the TPW process' 
                    open(54,file='run/urad_process.dat')
                    open(52,file='run/tpw_process.dat')
                    open(66,file='run/enip_process.dat')
                    IF(icetest) open(99,file='run/icetest.dat')
                    IF(icetest.or.fixload>0) open(24,file='run/venus.dat')
                ENDIF
                IF(fosslit) call initialize_visc(1)
                call cpu_time(t3)                
                exit
            ENDIF 
            
            ! Setting the length of next time step
            IF (ilibini==2) THEN
                IF (mod(nstp,10)==2.and.nstp>9) THEN
                    !call output(3)                    
                    solchange = abs(1.-urder*h/(real(ur2(1,0))-urold))
                    IF(use_urvalue) solchange = abs((real(ur2(1,0))-urold)/urold)
                    IF (solchange < uptol) THEN 
                        hnext = min(h*up, maxwellmin);
                        IF(trel/maxwellmin>MAXSTPINI) hnext = min(h*up, 10.*trel/MAXSTPINI);                    
                    ELSEIF (solchange > downtol) THEN
                        hnext = h/down;
                    ELSE 
                        hnext = h
                    ENDIF
                    urder = (real(ur2(1,mp))-urold)/h;
                ELSE
                    hnext = h; 
                    IF(nstp==2) urder = (real(ur2(1,mp))-urold)/h;
                ENDIF
                h = hnext; 
                urold = real(ur2(1,mp));
            ENDIF
            
            ! Advancing the solution for the rotation vector
            IF( (x+h-trel)*(trel-x) .gt. 0.) h=trel-x    
            SELECT CASE(ilibini)
            CASE(0,2)                
                ! Using R-K for advancing yode (otherwise yode=const), time step is ruled by ilibini
                IF(couple==1.and.wMIA==0) call rkck(yode,dydx,nvar,x,h,yode,yerr,derivs);
                x = x + h; 
                hnext = trel/MAXSTPINI; 
            CASE(1)
                ! Using R-K for advancing the solution yode and also the time x
                call rkqs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
                IF(couple==0) stop 'ilibini=1 is designed for coupled solutions only'
            END SELECT
            IF (ilibini/=2) h = hnext;
            
        ! COMPUTING THE STUDIED TPW PROCESS, Liouville equation is solved
        CASE(1) 
            solver = merge(ilibtest(0), ilibevo, (ref_curve.and.test) .or. x<rozjezd)
            IF(solver==0 .or. solver==5 .or. solver==7) THEN
                hnext = tmax/MAXSTP*cutstp; 
                h = hnext;          
            ENDIF
            ! Saving time step
            IF(tstzoom>0.0) dxsav2 = (tmax/kmax2) / (tstzoom*exp(-tstzoom*x/tmax))
            IF((abs(x-xsav2).gt.abs(dxsav2)) .or. ((x+h-tmax).gt.0.)) THEN
                IF(kount.lt.kmax2-1)THEN
                    kount = kount+1;  
                    xsav2 = x;
                    IF (crash/=1) call output(2)
                    IF (icetest.or.fixload>0) call output(10)
                    call output(11)
                ENDIF
            ENDIF
            IF(test) call output(6);
             
            ! Reducing the time step in order not to overshoot the requested time tmax
            IF((x+h-tmax).gt.0.) h=tmax-x
            IF (h==0.) exit            

            ! ADVANCING THE SOLUTION for the rotation vector IN TIME
            beta=eps; eps=merge(eps0,eps,x<rozjezd);
            SELECT CASE(solver)
            CASE(0)
                ! wMIA runs at a constant time-stepping, derivs is called with imprnt==1 (no ODE solvers, no trial steps)  
                IF (wMIA>0) THEN
                    x = x + h;
                ELSE    ! Adams-Bashforth 5th order scheme 
                    call ab5(yode,dydx,x,h,yscal);                     
                ENDIF
                hdid=h           
            CASE(1)     ! stepper Odeint: rkqs + rkck   
                call rkqs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)         
            CASE(2)     ! stepper Odeint: bsstep + mmid + pzextr        
                call bsstep(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)       
            CASE(3)
                call stiff(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,jac)
            CASE(4)
                call stifbs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,jac)
            CASE(7)     ! explicit Euler scheme
                yode = yode + dydx*h; 
                x = x + h; hdid = h;  
            END SELECT
            eps=beta;
            
            IF(hdid.eq.h)THEN;  nok=nok+1;  ELSE;   nbad=nbad+1;    ENDIF; 
            IF(x>rozjezd.and.nrhse>1) truedtmax=max(hdid,truedtmax)
            IF (x>=tmax) THEN
                print *,'right way to exit. nok: ',nok,' nbad: ',nbad
                print *,'nrhse: ',nrhse,' trialdtmax ',trialdtmax/yr
                IF(wMIA>0) print *,'Mean number of quasi-fluid iterations ',meanqfit(1)/meanqfit(2),', vs qfit= ', qfit
                IF(ilibevo>=5) print *,'nrhse: ',param(35),' trialdtmax ',param(33)/yr,' nstp ',param(34)
                IF(x>tmax) print *,'overshot',x,tmax  
                imprnt=1;   call derivs(tmax,yode,dydx);
                IF(test.and.(ref_curve.or.x==tmax)) call output(6); 
                call cpu_time(t3); 
                IF(ref_curve) THEN 
                    nstpref=nstp+1; ilibevo=ilibtest(1); eps=eps0test; tol=tol0test; ntest=2; 
                ELSE 
                    nstpcurr=merge(nstp+1,nstp,x==tmax); chyba=l2norm(); call output(7); 
                ENDIF 
                ref_curve=.false.; 
                call print_eige(Inertia*sclJ,'TOTALEND',1,x)
                exit
            ENDIF
            IF(x<rozjezd) nroz=nstp
            IF (((.not.ref_curve).and.test.and.x>rozjezd.and.nstp==(MAXSTP-nroz)/100.and.(x-rozjezd)<(tmax-rozjezd)/200).or.crash==1) THEN 
                IF(crash==0) crash=2;   call output(8);
                dispatch=dispatch+1;    crash=3;    exit;           
            ENDIF
             !IF(abs(hnext).lt.hmin) print *,'stepsize smaller than minimum in odeint'
             h=hnext
        END SELECT
      ENDDO ! End of the main time-stepping loop
            
      IF(.not.test.and.pr==1) hotovo=.true.
      IF(ref_curve.and.test.and.pr==1) stop 'nemam referencni krivku'
      IF(x<tmax.and.test.and.crash/=3.and.pr==1) call output(8);
      IF(pr==0.and.x/=0.) print *,'MAXSTPINI exceeded and trel not reached'
      IF(pr==1.and.x<tmax) print *,'MAXSTP exceeded and tmax not reached'  
      IF (fosslit.and.doubrel) THEN
        J0nolit = Inertia
        Vsurf0nolit = Vsurf(0)
        tlast = 0.; Memlast=0.; ylast=0.;
        selfg1=0.; selfg2=0.; Vpres=0.;
        add_Ifb = .true.; doubrel = .false.;
        print *,"Fosslit Doubrel: Running the relaxation stage again"
      ELSE
        pr=1;
      ENDIF
    
    ENDDO    
  !END SUBROUTINE odeint
      
      IF(tmax>0.) THEN
         close(52); close(54); close(66); 
         IF(icetest) close(99)
         IF(icetest.or.fixload>0) close(24)
      ENDIF 
      IF(test) THEN; close(26); close(33); ENDIF;      
      !IF (modd/=0) call zapis3D(N,Ntheta,Nphi)

      call cpu_time(t4)
      
      print *
      print *,'-----------------------------------------------------------'
      print *,'COMPUTATIONAL TIME ANALYSIS'
      print *,'initializing fields ',t2-t1
      print *,'computing solution (+output if on)  ',t3-t2
      print *,'Computing RHS vector ',tnapln-taasolve
      print *,'AAsolver:   ',taasolve,'   factorization: ',tfakt
      print *,'Liouville:  ',tliou    
      print *,'output:  ',t4-t3
      IF(test) print *,'dispatched  ', dispatch,'nstpref  ',nstpref, 'nroz  ',nroz

      !print '(1e27.7)',redge(:),real(ur2(:,0))
      CONTAINS

      FUNCTION rheo_par(r,mode)      ! setting viscosity inside the mantle
      REAL :: rheo_par,r            ! in Pa.s=N.s/m2
      INTEGER :: i,mode
      rheo_par = 0.0
      DO i=1,Nvinp
        IF (r <= r_vinp(i)) THEN
            SELECT CASE(mode)
            CASE(1)            
                rheo_par=visko_inp(i)    
            CASE(2)
                rheo_par=shear_inp(i)
            CASE(3)
                rheo_par=rho_inp(i)
            END SELECT
        ENDIF
      ENDDO
      END FUNCTION

      FUNCTION gfunc(r)
      REAL :: gfunc,r,m
        !units are N/kg
        gfunc = Gconst*mass(r)/r**2.               
        IF (model==1.and.layered) THEN
            m = 4./3.*pi*rhocore*rmin**3. + 4./3.*pi*(rhoc_nsync-rhocore)*rcore_nsync**3.
            m = m + 4./3.*pi*rho0(1)*(min(r,rmax)**3.-rmin**3.)
            print *,'Comparison: ', gfunc, Gconst*m/r**2., (gfunc - Gconst*m/r**2.) / gfunc
            gfunc = Gconst*m/r**2.
        ENDIF        
      END FUNCTION    
      
      FUNCTION rhofunc(r)
      REAL :: rhofunc,r
      INTEGER :: idl
        idl = int(int((rmax-r)/(dr/2))/2)+1
        SELECT CASE(idl)
        CASE(1)
            rhofunc = rho0(idl)+(rho0(idl+1)-rho0(idl))*(redge(idl)-r)/dr/2
        CASE(N)
            rhofunc = rho0(idl)
        CASE DEFAULT
            rhofunc = rho0(idl)+(rho0(idl+1)-rho0(idl))*(rcent(idl)-r)/dr
        END SELECT
      END FUNCTION

      FUNCTION mass(r)                      
      REAL :: mass,r,loss,drh,ko(Nvinp)
      INTEGER :: i,ir,imax=30000     
        ! The mantle is divided into imax layeres, loss is computed for ir layers (according to r)
        loss = 0. 
        drh = (rmax-rmin)/imax
        ir = int((rmax-r)/drh)                  
        DO i=1,ir
            loss = loss + 4.*pi*((rmax-i*drh)**2)*drh*rhofunc(rmax-i*drh)
        ENDDO
        mass = Mbody-loss

        IF(model==5.and.layered) THEN
        ! Overriding the above computed mass with an analytical formula
            DO i=1,Nm5l-1
                ko(i) = 4./3.*pi*(r_vinp(Nm5l-i+1)**3.0) * (rho_inp(Nm5l-i+1)-rho_inp(Nm5l-i))
            ENDDO
            mass = 0.0
            ir = Nm5l
            DO i=1,Nm5l-1
                IF(r>r_vinp(i+1)) THEN
                    mass = mass + ko(Nm5l-i)
                    ir = ir - 1
                ENDIF
            ENDDO
            mass = mass + 4./3.*pi*r**3.*rho_inp(ir)
        ENDIF
      END FUNCTION mass
      
      SUBROUTINE comp_Mbody()
      INTEGER :: ii     
        ! The core mass is already included in subroutine init_rhoGcore
        Mbody = Mbody + 4.*pi*integrate(rcent(2:N),rmin,rmax,rho0(N),rho0(1),rho0(2:N),2)
        Isph = Isph + 2./3.*integrate(rcent(2:N),rmin,rmax,rho0(N),rho0(1),rho0(2:N),4)       
        IF (model==5.and.layered) THEN
        ! Overriding the above Isph with an analytical formula, Mbody is not needed in function mass in this case
            Isph = 2./3.*(rmax**5.0) * rho_inp(1)/5.
            DO i=2,Nm5l
                Isph = Isph + 2./3.*(r_vinp(i)**5.0) * (rho_inp(i)-rho_inp(i-1))/5.    
            ENDDO
        ENDIF
        Ball=0.0
        DO ii=1,3
            Ball(ii,ii) = 4.*pi*Isph
        ENDDO
      END SUBROUTINE comp_Mbody
      
      SUBROUTINE initialize_redrec
      INTEGER :: k1
        !initializing vectors redge(N), rcent(N+1)
        DO k1=1,N
            redge(k1) = rmax - (k1-1)*dr
            rcent(k1) = redge(k1) + dr/2
        ENDDO
        rcent(N+1) = redge(N)-dr/2
      END SUBROUTINE initialize_redrec

      SUBROUTINE initialize_visc(proc)
      INTEGER, INTENT(IN) :: proc
      INTEGER :: k1
      REAL :: Tc,Z,Rgas,Eact,Tbot,Ttop,rphase,Aconst,press
        press = 0
        DO k1=2,N
            ! for g0 and rho0, indices 2,N mark the layer centre values
            press = press + dr*g0(k1+1)*rho0(k1+1)
        ENDDO
        press = press / 1.e6
        
        DO k1=1,N
            SELECT CASE(eta_mode)
            CASE(0)
                visko(k1) = etaref
            CASE(1)
                visko(k1) = rheo_par(redge(k1),1)
            CASE(2)
                visko(k1) = etamax**((redge(k1)-rmin)/(rmax-rmin)) / etaref**((redge(k1)-rmin)/(rmax-rmin)-1.)
            CASE(3)
                ! Pluto
                Tbot = 273.16
                ! Taylor expansion of the melting temperature at the zero pressure depth
                Tbot = -0.000147*press**2 - 0.0748*press + Tbot
                Ttop = 47.; ! 59. (Enceladus)
                Rgas = 8.314 ; Eact=59.e3; Aconst = 9e-8;
                rphase = rmin
                Z = log(Ttop/Tbot) / (rmax-rphase)
                Tedge(k1) = Tbot * exp( Z*rmax*(redge(k1)-rphase) / redge(k1) )
                ! Diffusion creep, Eq. 3.61 in Kihoulou's master thesis
                visko(k1) = min( Tedge(k1)*dgrain**2/(2.0*Aconst)*exp(Eact/(Rgas*Tedge(k1))), etamax)
            END SELECT
            ! Adding elastic LITHOSPHERE
            IF((redge(k1)>rmax-lithick).and..not.(fosslit.and.proc==0)) visko(k1) = etamax*1.e11
        ENDDO
        print '(a,f8.2,a,f8.2)',' Hydrostatic pressure at the bottom of the shell in MPa', press, ', Tbot = ', Tbot
        IF(eta_mode==1 .and. lithick>0.0) print *,'WARNING: lithick will override your visko_inp layers'
        IF(fosslit.and.lithick==0.0) print *,'WARNING: fosslit=.true. but lithosphere thickness is 0'
      END SUBROUTINE initialize_visc

      SUBROUTINE init_rhoGcore           !core treatment reported as troublesome for model==4, check!!
      INTEGER, PARAMETER :: Nrd=1000
      INTEGER :: k1,k2,k3,k4,ierr,ihelp,ihelpp,ni=2
      INTEGER :: idr, idrho, idG, ideta
      REAL :: hlp(10),rhoC(Nrd),rcore(0:Nrd),rmantle(0:Nrd),rhoM(Nrd),Grd(Nrd),etard(Nrd),dierr      
        ! initializes rcore and rmantle, needed only when the input files contains phase transitions (two values for the same radius)
        rcore=0.; rmantle=0.;   
        IF(model==3) THEN
            SELECT CASE(sourcef)
            CASE(0)
                ! prem.csv: 1.radius, 2.depth, 3.density, 4.Vpv, 5.Vph, 6.Vsv, 7.Vsh, 8.eta, 9.Q-mu, 10.Q-kappa(/viscosity for readeta)
                open(55,file='prem.csv');  
                idr=1; idrho=3; idG=9; ideta=10
            CASE(1)
                open(55,file='venus.csv');  
                idr=1; idrho=2; idG=3; ideta=4;
                readeta=.true.
                print *,'Reading venus.csv, readeta is set to .true.'
            END SELECT
                
            k1=1; k2=1; ierr=1;
            DO WHILE(ierr>=0)
                SELECT CASE(sourcef)
                CASE(0)
                    read (55,fmt=*,iostat=ierr) hlp(1),hlp(2),hlp(3),hlp(4),hlp(5),hlp(6),hlp(7),hlp(8),hlp(9),hlp(10)
                    hlp(1)=hlp(1)*1000. 
                    hlp(3)=hlp(3)*1000.
                    hlp(9)=hlp(3)*hlp(6)*hlp(6)*1.e6
                CASE(1)
                    read (55,fmt=*,iostat=ierr) hlp(1),hlp(2),hlp(3),hlp(4)
                END SELECT

                IF (ierr==0) THEN
                    IF(hlp(idr)>rmin) THEN
                        rmantle(k1)=hlp(idr)
                        !dealing with double inputs corresponding to phase transitions in PREM
                        if(rmantle(k1)==rmantle(k1-1)) k1=k1-1      
                        rhoM(k1)=hlp(idrho)
                        Grd(k1)=hlp(idG)
                        etard(k1)=hlp(ideta)
                        k1=k1+1
                    ELSE                    
                        rcore(k2)=hlp(idr);
                        if(rcore(k2)==rcore(k2-1)) k2=k2-1
                        rhoC(k2)=hlp(idrho);
                        k2=k2+1
                    ENDIF
                ENDIF
            ENDDO       
            close(55)
            k1=k1-1; IF(k1==0) stop 'PREM model not read';

            ! Using polint (from Numerical recipes) to interpolate prem.csv quantities onto the grid
            ! For ni=2 this gives simply linear interpolation between each two prem.csv points
            ! Special treatment of the first point is needed for rho0, because rho0(1) is at the cell edge
            call polint(rmantle(1:ni),rhoM(1:ni),ni,redge(1),rho0(1),dierr);
            DO k3=2,N
                DO k4=1,k1;  IF(rcent(k3)>rmantle(k4)) exit;  ENDDO;
                IF(k4<k1/2) call polint(rmantle(k4-1:k4+ni-2),rhoM(k4-1:k4+ni-2),ni,rcent(k3),rho0(k3),dierr);
                IF(k4>=k1/2) call polint(rmantle(k4-ni+1:k4),rhoM(k4-ni+1:k4),ni,rcent(k3),rho0(k3),dierr);
            ENDDO
            ! Doing extrapolation in if rmax>rmantle(1)
            IF( redge(1)>rmantle(1) ) print *,'WARNING: rmax in the input file is smaller than rmax specified in param.in'
            DO k3=1,N
                DO k4=2,k1;  IF(redge(k3)>rmantle(k4)) exit;  ENDDO;                
                IF(k4<k1/2) call polint(rmantle(k4-1:k4+ni-2),Grd(k4-1:k4+ni-2),ni,redge(k3),shear(k3),dierr);
                IF(k4>=k1/2) call polint(rmantle(k4-ni+1:k4),Grd(k4-ni+1:k4),ni,redge(k3),shear(k3),dierr);
                IF(k4<k1/2.and.readeta) call polint(rmantle(k4-1:k4+ni-2),etard(k4-1:k4+ni-2),ni,redge(k3),visko(k3),dierr);
                IF(k4>=k1/2.and.readeta) call polint(rmantle(k4-ni+1:k4),etard(k4-ni+1:k4),ni,redge(k3),visko(k3),dierr);
            ENDDO
        ENDIF     ! end of reading PREM
      
        ! CORE: initialization of its Mbody and Isph contributions
        IF (model/=3.or.homocore) rhoC = rhocore;
        IF (model==4) rhoC = homorho;
        IF (model==3.and.k2==1) stop 'Core not read for the input data file'
        IF (model/=3) THEN
            k2 = Nrd/2 
            rcore(1:k2) = (/(rmin*(1-real(k3)/real(k2)), k3=1,k2)/)
        ENDIF      

        rho0(N+1) = rhoC(1)-rho0(N);
        IF (homocore) THEN  ! Core with a uniform density
            Isph = 2./15.*rmin**5. * rhoC(1)
            Mbody = 4./3.*pi*rmin**3. * rhoC(1)
            IF (rcore_nsync>0.) THEN
                print *,'Assuming rotationally DECOUPLED CORE with radius ',rcore_nsync,' [m] and density ',rhoc_nsync,' [kg/m3]'
                IF(.not.bot_noslip) print *,'overlain by ',rmin-rcore_nsync,'[m] thick layer with density ',rhoC(1),' [kg/m3]'
                Isph = Isph - 2./15.*rhoC(1)*rcore_nsync**5.
                Mbody = Mbody + 4./3.*pi*(rhoc_nsync-rhoC(1))*rcore_nsync**3.
            ENDIF
        ELSE
            Isph = 2./3.*integrate(rcore(2:k2-1),0.,rmin,rhoC(k2),rho0(N)+rho0(N+1),rhoC(2:k2-1),4)
            Mbody = 4.*pi*integrate(rcore(2:k2-1),0.,rmin,rhoC(k2),rho0(N)+rho0(N+1),rhoC(2:k2-1),2)
            IF(rcore_nsync>0.) stop 'Decoupling is implemented when the core is uniform (homocore=.true.)'
        ENDIF
        ! Excluding the core from inertia tensor (only the spherical part, CMB topography is still included)
        IF(shell_only) Isph=0.;
      END SUBROUTINE init_rhoGcore
        
      SUBROUTINE initialize_g
      INTEGER :: k1
        IF (konst_g) THEN
            g0 = g0const
        ELSE
            g0(1) = gfunc(redge(1))
            DO k1=2,N
                g0(k1) = gfunc(rcent(k1)) 
            ENDDO
            g0(N+1) = gfunc(redge(N))
        END IF
        print *,"g_surf ",g0(1), "g_CMB ",g0(N+1)
      END SUBROUTINE initialize_g
      
      SUBROUTINE write_rprofs
      INTEGER :: k1
        open(22,file='run/gridcentres.dat')
        open(23,file='run/gridedges.dat')
        DO k1=1,N
            write(22,'(8e30.17)') merge(redge(1), rcent(k1), k1==1), g0(k1), rho0(k1)
            write(23,'(8e30.17)') redge(k1), shear(k1), visko(k1), Tedge(k1)
        ENDDO
        write(22,'(8e30.17)') redge(N), g0(N+1), rho0(N+1)+rho0(N) 
        close(22); close(23);
      END SUBROUTINE write_rprofs
      
      ! WARNING: going from pr=0 to pr=1 (and restarting) is not smooth: dydx gets overridden, dJdt starts from zero
      SUBROUTINE readstate      
      INTEGER :: ierr,k1,k3n,k5n,k2
      REAL :: hlpr(6)
      COMPLEX :: hlp(0:jmax)
      open(47,file='inistate.dat');  
      k1=1; 
      ierr=1;
      k3n = NN + 3*N
      DO WHILE(ierr>=0)
        hlpr = 0.0
        IF(k1<=k3n) read(47,fmt=*,iostat=ierr) hlpr(1),hlpr(2),hlpr(3),hlpr(4),hlpr(5),hlpr(6)
        IF(k1> k3n) read(47,fmt=*,iostat=ierr) hlpr(1)
        IF (ierr==0) THEN 
            DO k2=0,jmax 
                hlp(k2) = CMPLX(hlpr(2*k2+1),hlpr(2*k2+2))
                IF(k1<=NN             )      ylast(k1,k2)=hlp(k2)
                IF(k1>NN .and. k1<=k3n)      Memlast(k1-NN,k2)=hlp(k2)
            ENDDO
            IF(k1>k3n    .and. k1<=k3n+1)    k2Te=hlpr(1)
            IF(k1>k3n+1  .and. k1<=k3n+2)    k2Tf=hlpr(1)
            IF(k1>k3n+2  .and. k1<=k3n+5)    dydx(k1-k3n-2)=hlpr(1)
            IF(k1>k3n+5  .and. k1<=k3n+8)    yode(k1-k3n-5)=hlpr(1)
            IF(k1>k3n+8  .and. k1<=k3n+17)   Inertia(mod(k1-k3n-9,3)+1,(k1-k3n-9)/3+1)=hlpr(1)
            IF(k1>k3n+17 .and. k1<=k3n+18)   tlast=hlpr(1)
            IF(k1>k3n+18 .and. k1<=k3n+19)   dJdtw0(3)=hlpr(1)
            k1=k1+1
        ENDIF
      ENDDO      
      close(47)
      END SUBROUTINE readstate
      
      FUNCTION l2norm()
      INTEGER, PARAMETER :: nstupen=4
      INTEGER :: i,refi,l
      REAL :: ctverec,ipol,ipolace(3),x(nstupen),y(nstupen),dy,ipolctverec,l2norm(2)          
          refi=1; ctverec=0.;  ipolctverec=0.;
      
          DO i=1,nstpcurr        
            DO WHILE (refer(refi+4,1)<current(i,1)); refi=refi+1; ENDDO;   
            IF (refi+4>nstpref) stop 'jsme za'
            IF (current(i,1)>refer(refi+4,1).or.current(i,1)<refer(refi,1)) stop 'extrapolace'
            DO l=2,l2imax
                x=refer(refi:refi+nstupen-1,1); y=refer(refi:refi+nstupen-1,l)
                call polint(x,y,nstupen,current(i,1),ipol,dy);
                SELECT CASE(lnorm)
                CASE(0)
                    ctverec=max(abs(current(i,l)-ipol),ctverec)
                    ipolctverec=max(abs(dy),ipolctverec)
                CASE(1)
                    ctverec=ctverec+abs(current(i,l)-ipol)
                    ipolctverec=ipolctverec+abs(dy)
                CASE(2)
                    ctverec=ctverec+(current(i,l)-ipol)*(current(i,l)-ipol)
                    ipolctverec=ipolctverec+dy*dy
                END SELECT
                !write (59,'(10e30.17)') current(i,1)/yr,ipol
            ENDDO
            ENDDO
        
          IF (lnorm>0) THEN
            ctverec=ctverec/real(nstpcurr*(l2imax-1)); ipolctverec=ipolctverec/real(nstpcurr*(l2imax-1));
          ENDIF;
          l2norm(1)=ctverec; l2norm(2)=ipolctverec;
      END FUNCTION l2norm
        
      SUBROUTINE output(i)
      COMPLEX :: trakce_rad, trakce_tan, phipot0
      REAL :: rsq, LOD0, drivepot, Qbulge(3,3), loadbg(3)
      REAL, SAVE :: omold(3), alongtrack=0., normaldir, eigen(3), eigenvec(3,3)
      INTEGER :: i,k1,iblb,ihh
      CHARACTER(100) :: labiter
      
      IF(i==1.or.i==2) THEN
        ! Without tides, adist is zero and thus absH reduces to I.w; Note that physical value of absH is scaled by sclJ*sclOmega
        absH = absv(matmul(Inertia,yode)) + absv(yode)*Mbody/(1.+1./hostmassf)*adistsq / sclJ
        rsq = rmax**2.0
        drivepot = (real(phi(2,mp,yode,'rot')+phi(2,mp,tidax,'tid'))*rsq + selfg*real(Vsurf(mp))) / g0(1)
        write(54,'(15e26.16)') x/yr, real(ur2(1,mp)), real(ur2(N,mp)), &
         Erot-Erot0, EnG(1)-EnG0, Eel-Eel0, Edisip, Etid-Etid0, Esrc(1)+Esrc(2), absH, Edef, k2t(mp), drivepot
        write(52,'(10e26.16)') x/yr, yode/omg0, merge(tol,eps,ilibevo>=5), absH/absH0, tidax/omg0
        ! Term-by-term Energy analysis
        write(66,'(22e26.16)') x/yr, -SuTrEn, SEng-SEng0, -CMBTopoEn, SEng_CMB-SEng_CMB0, -BoFoEn, -0.5*(Ssg-Ssg0), &   ! 1-7
         -CMBSgPrEn, -0.5*(Ssg_CMB-Ssg_CMB0), -BoRoEn-CMBRoPrEn, Erot-Erot0, Eellost, Eel-Eel0, Edisip, BougrEn, &      ! 8-15
         EdrhoV-EdrhoV0, EdrhoPhi-EdrhoPhi0, Esgrav-Esgrav0, -BoTiEn-CMBTiPrEn, Etid-Etid0, Esrc(1), Esrc(2)            ! 16-22
        IF(x==0.) THEN
            print *,'Erot, Etid, EnG, Eel: ', Erot, Etid, EnG(1), Eel
            print *
        ENDIF
      ENDIF;
      
      SELECT CASE(i)  
        CASE(1)          
            IF(x==0.) print *,'   time [yr],               driving potential [m],            ur2 [m],              dt [yr],               |w|/|w0|'
            print '(6e26.17)', x/yr, drivepot, real(ur2(1,mp)+extload(mp)), h/yr, absv(yode)/omg0
        CASE(2)                                     
            call eige_problem(Inertia,eigen,eigenvec)
            LOD0 = 2*pi/(omg0*sclOmega)
            ! This output (icetest.dat) is used for the benchmark comparison with Spada et al., 2011
            IF (icetest) THEN
                phipot0 = phi(2,0,(/0.,0.,1./),'rot')
                write(99,'(20e30.20)') x/yr, yode(1:2)/omg0, LOD0*(omg0/yode(3)-1.), real(ur2(1,0)-ur20), real(ur2(1,1)), imag(ur2(1,1)),&
                  real(rsq*phi(2,1,yode,'rot')), imag(rsq*phi(2,1,yode,'rot')), real(Vsurf(1)), imag(Vsurf(1)), real(rsq*(phi(2,0,yode,'rot')-phipot0)),&
                  imag((rsq)*(phi(2,0,yode,'rot')-phipot0)), real(Vsurf(0)-Vsurf0), imag(Vsurf(0)-Vsurf0), eigenvec(1:2,3)
            ENDIF            
            if(x==0.) print *,'time [yr], dt [yr], colatitude, longitude, |w|/|w0|'
            print '(7e20.10)', x/yr, h/yr, col_lon(yode), absv(yode)/omg0
        CASE(3)
            print '(6e26.17)',x/yr,h/yr,abs(urder),abs((real(ur2(1,0))-urold)/h),1.-(urder*h/(real(ur2(1,0))-urold)),urold  
        CASE(4)
            open(53,file='inistate.dat')
            write(53,'(2e30.20)') TRANSPOSE(ylast(:,0:jmax)), TRANSPOSE(Memlast(:,0:jmax))
            write(53,'(1e30.20)') k2Te, k2Tf, dydx, yode, Inertia, tlast, dJdtw0(3)
            close(53)
        CASE(5)
            r1=rmax+real(ur2(1,0)+extload(0))*sqrt(5/pi)/2;
            r2=rmax-real(ur2(1,0)+extload(0))*sqrt(5/pi)/4;    
            print *,'(w-w0)/w0, (H-H0)/H0 [percent] ', (absv(yode)/omg0-1.0)*100., (absH/absH0-1.0)*100.                
            print *,'bilance/maxclen'
            print *,(SuTrEn + CMBTopoEn + CMBSgPrEn + CMBRoPrEn + CMBTiPrEn - (Eel-Eel0+Edisip) + BoFoEn + BoRoEn + BoTiEn)&
                /max(abs(SuTrEn),abs(CMBTopoEn),abs(CMBSgPrEn),abs(CMBRoPrEn),abs(Edisip),abs(Eellost),abs(BoFoEn),abs(BoRoEn),abs(BoTiEn))
            print *,'Inertia diagonal (A,B,C): '
            print *, Inertia(1,1), Inertia(2,2), Inertia(3,3)
            print *,'Inertia non-spherical: '
            print *, Inertia(1,1)-Ball(1,1)/sclJ, Inertia(2,2)-Ball(2,2)/sclJ, Inertia(3,3)-Ball(3,3)/sclJ
        CASE(6)
            IF (nstp<=kref) THEN
                iblb=merge(1,0,x>=tmax)
                IF (ref_curve) THEN; refer(nstp+iblb,1)=x; refer(nstp+iblb,2:4)=yode; 
                ELSE; current(nstp+iblb,1)=x; current(nstp+iblb,2:4)=yode;  ENDIF;
            ELSE
                print *,'There is no space for writing the reference curve', ref_curve
            ENDIF
        CASE(7)
            write(26,'(10e30.17)') real(ilibevo),merge(tol,eps,ilibevo>=5),chyba,t3-t2,real(nstpcurr),truedtmax/yr,real(nrhse)
        CASE(8)         
            write (33,*) crash,merge(tol,eps,ilibevo>=5),ilibevo
        CASE(9)            
            open(35,file="profiler.dat"); open(38,file="profiles.dat");
            DO k1=1,N
              trakce_rad = a0*(a1*ylast(6*k1-3,mp)+a2*ylast(6*k1-2,mp)+a3*ylast(6*k1-1,mp))+b0*(b1*ylast(6*k1-3,mp)+b3*ylast(6*k1-1,mp)+b4*ylast(6*k1,mp))
              trakce_tan =-b0*(a1*ylast(6*k1-3,mp)+a2*ylast(6*k1-2,mp)+a3*ylast(6*k1-1,mp))+a0*(b1*ylast(6*k1-3,mp)+b3*ylast(6*k1-1,mp)+b4*ylast(6*k1,mp))
              write(35,'(10e30.17)') redge(k1),real(trakce_rad),real(trakce_tan)
              ! outputting [\rho_0] as it is computed in equation of motion
              if(k1>1) write(38,'(10e30.17)') rcent(k1), (rho0r(k1)-rho0r(k1-1)), g0(k1), real(PhiIn(k1,0)), real(PhiEx(k1,0)), rho0(k1)
            ENDDO
            close(35); close(38);
        CASE(10)
            if(x>0.) alongtrack = alongtrack + angledeg(yode,omold)
            normaldir = angledeg(yode,Venload)
            omold = yode
            write(24,'(5e26.16)') x/yr, normaldir, alongtrack, col_lon(yode)
        CASE(11)
            IF (tides) THEN
                Qbulge = Umatrix(yode/absv(yode), tidax/absv(tidax))
                loadbg = matmul(TRANSPOSE(Qbulge), loadax)
                open (92,file='run/tidax.dat',Access='append')
                write(92,'(4e24.14,a,i0)') x/yr,absv(tidax),col_lon(tidax),'  iteration  ',qfit_used
                close(92)
            ELSE
                loadbg = loadax
            ENDIF
            open (91,file='run/omega.dat',Access='append')
            write(91,'(7e24.14,a,i0)') x/yr,absv(yode),col_lon(yode),lamp,col_lon(loadbg),'  iteration  ',qfit_used
            close(91)
            write(labiter,"(a,i0)") "iteration  ", qfit_used
            call print_eige(Inertia,trim(labiter),2,x)
      END SELECT
            
      END SUBROUTINE output
      
      SUBROUTINE layered_model()
      IMPLICIT NONE
      INTEGER :: k1,iblb
        DO k1=1,N               
            visko(k1)=rheo_par(redge(k1),1)
            shear(k1)=rheo_par(redge(k1),2)
            rho0(k1)=rheo_par(rcent(k1),3)
            rho0r(k1)=rheo_par(redge(k1),3)
        ENDDO
        rho0(N+1)=rhocore-rho0(N); 
        ! For the initializations above, r_vinp(1) should be > rmax+dr, and r_vinp(Nm5l) should be < rmin
        ! For later operations it is convenient to have these set respectively to rmax and rmin
        r_vinp(1) = rmax
        DO k1=1,Nvinp
            IF(r_vinp(k1) < rmin) exit
        ENDDO
        Nm5l = k1
        IF(layered) print *,'Using optimization for LAYERED models, Nm5l = ', Nm5l
        r_vinp(Nm5l) = rmin
        ! visko_inp(Nm5l) is needed in case initialize_visc is called after r_vinp(Nm5l) is set to rmin
        visko_inp(Nm5l) = visko_inp(Nm5l-1)
        rho_inp(Nm5l) = rhocore

        iblb=1        
        DO k1=2,N
            if(abs(rho0r(k1)-rho0r(k1-1))>0.) then; sjmp(iblb)=k1; iblb=iblb+1; endif;        
        ENDDO        
      END SUBROUTINE layered_model
      
      SUBROUTINE load_with_cap(colat,long,h,rhoi,alpha,time,load)
      IMPLICIT NONE
      REAL, INTENT(IN) :: colat, long, h, rhoi, alpha, time
      COMPLEX, INTENT(OUT) :: load(0:jmax)
        ! degree, order, latitude, longitude
        load(0)=icecap(2,0,colat,long,h,rhoi,alpha);
        load(1)=icecap(2,1,colat,long,h,rhoi,alpha);
        load(2)=icecap(2,2,colat,long,h,rhoi,alpha);
        ! load is treated as a surface topography loading defined in meters
        load = load/rho0(1)
        IF(sinsmooth) THEN
            load=load*merge(1., sin(time/tgrowth*pi/2.), time>=tgrowth)
        ELSE                
            load=load*merge(1., time/tgrowth, time>=tgrowth)                
        ENDIF
      END SUBROUTINE load_with_cap

      FUNCTION IfromC(C20, C21, C22, S21, S22)
      REAL, INTENT(IN) :: C20, C21, C22, S21, S22
      REAL :: IfromC(3,3)
      REAL :: VtoJcst, PhiR20, PhiR21, PhiR22, PhiI21, PhiI22
         VtoJcst=(5.*rmax**3.)/(2*pi*Gconst)
         PhiR20 = C20*Gconst*Mbody/rmax/ sqrt((2.*2+1.)*fac(2-0)/fac(2+0)/4.*pi)
         PhiR21 = C21*Gconst*Mbody/rmax/ sqrt((2.*2+1.)*fac(2-1)/fac(2+1)/4.*pi) /2.0
         PhiR22 = C22*Gconst*Mbody/rmax/ sqrt((2.*2+1.)*fac(2-2)/fac(2+2)/4.*pi) /2.0
         PhiI21 = S21*Gconst*Mbody/rmax/ sqrt((2.*2+1.)*fac(2-1)/fac(2+1)/4.*pi) /2.0
         PhiI22 = S22*Gconst*Mbody/rmax/ sqrt((2.*2+1.)*fac(2-2)/fac(2+2)/4.*pi) /2.0
         IfromC(1,1)=(sqrt(pi/5)*VtoJcst/3*PhiR20-VtoJcst*sqrt(2*pi/15)*PhiR22)
         IfromC(2,2)=(sqrt(pi/5)*VtoJcst/3*PhiR20+VtoJcst*sqrt(2*pi/15)*PhiR22)
         IfromC(3,3)=(-2*sqrt(pi/5)*VtoJcst/3*PhiR20)
         IfromC(1,3)=VtoJcst*sqrt(2*pi/15)*PhiR21
         IfromC(2,3)=-VtoJcst*sqrt(2*pi/15)*PhiI21
         IfromC(1,2)=VtoJcst*sqrt(2*pi/15)*PhiI22
         IfromC(2,1)=IfromC(1,2); IfromC(3,1)=IfromC(1,3); IfromC(3,2)=IfromC(2,3); 
      END FUNCTION IfromC

      SUBROUTINE comp_Iload(Iload)
      IMPLICIT NONE
      REAL :: VtoJcst, C20, C21, C22, S21, S22, C20foss, C21foss, C22foss, S21foss, S22foss
      REAL :: axis(3), eigen(3), eigev(3,3), theta, zunit(3)=(/0.,0.,1./), twistax(3)
      REAL, INTENT(OUT) :: Iload(3,3)
      INTEGER :: indice     
        Iload=0.
        SELECT CASE(fixload)
        CASE(1)             
            ! See the master thesis, Chapter 3.3., only there is a mistake in the general definition of Inertia
            Iload(1,1) = (cos(cap_col*pi/180.))**2.0
            Iload(2,2) = 1.0
            Iload(3,3) = (sin(cap_col*pi/180.))**2.0
            Iload(1,3) = -sin(cap_col*pi/180.)*cos(cap_col*pi/180.)
            Iload(3,1) = Iload(1,3)
            Iload = Iload * (Mpoint * rmax**2.0)
        CASE(2)
            !Moon from Keane and Matsuyama. Using all mass anomalies (spherical disks) as the load:
            Iload = IfromC(-77.e-6, -56.e-6, -16.5e-6, -15.6e-6, 0.3e-6)            
            !Fossil buldge
            Ifoss = IfromC(-126.2e-6, 56.e-6, 38.9e-6, 15.6e-6, -0.3e-6)

            call print_eige(Ifoss,'FOSSILBEFORE',1,0.)
            call print_eige(Iload,'LOADBEFORE',1,0.)
            call print_eige(Iload+Ifoss,'TOTALBEFORE',1,0.)
            call eige_problem(Ifoss, eigen, eigev)
            indice = MAXLOC((eigen), DIM=1)              
            theta = angledeg(eigev(:,indice),zunit)
            print *,'Rotating the fossil bulge and impact basins by ', theta,' degrees'
            axis = cross_product( eigev(:,indice), zunit )
            axis = axis / absv(axis)
            open(94,file='run/principal.dat',Access='append')
                write(94,'(7e24.14,a,a)') 0.,0.,0.,axis(1),axis(2),axis(3),theta, '  AXIS'
            close(94)
            call rotIbyVec(axis, theta, Ifoss)
            call rotIbyVec(axis, theta, Iload)            
            call print_eige(Ifoss,'FOSSILROT',1,0.)
            call print_eige(Iload,'LOADROT',1,0.)
            call print_eige(Iload+Ifoss,'TOTALROT',1,0.)

            Iload = (Iload + Ifoss)*SPAf
            !stop 'Iload called'
        CASE(3)
            print *,'Icetest loading ON, but isostatic compensation is disabled (SURFACE TRACTION form the LOAD is turned OFF)'
            print *,'Iload is not hard-wired: it is a surface topography that GENERATES SELF-GRAVITY. Total mass ice00: ', ice00
            IF(.not.icetest) stop 'WARNING: Fixload=3 means that extload will serve as loading, switch icetest ON'
        CASE(4)
            ! Pluto, Sputnik Planitia
            Iload = IfromC(1379.e-6, 0.0, 0.0, 0.0, 0.0)
            Iload = Iload*SPAf
            theta = cap_col
            twistax = (/cos(cap_lon*pi/180.), sin(cap_lon*pi/180.), 0./)
            print *,'Pluto: Sputnik Planitia defined as inertia tensor contribution, C20 = 1379.e-6, SPAf ', SPAf
            print *,'Inertia tensor contribution of the load:'
            print '(3e16.6)', Iload
            print *,'Rotating the imposed load by ', theta,' degrees'
            call rotIbyVec(twistax, theta, Iload)
        CASE(5)
            call load_with_cap(cap_col,cap_lon,hice,rhoice,cap_width,tgrowth,extload)
            Iload = Inertia_Tensor(4.*pi*Gconst/5.*rho0(1)*rmax*extload)
            print *,'Fixload=5: load_with_cap is used, but SURFACE TRACTION and SELF-GRAVITY of the LOAD are OFF'
            IF(icetest) stop 'Fixload=5 means that cap is treated as Iload, switch icetest OFF'
            extload = 0.
        END SELECT
      END SUBROUTINE comp_Iload

      SUBROUTINE rotIbyVec(vector, theta, Irot)
      ! Rotates the principal axes of I around vector by angle theta
      IMPLICIT NONE
      REAL, INTENT(IN) :: vector(3), theta
      REAL, INTENT(INOUT) :: Irot(3,3)
      REAL, DIMENSION(3,3) :: Rot, diag, eigev, eigevinv
      REAL :: thetarad, c, s, u, v, w, eigen(3), v2amp, vamp
        thetarad = theta * pi/180.          
        c = cos(thetarad)
        s = sin(thetarad)
        u = vector(1)
        v = vector(2)
        w = vector(3)
        v2amp = (u**2+v**2+w**2)
        vamp = sqrt(v2amp)
        
        Rot(1,1) = (u**2 + c*(v**2+w**2))
        Rot(1,2) = (u*v*(1-c) - w*vamp*s)
        Rot(1,3) = (u*w*(1-c) + v*vamp*s)
        Rot(2,1) = (u*v*(1-c) + w*vamp*s)
        Rot(2,2) = (v**2 + c*(u**2+w**2))
        Rot(2,3) = (v*w*(1-c) - u*vamp*s)
        Rot(3,1) = (u*w*(1-c) - v*vamp*s)
        Rot(3,2) = (v*w*(1-c) + u*vamp*s)
        Rot(3,3) = (w**2 + c*(u**2+v**2))
        Rot = Rot / v2amp
        
        call eige_problem(Irot,eigen,eigev)
        diag = 0.        
        diag(1,1) = eigen(1)
        diag(2,2) = eigen(2)
        diag(3,3) = eigen(3)        
        eigev(:,1) = matmul(Rot,eigev(:,1))
        eigev(:,2) = matmul(Rot,eigev(:,2))
        eigev(:,3) = matmul(Rot,eigev(:,3))

        call inverseI(eigev, eigevinv)
        Irot = matmul(matmul(eigev, diag), eigevinv)
      END SUBROUTINE rotIbyVec

      SUBROUTINE inverseI(A,Ainv)
      IMPLICIT NONE
      REAL, INTENT(IN) :: A(3,3)
      REAL :: Ainv(3,3), Det, Com(3,3)
        Det = A(1,1)*A(2,2)*A(3,3) &
            - A(1,1)*A(2,3)*A(3,2) &
            - A(1,2)*A(2,1)*A(3,3) &
            + A(1,2)*A(2,3)*A(3,1) &
            + A(1,3)*A(2,1)*A(3,2) &
            - A(1,3)*A(2,2)*A(3,1)

        Com(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))
        Com(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
        Com(1,3) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))
        Com(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
        Com(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))
        Com(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
        Com(3,1) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))
        Com(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
        Com(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))

        Ainv = TRANSPOSE(Com) / Det
      END SUBROUTINE inverseI
      
      SUBROUTINE LIve()
        print *
        print *,'<------------/-----------'        
        print *,'<--------xxx/------------'
        print *,'<     LIOUSHELL          '
        print *,'<----xooo~/~ooox-----_---'
        print *,'<------xo/ooox------/ \--'
        print *,'<-------/xxx--------\o/--'
        print *,'<------/-------------^-VP'
        print *
      END SUBROUTINE LIve

      END PROGRAM modul_verze
