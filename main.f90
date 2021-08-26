        ! LiouVEL: code for computing polar wander and viscoelastic deformation of planetary mantles
        ! For implementation documentation see Patocka et al. (2017) and Mgr thesis at http://geo.mff.cuni.cz/users/patocka/

        ! tstart is the time for which the body is initially rotated in order to reach lithostatic equilibrium
        ! tmax is the time of the studied process (total time is tstart+tmax)
        ! MAXSTPINI: maximum time steps of initial relaxation; 
        ! MAXSTP: maximum timesteps for the studied process.
        ! model: defines the radial structure of the body that we work with
        ! 1-homogeneous, 2-constant gradient of density, 3-PREM (input file needed), 4-asteroid, 5-discrete layers
        ! loadstate: loading initial state of the body from a file - only the TPW process is then computed (from tstart to tstart+tmax)
        ! savestate: saves final deformation of the body into inistate.dat, to be loaded in a TPW simulation
        ! couple: applies only during the initial relaxation - 3rd component of the Liouville equation can be switched on (1) or off (0)

        ! ilibini: determines the time stepping criterion in the relaxation phase (pr=0) 
        ! 0-constant dt, 1-Runge-Kutta 5th order used for advancing in time, 2-dt computed from change in displacement vector
        ! ilibproc: determines the solver (and time stepping) of the Liouville equation (evolutionary ODE)
        ! 0-Adams-Bashforth 5th order, 1-Runge-Kutta 5th order, ... for other see inside the code, 5 and 6 need MKL library
        ! memschema: 0-explicit euler, 1-implicit crank-nicholson, 2-crank-nicholson computed iteratively, 3-implicit euler  
        ! ilib: solver for the main system AAx=y
        ! ilibb: solver for the inversion of J on the LHS of Liouville eq 
      
      MODULE mConst
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: opt=0, test_reseni=1*opt, modd=1*opt, N=460, NN=6*N+2, jmin=2, jmax=2
      INTEGER, SAVE :: ilib=4, ilibb=1, ilibini=2, ilibproc=0, MAXSTPINI=1000000, MAXSTP=100000
      INTEGER, SAVE :: selfg=1, savestate=0, loadstate=0, tcross=0 
      INTEGER, SAVE :: couple=0, model=1, memschema=1, qfit=5, sgit=5, eta_mode=1
      ! yr is used just for output conversion, sidday is the rotation period of studied body
      REAL, SAVE :: sidday=86164.1
      REAL, SAVE :: tstart, tmax, tol0=1.e-11, eps0=1.e-9
      INTEGER, PARAMETER :: Nphi=100, Ntheta=50, kl=7, ku=7, nvar=3, Nvinp=6, nintp=3
      INTEGER, SAVE :: kmax=50000, kmax2=500000, fixload=0, icemode=0
      REAL,SAVE :: rmax, rmin, dr
      REAL,PARAMETER :: G=6.6732e-11, TINY=1.e-20, pi=3.141592653589793, yr=31536000., hmin=pi*1.e4
      REAL, SAVE :: rhocore=10750., etaref=1.e21, sclJ, homorho=3000., Gref=7.e10, etamax=1.e29, g0const=9.8
      REAL, SAVE :: tgrowth=0.0, kick_ang=0., kick_amp=1., Mpoint=0., SPAf=1., sclOmega, absH0, omega0(3)
      REAL, SAVE :: hice=1500., rhoice=931., cap_width=10., cap_lon=75., cap_col=25.
      REAL, SAVE :: r_vinp(Nvinp), rho_inp(Nvinp), visko_inp(Nvinp), shear_inp(Nvinp)
      LOGICAL, SAVE :: sgres=.false., order_analyse=.false., iceslow=.false., isotest=.false., layered=.false.
      LOGICAL, SAVE :: konst_g=.false., ice00=.false., icetest=.false., CAkor=.false., sinsmooth=.true.
      LOGICAL, SAVE :: homocore=.true., core_wobble=.true., extended=.false., write_prof=.true., readeta=.false.
      LOGICAL, SAVE :: PW_Hu=.true., Hu_Venus=.false.
      INTEGER,SAVE :: crash, ilibtest(0:6)=(/0,0,1,2,4,5,6/), quasi_fluid=0
      COMPLEX :: extload(0:jmax)=-1000.
      INTEGER, PARAMETER :: kref=3000000, Npicsavmax=30, mp=0
      CHARACTER(20) :: rule='simpson' !'trapezo' 'simpson'
      
      namelist /switches/ loadstate, savestate, model, konst_g, eta_mode, isotest, layered, couple
      namelist /params/ MAXSTPINI, MAXSTP, sgit, kmax, kmax2, rhocore, rmin, rmax, &
         etaref, homorho, Gref, etamax, g0const, sidday
      namelist /proces/ tstart, tmax, iceslow, icetest, tgrowth, kick_ang, Mpoint, quasi_fluid, core_wobble, &
         sinsmooth, kick_amp, cap_lon, cap_col, hice, rhoice, cap_width, fixload, SPAf, PW_Hu, Hu_Venus, &
         CAkor, ice00, icemode
      namelist /solvers/ ilib, ilibb, ilibini, ilibproc, eps0
      namelist /model_prof/ r_vinp, rho_inp, shear_inp, visko_inp, readeta

      !eta_mode     0... const viscosity;   1... list of viscosities;   2... exponential increase from etaref to etamax
                  
      END MODULE
    
      MODULE mShared
      USE mConst
      IMPLICIT NONE
      
      COMPLEX :: selfg1(2:N,0:jmax),selfg2(2:N,0:jmax),Memlast(3*N,0:jmax),ylast(NN,0:jmax),rotemp
      COMPLEX :: selfg1Tr(0:jmax,2),selfg2Tr(0:jmax,2),Vpres(0:jmax),Vsurf(0:jmax),PhiIn(2:N,0:jmax),PhiEx(2:N,0:jmax)
      REAL :: tlast,redge(N),rcent(N+1),rho0(N+1),rho0r(N),g0(N+1),A(6,6),B(6,6),Venload(3)
      REAL :: sclBC=1.e-7,sclRK=1.e0,sclREO=1.e0,shear(N),visko(N),Tedge(N),k2t(0:jmax)
      INTEGER :: j,pr,indx(NN),nrot,nstp,nrhse,imprnt,npicsav,sjmp(Nvinp)
      INTEGER :: ldab=2*kl+ku+1,ipiv(NN),info,lwork=-1,lworkj=-1,liwork=-1
      REAL :: w(NN),wmin,wmax,rx(NN),ix(NN),al(NN,kl),nrd,Inertia(3,3),J0(3,3),Ih(3),Iload(3,3), Ifoss(3,3)
      REAL :: eigenumber(3),eigevec(3,3),d,trialdtmax,rozjezd=0.*yr,Ball(3,3)
      REAL :: tfakt,taasolve,tliou,tnapln,th1,th2,ur20,k2Le,k2Te,k2Tf,crotfaktor,m3jump,dc33
      REAL :: a0,b0,a1,a2,a3,b1,b3,b4
      LOGICAL :: write_totalE=.true. 
      COMPLEX, ALLOCATABLE :: ymidp(:,:)
      
      REAL :: Edisip,ErotHarm,Egrav,Edef,Eellost,Eel,Erot,EdrhoPhi,EdrhoV
      REAL :: Erot0,Egrav0,Eel0,ErotHarm0,EnG0,EdrhoPhi0,EdrhoV0,dJdtw0(3)
      REAL :: SEng,SEng_CMB,Ssg,Ssg_CMB,SEng0,SEng_CMB0,Ssg0,Ssg_CMB0
      REAL :: BoFoEn,BoRoEn,SuTrEn,CMBTopoEn,CMBSgPrEn,CMBRoPrEn,BougrEn
      COMPLEX :: Ek(0:jmax),ur2(N,0:jmax),ur(N+1,0:jmax),ycheck(NN),ut2(N,0:jmax),Vsurf0
      
      REAL,ALLOCATABLE :: v(:,:),yhelp(:,:,:),AA(:,:),AAcheck(:,:),ab(:,:),abini(:,:),refer(:,:),current(:,:)
      REAL, ALLOCATABLE :: D3vx(:,:,:), D3vy(:,:,:), D3vz(:,:,:), D3vrad(:,:,:)
      REAL,ALLOCATABLE :: work(:),workj(:)
      INTEGER, ALLOCATABLE :: iwork(:)
      
      END MODULE    
    
      INCLUDE 'mProc.f90'
      
      PROGRAM modul_verze
      USE mConst
      USE mShared
      USE nr        ! Contains ludcmp/ludskb solvers and RK for evolutionary ODE
      USE mProc   
      
      IMPLICIT NONE
      
      INTEGER :: m,k1,k2,k3
      REAL :: t1,t2,t3,t4,csdef,csref,alfa,beta,absH,lovangl
      
      REAL, PARAMETER :: up=1.01, down=1.0, uptol=0.04, downtol=0.07 
      INTEGER :: i,nsav,ido,nbad,nok,kount,nstpcurr,nstpref,ilibevo,cutstp,solver
      REAL :: x,dxsav,dxsav2,h,hdid,hnext=2.*hmin,xsav,xsav2,dydx(nvar),yode(nvar),yscal(nvar),yerr(nvar)
      REAL :: urold, urder, trelaxmax, trelaxmin, chyba(2), solchange
      REAL :: Mearth=0
      
      REAL :: param(50),eps,tol,truedtmax
      REAL :: r1,r2,alfa0,r10,r20,alfa00
      COMPLEX :: ic=(0.,1.)
      REAL, PARAMETER :: epsrange=1.e6,epsjump=1.4,eps0test=1.e-10,tolrange=1.e5,toljump=epsjump,tol0test=1.e-9
      INTEGER :: ntest,nstp_proceed=0,l2imax=3,lnorm=0,dispatch=0,nroz
      LOGICAL :: ref_curve, hotovo=.false., test=.false., use_urvalue=.false.
      
      open(1,file='param.in',status='old')
      read(1,switches); rewind(1);
      read(1,params); rewind(1);
      read(1,proces); rewind(1);
      read(1,solvers); rewind(1);
      read(1,model_prof); rewind(1);
      close(1)
      sclOmega=2.*pi/sidday
      call LIve()
      
      IF(quasi_fluid>0) THEN
        print *
        print *,'Running the QUASI_FLUID approximation, option ', quasi_fluid
        IF(qfit<=1) STOP 'Set number of quasi-fluid iterations to more than 1'
        IF(ilibproc/=0) print *,'WARNING: You will be changing omega inside derivs! NOT TESTED!'        
      ENDIF              
      IF(isotest) THEN
        print *,'Running ISOTEST, setting couple=0 '
        couple=0
        IF(ilibini==1) THEN; ilibini=2; print *,'Overriding ilibini'; ENDIF;
        IF(tmax/=0.) print *,'WARNING: tmax/=0., but ISOTEST does not work in process pr=1'
      ENDIF
      IF(iceslow.and.tgrowth==0.) print *,'WARNING: You should set tgrowth if you want iceslow'
      dr=(rmax-rmin)/(N-1)
      tstart=tstart*yr; tmax=tmax*yr; tgrowth=tgrowth*yr;

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
      IF (test_reseni/=0) allocate(AAcheck(NN,NN))
      allocate(yhelp(NN,2,0:jmax))
      IF(ilib/=2.or.opt/=0) allocate(AA(NN,NN))
      IF(memschema==4) allocate(ymidp(NN,0:jmax))                

      ! Viscosity model
      ! visko_inp=1.e19*(/100000.,1.,10.,100.,1000.,1./)      
      ! r_vinp=(/637.1,627.0,607.0,572.0,537.0,358.0/)*1.e4

      !initialization of rcent, redge, rho0, shear   
      call initialize_rprof                 
      IF (model==1.or.model==2) THEN; rho0=homorho;  rho0(N+1)=rhocore-rho0(N);  shear=Gref;   ENDIF;
      IF (model==2) THEN; rho0(2:N)=homorho*(/(1.+real(k1-2)/(2*N), k1=2,N)/); rho0(N+1)=rhocore-rho0(N); ENDIF;
      ! Asteroid - no core (noslip on the inside)  
      IF (model==4) THEN; rho0=homorho; rho0(N+1)=0.; shear=Gref; visko=etaref;  ENDIF;    

      ! rho0 on grid edges for non-layered models
      rho0r(:)=0.5*(rho0(1:N)+rho0(2:)); 
      rho0r(N)=rho0(N);   

      ! For model==5 the above initializations are overridden
      IF(model==5) call layered_model()
      IF(layered) print *,'Using optimization for LAYERED models'
      call Mbody
      call initialize_g
      call comp_Iload(Iload)
      IF(fixload==1) cap_lon=0.0
      ! Unit sphere coordinates of prescribed loading
      Venload = (/cos(cap_lon*pi/180.)*sin(cap_col*pi/180.), sin(cap_lon*pi/180.)*sin(cap_col*pi/180.), cos(cap_col*pi/180.)/)
      IF(fixload>0) print *,'Load axis ', Venload
      call manexp(4.*pi*Ih(3),i)
      sclJ = 10.**(i-1)

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
      alfa=pi/8; beta=pi/4; trelaxmax=maxval(visko/shear); trelaxmin=minval(visko/shear);        
      a0=sqrt(2./5.); a1=-sqrt(2./15.); a2=sqrt(1./3.); a3=-sqrt(7./30.);
      b0=-sqrt(3./5.); b1=sqrt(3./15.); b3=sqrt(1./35.); b4=-sqrt(4./7.);      

      omega0=(/0.,0.,1./)
      IF(isotest) omega0=0.
      IF(.not.isotest) extload=0.      
      print '(A,2e20.10)',' max(trelax)/tstart, min(trelax)/tstart: ',trelaxmax/tstart, trelaxmin/tstart

      call faktorizace(0.,1,jmax)
            
    !SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)    
      dxsav = tstart/kmax;  dxsav2 = tmax/kmax2;
      xsav = -2.*dxsav;  xsav2 = -2.*dxsav2;     
      tol=tol0; eps=eps0; cutstp=1; ref_curve=.true.
      ilibevo = merge(ilibtest(0),ilibproc,test)
      
    DO WHILE (.not.hotovo)   
        
      IF (test) print *,'newrun. eps, tol:  ',eps,tol,'ilibevo, cutstp   ',ilibevo,cutstp
      call cpu_time(t2)
      tfakt=0.; taasolve=0.; tliou=0.; tnapln=0.;
      x=0.;  
      yode = omega0; 
      h=2.*hmin; 
      IF(ilibini==0) h=tstart/MAXSTPINI;
      crash=0;  truedtmax=0.;   hdid=0.;    nok=0;   nbad=0;   kount=0;   nsav=0;   trialdtmax=0.;
      nrhse=0;  param=0.;   ido=1;  
      param(4) = MAXSTP*merge(10,1,ref_curve);
      Edisip=0.; Eellost=0.; SuTrEn=0.; CMBTopoEn=0.; CMBSgPrEn=0.; CMBRoPrEn=0.; BoFoEn=0.; BoRoEn=0.; BougrEn=0.;

      ! Loading deformation of the body from file inistate.dat
      IF (loadstate==1.or.(test.and.(.not.ref_curve))) THEN; 
          call readstate; x=0.; 
          print '(a,e15.5)',' LOADING from time [years] ', tlast/yr  
          IF(tmax==0.) THEN; 
            print *,'WARNING: running stage pr=0 again, because loadstate=1 and tmax=0' 
          ELSE; 
            pr=1;            
            IF(tlast/=tstart) print *,'WARNING: incorrect load time'
          ENDIF 
          tlast = 0.; h = 2.*hmin;   
          Vpres(1:jmax) = 0.;    
          selfg1(:,1:jmax) = 0.; selfg2(:,1:jmax) = 0.;     
          Vsurf0 = Vsurf(0);
          omega0 = yode; J0 = Inertia; 
          absH0 = absv(matmul(J0,yode))
          ur20 = real(0.5*(sqrt(2./5.)*(ylast(1,0)+ylast(7,0))-sqrt(3./5.)*(ylast(2,0)+ylast(8,0))))
          print *,'Loaded degree two deformation in [m]: ', ur20
          print '(a,e20.10)',' Loaded state yode(3) = ', yode(3)
          print *,'J0(3,3)-J0(1,1): ', J0(3,3)-J0(1,1), ', J0:'
          print '(3f16.11)', J0
          print *,'Euler wobble period [yr]', sidday*sqrt(J0(1,1)*J0(2,2)/( (J0(3,3)-J0(1,1))*(J0(3,3)-J0(2,2)) )) / yr
      ENDIF  
      IF(model==5.and.CAkor.and.pr==1) THEN
        IF(sclJ /= 1.e37) print *,'WARNING: Creating C-A and C correction via different sclJ than expected' 
        DO k3=1,3
            Ball(k3,k3) = Ball(k3,k3) + (8.0394 - J0(3,3))*1.e37; 
            J0(k3,k3) = J0(k3,k3) + (8.0394 - J0(3,3));
        ENDDO
        absH0 = absv(matmul(J0,yode))
        print *,"USING CAkor to fit SPADA et al., 2011. Model C-A vs required: ", (J0(3,3)-J0(1,1))*1.e2, 2.6947
      ENDIF
      
      IF (.not.ref_curve) THEN
        IF (ilibevo==0) THEN; cutstp=cutstp+1;  ELSEIF(ilibevo<5) THEN;   eps=eps*epsjump; 
        ELSE;  tol=tol*toljump; cutstp=cutstp+1; ENDIF;
        IF((eps/eps0test>epsrange).or.(tol/tol0test>tolrange).or.(ilibevo==0.and.cutstp==20)) THEN 
            eps=eps0test;   tol=tol0test;   cutstp=1;
            IF (ntest<=size(ilibtest)-1) THEN; ilibevo=ilibtest(ntest); ntest=ntest+1; ELSE; hotovo=.true.; ENDIF;
        ENDIF;        
      ENDIF
      IF(ilibini==2.and.pr==0) h=trelaxmin/1000;

      IF(icetest .and. (isotest.or.pr==1)) THEN
        call load_with_ice(cap_col,cap_lon,hice,rhoice,cap_width,extload)
        print *,'Loading the body with a surface ice cap. Degree two components in [m] ' 
        print '(2f15.7)', extload
      ENDIF
      ! Employing an instantaneous change in |w| due to load emplacement
      IF(pr==1 .and. (icetest.or.fixload==1) .and. model==5 .and. quasi_fluid==0) THEN
            IF(ice00) THEN
                DO k3=1,3
                    Ball(k3,k3) = Ball(k3,k3) + 8.*pi/3.*legendre(0,hice,rhoice,cap_width)*rmax**4.
                ENDDO
            ENDIF
            if(ice00) print *, 'Icecap mass [kg]: ', 4*pi*legendre(0,hice,rhoice,cap_width)*rmax**2.
            crotfaktor = 2./3.*(J0(3,3)-J0(1,1))/J0(3,3)*k2Te/k2Tf
            dc33 = -4./3.*sqrt(pi/5.)*rmax**4*real(rho0(1)*extload(0))*(1.+k2Le) 
            if(ice00) dc33 = dc33 + 8.*pi/3.*legendre(0,hice,rhoice,cap_width)*rmax**4.
            m3jump = -dc33/(J0(3,3)*sclJ)/(1.-crotfaktor)
            IF(CAkor) THEN
                print *,'discrepancy in m3jump due to the value of C [ms]: ',-(-dc33/(J0(3,3)*sclJ) - m3jump)*sidday*1000.
                ! Spada's dc33 - the sign is corrected, see Patocka et al., 2018
                dc33 = -4./15.*pi*rmax**4.*legendre(2,hice,rhoice,cap_width)*(3.*(cos(25.*pi/180.))**2.-1.)*(1.-0.24398316)  
                print *,'discrepancy in m3jump due to both the value of C and k2Le [ms]: ',-(-dc33/(J0(3,3)*sclJ) - m3jump)*sidday*1000.
                print *,'discrepancy in m3jump due to neglecting crotfaktor [ms]: ',-dc33/(J0(3,3)*sclJ)*crotfaktor*sidday*1000.
                m3jump = -dc33/(J0(3,3)*sclJ)
                if(.not.iceslow) yode(3) = yode(3) + m3jump
            ELSE
                ! Directly computing the Inertia tensor after loading to avoid any use of Love numbers
                imprnt=1;  call derivs(x,yode,dydx);  imprnt=0;
                call print_eige(Inertia - J0, 'LOAD PROPERTIES',.false.)
                print *,'ur20 (t=0)', real(0.5*(sqrt(2./5.)*(ylast(1,0)+ylast(7,0))-sqrt(3./5.)*(ylast(2,0)+ylast(8,0))))
                m3jump = yode(3) * (absH0 / absv(matmul(Inertia,yode)) - 1.0)
                if(.not.iceslow) yode(3) = yode(3) * absH0 / absv(matmul(Inertia,yode))
            ENDIF
            print *,'Employed jump in omega(3) ', m3jump/yode(3)*100, ' [percent]'
            print *,'Free (Chandler) wobble period (using Earth Love numbers) [yr]: ',&
             2.*pi /( (J0(3,3)-J0(1,1))/J0(1,1) * (k2Tf-k2Te)/k2Tf * sclOmega )/yr
      ENDIF

      SELECT CASE(pr)
      CASE(0)
        nstp_proceed = MAXSTPINI+1 
      CASE(1)
        nstp_proceed = merge(10*MAXSTP, MAXSTP+1+int(10*MAXSTP*rozjezd/tmax), ref_curve.and.test)
      END SELECT

      ! MAIN TIME-STEPPING LOOP
      DO nstp=1,nstp_proceed        
        IF(icetest.and.iceslow.and.pr==1) THEN
            call load_with_ice(cap_col,cap_lon,hice,rhoice,cap_width,extload)
            IF(sinsmooth) THEN
                extload=extload*merge(1.,sin(x/tgrowth*pi/2.),x>=tgrowth)
            ELSE                
                extload=extload*merge(1.,x/tgrowth,x>=tgrowth)                
            ENDIF
        ENDIF
        
        IF ( (pr==1 .or. loadstate==1) .and. x==0.) THEN
            yode(3) = yode(3)*kick_amp
            yode = yode(3)*(/sin(kick_ang)**2.,sin(kick_ang)*cos(kick_ang),cos(kick_ang)/)
            print '(a,f6.1,a,f6.1,a)',' Kick: |w|=|w| x ', kick_amp, ', kicking w by ', kick_ang, ' deg'
            print '(a,e12.5,a)', ' Initial change of LOD by: ', (1.-yode(3)/absv(omega0))*sidday*1000., ' [ms]'
        ENDIF           
        
        imprnt=1; 
        ! Computing viscoelastic deformation for the given loading
        call derivs(x,yode,dydx);  
        imprnt=0;
        
        ! setting a measure for computing the accuracy of obtained solution
        ! yscal=abs(yode)+abs(h*dydx)+TINY
        yscal=(/1.,1.,1./)*sqrt(dot_product(yode,yode)+dot_product(h*dydx,h*dydx))+TINY
        
        SELECT CASE(pr)
        ! Initial hydrostatic relaxation of the body, rotation amplitude changes only if couple=1
        CASE(0)
            IF ((abs(x-xsav).gt.abs(dxsav)).or.nstp>nsav+9) THEN
                IF (kount.lt.kmax-1) THEN                     
                    kount=kount+1;  xsav=x;     nsav=nstp;
                    IF (x==0.) THEN
                        r10 = rmax+real(ur2(1,0)+extload(0))*sqrt(5/pi)/2;
                        r20 = rmax-real(ur2(1,0)+extload(0))*sqrt(5/pi)/4; 
                        alfa00 = 2*acos(r10/r20)/sqrt(r20*r20-r10*r10); 
                        absH = absv(matmul(Inertia,yode))
                        absH0 = absH;
                        !call output(9)
                    ENDIF;                     
                    call output(1);
                ELSE
                    print *,'maximum number of records, kmax, exceeded'
                ENDIF
            ENDIF

            ! Initial relaxation has reached requested time
            IF (x==tstart) THEN
                ! WARNING: transition from pr=0 to pr=1 is not smooth because dJdt is set to zero
                IF (savestate==1) call output(4)
                tlast=0.; x=0.; kount=0; h=2*hmin; 
                IF(ilibevo==0) h=tmax/MAXSTP;
                call output(5)
                Edisip = 0.; Eellost = 0.; Vsurf0 = Vsurf(0); 
                omega0 = yode; J0 = Inertia; 
                absH0 = absv(matmul(J0,yode))
                ur20 = real(0.5*(sqrt(2./5.)*(ylast(1,0)+ylast(7,0))-sqrt(3./5.)*(ylast(2,0)+ylast(8,0))))
                close(54); close(66);
                IF(tmax==0.) THEN
                    hotovo=.true.
                ELSE
                    print *,'starting the (second) process' 
                    open(54,file='run/urad_process.dat')
                    open(52,file='run/tpw_process.dat')
                    open(66,file='run/enip_process.dat')
                    IF(icetest) open(99,file='run/icetest.dat')
                    IF(icetest.or.fixload>0) open(24,file='run/venus.dat')
                ENDIF
                call cpu_time(t3)                
                exit
            ENDIF 
            
            ! Setting the length of next time step
            IF (ilibini==2) THEN
                IF (mod(nstp,10)==2.and.nstp>9) THEN
                    !call output(3)                    
                    solchange = abs(1.-urder*h/(real(ur2(1,0))-urold))
                    IF(use_urvalue) solchange = abs((real(ur2(1,0))-urold)/urold)
                    IF(solchange<uptol) THEN 
                        hnext=min(h*up, trelaxmin);
                        IF(tstart/trelaxmin>MAXSTPINI) hnext=min(h*up, 10.*tstart/MAXSTPINI);                    
                    ELSEIF(solchange>downtol) THEN
                        hnext=h/down;
                    ELSE 
                        hnext=h
                    ENDIF
                    urder=(real(ur2(1,mp))-urold)/h;
                ELSE
                    hnext=h; 
                    IF(nstp==2) urder=(real(ur2(1,mp))-urold)/h;
                ENDIF
                h=hnext; urold=real(ur2(1,mp));
            ENDIF
            
            ! Advancing the solution for the rotation vector
            IF( (x+h-tstart)*(tstart-x) .gt. 0.) h=tstart-x    
            SELECT CASE(ilibini)
            CASE(0,2)                
                IF(couple==1) THEN      
                ! Using R-K for advancing the solution, but time step is ruled by ilibini    
                    call rkck(yode,dydx,nvar,x,h,yode,yerr,derivs);
                    x=x+h; 
                    hnext=tstart/MAXSTPINI; 
                ELSE
                ! Rotation vector remains unaltered, only time is advanced                                    
                    x=x+h;
                    hnext=tstart/MAXSTPINI;
                ENDIF;
            CASE(1)
                ! Using R-K for advancing the solution and determining the next time step
                IF(couple==1) call rkqs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
            END SELECT    
            IF (ilibini/=2) h=hnext;
            
        ! COMPUTING THE STUDIED TPW PROCESS, Liouville equation is solved
        CASE(1) 
            solver = merge(ilibtest(0), ilibevo, (ref_curve.and.test) .or. x<rozjezd)
            IF(solver==0 .or. solver==5 .or. solver==7) THEN
                hnext = tmax/MAXSTP*cutstp; 
                h = hnext;          
            ENDIF
            ! Saving time step
            IF(abs(x-xsav2).gt.abs(dxsav2)) THEN
                IF(kount.lt.kmax2-1)THEN                     
                    kount = kount+1;  
                    xsav2 = x;
                    IF (crash/=1) call output(2)
                    IF (icetest.or.fixload>0) call output(10)
                ENDIF
            ENDIF
            IF(test) call output(6);
             
            ! Reducing the time step in order not to overshoot the requested time tmax
            IF((x+h-tmax).gt.0.) h=tmax-x
            IF (h==0.) exit            

            ! Advancing the solution for the rotation vector
            beta=eps; eps=merge(eps0,eps,x<rozjezd);
            SELECT CASE(solver)
            CASE(0)       
                if(quasi_fluid>0) then
                    x = x + h;
                else   ! Adams-Bashforth 5th order scheme 
                    call ab5(yode,dydx,x,h,yscal);                     
                endif
                hdid=h           
            CASE(1)     ! stepper Odeint: rkqs + rkck   
                call rkqs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)         
            CASE(2)     ! stepper Odeint: bsstep + mmid + pzextr        
                call bsstep(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)       
            CASE(3)
                call stiff(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,jac)
            CASE(4)
                call stifbs(yode,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,jac)
            CASE(5)
                alfa=x;
                !call divpag(ido,nvar,fcn,fcnj,Inertia,x,x+h,tol,param,yode)
                hdid=x-alfa;
            CASE(6)
                param(8)=1.; alfa=x;
                DO WHILE (x<tmax)
                    !call divpag(ido,nvar,fcn,fcnj,Inertia,x,tmax,tol,param,yode)
                    IF (ido==5) exit
                ENDDO
                hnext=param(31); hdid=x-alfa;
            CASE(7)
                ! Explicit Euler scheme (not used in the current version)    
                yode = yode + dydx*h; 
                x = x + h; hdid = h;  
            END SELECT
            eps=beta;
            
            IF(hdid.eq.h)THEN;  nok=nok+1;  ELSE;   nbad=nbad+1;    ENDIF; 
            IF(x>rozjezd.and.nrhse>1) truedtmax=max(hdid,truedtmax)
            IF(x>=tmax) THEN
                print *,'right way to exit. nok: ',nok,' nbad: ',nbad
                print *,'nrhse: ',nrhse,' trialdtmax ',trialdtmax/yr
                IF (ilibevo>=5) print *,'nrhse: ',param(35),' trialdtmax ',param(33)/yr,' nstp ',param(34)
                IF (x>tmax) print *,'overshot',x,tmax  
                imprnt=1;   call derivs(tmax,yode,dydx);
                IF(test.and.(ref_curve.or.x==tmax)) call output(6); 
                call cpu_time(t3); 
                IF(ref_curve) THEN 
                    nstpref=nstp+1; ilibevo=ilibtest(1); eps=eps0test; tol=tol0test; ntest=2; 
                ELSE 
                    nstpcurr=merge(nstp+1,nstp,x==tmax); chyba=l2norm(); call output(7); 
                ENDIF 
                ref_curve=.false.; 
                call print_eige(Inertia*sclJ,'TOTALEND',.true.)
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
      ENDDO 
            
      IF (.not.test.and.pr==1) hotovo=.true.
      IF (ref_curve.and.test.and.pr==1) stop 'nemam referencni krivku'
      IF (x<tmax.and.test.and.crash/=3.and.pr==1) call output(8);
      !IF (ido/=1) THEN; ido=3; call divpag(ido,nvar,fcn,fcnj,Inertia,x,x+h,tol,param,yode); ENDIF;
      IF(pr==0.and.x/=0.) print *,'MAXSTPINI exceeded and tstart not reached'
      IF(pr==1.and.x<tmax) print *,'MAXSTP exceeded and tmax not reached'  
      pr=1;
    
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
      print *,'dispatched  ', dispatch,'nstpref  ',nstpref, 'nroz  ',nroz

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
        gfunc=G*mass(r)/(r*r)               
        IF (model==1) THEN
            m=4.*pi*rmin**3.*rhocore/3.+4.*pi*rho0(1)/3.*(min(r,rmax)**3.-rmin**3.)
            gfunc=G*m/r**2.
        ENDIF        
      END FUNCTION    
      
      FUNCTION rhofunc(r)
      REAL :: rhofunc,r
      INTEGER :: kde
        kde=int(int((rmax-r)/(dr/2))/2)+1
        IF (kde==1) rhofunc=rho0(kde)+(rho0(kde+1)-rho0(kde))*(redge(kde)-r)/dr/2
        IF (kde/=1.and.kde/=N) rhofunc=rho0(kde)+(rho0(kde+1)-rho0(kde))*(rcent(kde)-r)/dr
        IF (kde==N) rhofunc=rho0(kde)
      END FUNCTION

      ! The mantle is divided into imax layeres, ubytek proto pocitame jen pro ir vrstev, dle aktualniho r
      FUNCTION mass(r)                      
      REAL :: mass,r,ubytek,drh, ko1,ko2,ko3,ko4,ko5            
      INTEGER :: i,ir,imax=30000     
        ubytek=0.; 
        drh=(rmax-rmin)/imax;      
        ir=int((rmax-r)/drh)                  
        DO i=1,ir
            ubytek=ubytek+4.*pi*((rmax-i*drh)**2)*drh*rhofunc(rmax-i*drh)
        ENDDO
        mass=Mearth-ubytek
        IF(model==5.and.layered) THEN
            ko1=4./3.*pi*rmin**3.*(rhocore-rho_inp(4));          
            ko2=4./3.*pi*r_vinp(4)**3.*(rho_inp(4)-rho_inp(3));
            ko3=4./3.*pi*r_vinp(3)**3.*(rho_inp(3)-rho_inp(2));  
            ko4=4./3.*pi*r_vinp(2)**3.*(rho_inp(2)-rho_inp(1));
            IF(r>r_vinp(2))     THEN; mass= ko1 + ko2 + ko3 + ko4 + 4./3.*pi*r**3.*rho_inp(1);
            ELSEIF(r>r_vinp(3)) THEN; mass= ko1 + ko2 + ko3       + 4./3.*pi*r**3.*rho_inp(2);
            ELSEIF(r>r_vinp(4)) THEN; mass= ko1 + ko2             + 4./3.*pi*r**3.*rho_inp(3);
            ELSEIF(r>rmin     ) THEN; mass= ko1                   + 4./3.*pi*r**3.*rho_inp(4);
            ELSE;                     mass= 4./3.*pi*r**3.*rhocore;                     ENDIF;
        ENDIF
      END FUNCTION mass
      
      SUBROUTINE Mbody
      INTEGER :: i     
        Mearth=Mearth + 4.*pi*integrate(rcent(2:N),rmin,rmax,rho0(N),rho0(1),rho0(2:N),2)
        Ih(3)=Ih(3) + 2./3.*integrate(rcent(2:N),rmin,rmax,rho0(N),rho0(1),rho0(2:N),4)       
        IF(model==5.and.layered) THEN
            Ih(3)=2./3.*rmin**5.*(rhocore-rho_inp(4))/5.
            Ih(3)=Ih(3) + 2./3.*r_vinp(4)**5.*(rho_inp(4)-rho_inp(3))/5.
            Ih(3)=Ih(3) + 2./3.*r_vinp(3)**5.*(rho_inp(3)-rho_inp(2))/5.
            Ih(3)=Ih(3) + 2./3.*r_vinp(2)**5.*(rho_inp(2)-rho_inp(1))/5.
            Ih(3)=Ih(3) + 2./3.*rmax**5.*rho_inp(1)/5.
        ENDIF
        Ball=0.0
        DO i=1,3
            Ball(i,i) = 4*pi*Ih(3)       
        ENDDO
      END SUBROUTINE Mbody
      
      SUBROUTINE initialize_rprof           !core treatment reported as troublesome for model==4, check!!
      INTEGER, PARAMETER :: Nrd=1000
      INTEGER :: k1,k2,k3,k4,ierr,ihelp,ihelpp,ni=2
      REAL :: hlp(11),rhoC(Nrd),rcore(0:Nrd),rmantle(0:Nrd),rhoM(Nrd),Grd(Nrd),etard(Nrd),dierr
      REAL :: Tc,Z,Rgas,Eact,Tbot,Ttop
      
        !initializing vectors redge(N), rcent(N+1), and visko(N)
        DO k1=1,N
            redge(k1)=rmax-(k1-1)*dr
            rcent(k1)=redge(k1)+dr/2
            SELECT CASE(eta_mode)
            CASE(0)
                visko(k1)=etaref
            CASE(1)
                visko(k1)=rheo_par(redge(k1),1)
            CASE(2)
                visko(k1)=etamax**((redge(k1)-rmin)/(rmax-rmin)) / etaref**((redge(k1)-rmin)/(rmax-rmin)-1.)
            CASE(3)
                Tbot=273.; Ttop=59.;
                Rgas=8.314 ; Eact=59.e3;
                Z = log(Ttop/Tbot) / (rmax-rmin)
                Tedge(k1) = Tbot * exp( Z*rmax*(redge(k1)-rmin) / redge(k1) )
                visko(k1) = min( etaref*exp(Eact/(Rgas*Tbot)*(Tbot/Tedge(k1)-1.0)), etamax )
            END SELECT            
        ENDDO
        rcent(N+1)=redge(N)-dr/2; 
      
        ! initializes rcore and rmantle, needed only when the input files contains phase transitions (two values for the same radius)
        rcore=0.; rmantle=0.;   
        IF(model==3) THEN
            open(55,file='prem.csv');  
            k1=1; k2=1; ierr=1;
            DO WHILE(ierr>=0)
            IF(readeta) THEN    
                read (55,fmt=*,iostat=ierr) hlp(1),hlp(2),hlp(3),hlp(4),hlp(5),hlp(6),hlp(7),hlp(8),hlp(9),hlp(10),hlp(11)
            ELSE
                read (55,fmt=*,iostat=ierr) hlp(1),hlp(2),hlp(3),hlp(4),hlp(5),hlp(6),hlp(7),hlp(8),hlp(9),hlp(10)
            ENDIF
            hlp(1)=hlp(1)*1000.; hlp(3)=hlp(3)*1000.;
            IF (ierr==0) THEN
                IF(hlp(1)>rmin) THEN
                    rmantle(k1)=hlp(1)
                    !dealing with double inputs corresponding to phase transitions in PREM
                    if(rmantle(k1)==rmantle(k1-1)) k1=k1-1      
                    rhoM(k1)=hlp(3)
                    Grd(k1)=hlp(3)*hlp(6)*hlp(6)*1.e6
                    if(readeta) etard(k1)=hlp(11)
                    k1=k1+1
                ELSE                    
                    rcore(k2)=hlp(1);
                    if(rcore(k2)==rcore(k2-1)) k2=k2-1
                    rhoC(k2)=hlp(3);
                    k2=k2+1
                ENDIF
            ENDIF
            ENDDO       
            close(55)
            k1=k1-1; IF(k1==0) stop 'PREM model not read';

            ! Using polint (from Numerical recipes) to interpolate prem.csv quantities onto the grid
            ! For ni=2 this gives simply linear interpolation between each two prem.csv points
            ! Special treatment of the first point is needed for rho0, because rho0(1) is at cell edge
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
      
        IF (model/=3.or.homocore) rhoC=rhocore;
        IF (model==4) rhoC=homorho;
        IF (model==3.and.k2==1) stop 'core not read for PREM model'
        IF (model/=3) THEN
            k2=Nrd/2; rcore(1:k2)=(/(rmin*(1-real(k3)/real(k2)), k3=1,k2)/)
        ENDIF
      
        rho0(N+1)=rhoC(1)-rho0(N);
        IF(homocore) THEN
            Mearth=4.*pi*rmin**3.*(rho0(N+1)+rho0(N))/3.
            Ih(3)=2./3.*rmin**5.*(rho0(N+1)+rho0(N))/5.
        ELSE
            Ih(3) = 2./3.*integrate(rcore(2:k2-1),0.,rmin,rhoC(k2),rho0(N)+rho0(N+1),rhoC(2:k2-1),4)
            Mearth = 4.*pi*integrate(rcore(2:k2-1),0.,rmin,rhoC(k2),rho0(N)+rho0(N+1),rhoC(2:k2-1),2)        
        ENDIF
        if(.not.core_wobble) Ih(3)=0.;      !excluding the core from the wobbling process      
      END SUBROUTINE initialize_rprof
        
      SUBROUTINE initialize_g
      INTEGER :: k1
        IF (konst_g) THEN
            g0=g0const
        ELSE
            g0(1)=gfunc(redge(1))
            DO k1=2,N
                g0(k1)=gfunc(rcent(k1)) 
            ENDDO
            g0(N+1)=gfunc(redge(N))
        END IF
        print *,"g_surf ",g0(1), "g_CMB ",g0(N+1)

        IF (write_prof) THEN     
          open(22,file='run/gridcentres.dat'); 
          open(23,file='run/gridedges.dat');
          DO k1=1,N
              write(22,'(8e30.17)') merge(redge(1), rcent(k1), k1==1), g0(k1), rho0(k1)
              write(23,'(8e30.17)') redge(k1), shear(k1), visko(k1), Tedge(k1)
          ENDDO
          write(22,'(8e30.17)') redge(N), g0(N+1), rho0(N+1)+rho0(N) 
          close(22); 
          close(23);
        END IF          
      END SUBROUTINE initialize_g      
      
      ! WARNING: going from pr=0 to pr=1 (and restarting) is not smooth: dydx gets overridden, dJdt starts from zero
      SUBROUTINE readstate      
      INTEGER :: ierr,k1
      REAL :: hlp
      open(47,file='inistate.dat');  
      k1=1; 
      ierr=1;
      DO WHILE(ierr>=0)
        read (47,fmt=*,iostat=ierr) hlp
        IF (ierr==0) THEN 
            IF (k1<=NN) ylast(k1,0)=hlp
            IF (k1>NN       .and.k1<=NN+3*N)    Memlast(k1-NN,0)=hlp
            IF (k1>NN+3*N   .and.k1<=NN+3*N+3)  dydx(k1-NN-3*N)=hlp
            IF (k1>NN+3*N+3 .and.k1<=NN+3*N+6)  yode(k1-NN-3*N-3)=hlp
            IF (k1>NN+3*N+6 .and.k1<=NN+3*N+15) Inertia(mod(k1-NN-3*N-6-1,3)+1,(k1-NN-3*N-6-1)/3+1)=hlp
            IF (k1>NN+3*N+15.and.k1<=NN+3*N+16) tlast=hlp
            IF (k1>NN+3*N+16.and.k1<=NN+4*N+15) selfg1(k1+1-(NN+3*N+16),0)=hlp
            IF (k1>NN+4*N+15.and.k1<=NN+5*N+14) selfg2(k1+1-(NN+4*N+15),0)=hlp
            IF (k1>NN+5*N+14.and.k1<=NN+5*N+15) Vpres(0)=hlp
            IF (k1>NN+5*N+15.and.k1<=NN+5*N+16) k2t(0)=hlp
            IF (k1>NN+5*N+16.and.k1<=NN+5*N+17) Vsurf(0)=hlp
            IF (k1>NN+5*N+17.and.k1<=NN+5*N+18) dJdtw0(3)=hlp
            k1=k1+1
        ENDIF
      ENDDO      
      close(47)
      ylast(:,1:2)=0.; Memlast(:,1:2)=0.;
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
      COMPLEX :: trakce_rad,trakce_tan
      REAL :: colatitude, yodescl(3), rmaxsq, LOD0, drivepot, longitude
      REAL, SAVE :: omold(3), alongtrack=0., normaldir
      INTEGER :: i,k1,iblb,ihh
      
      IF(i==1.or.i==2) THEN
        absH = absv(matmul(Inertia,yode))
        drivepot = (real(Phi(2,mp,yode*SclOmega))*rmaxsq + selfg*real(Vsurf(mp)))/g0(1)
        rmaxsq = rmax**2.0
        write(54,'(15e26.16)') x/yr, real(ur2(1,mp)), real(ur2(N,mp)), Erot-Erot0, EnG(1)-EnG0, Eel-Eel0, Edisip, absH, Edef, k2t(mp), drivepot
        write(52,'(10e26.16)') x/yr, yode/absv(omega0), merge(tol,eps,ilibevo>=5), absH/absH0
        write(66,'(20e26.16)') x/yr, -SuTrEn, SEng-SEng0, -CMBTopoEn, SEng_CMB-SEng_CMB0, -BoFoEn, -0.5*(Ssg-Ssg0), -CMBSgPrEn, -0.5*(Ssg_CMB-Ssg_CMB0),&
                -BoRoEn-CMBRoPrEn, Erot-Erot0, Eellost, Eel-Eel0, Edisip, BougrEn, EdrhoV-EdrhoV0, EdrhoPhi-EdrhoPhi0, Egrav-Egrav0           
        
        IF(x==0.) THEN
            print *,'At t=0: Erot, EnG, Eel ',Erot,EnG(1),Eel
            print *
        ENDIF
      ENDIF;
      
      SELECT CASE(i)  
        CASE(1)          
            trakce_rad = a0*(a1*ylast(3,mp)+a2*ylast(4,mp)+a3*ylast(5,mp)) + b0*(b1*ylast(3,mp)+b3*ylast(5,mp)+b4*ylast(6,mp))
            trakce_tan =-b0*(a1*ylast(3,mp)+a2*ylast(4,mp)+a3*ylast(5,mp)) + a0*(b1*ylast(3,mp)+b3*ylast(5,mp)+b4*ylast(6,mp))                       
           
            if(x==0.) print *,'   time [yr],     driving potential [m],        ur2 [m],        dt [yr]'
            print '(6e26.17)', x/yr, drivepot, real(ur2(1,mp)+extload(mp)), h/yr
        CASE(2)                                     
            call eige_problem(Inertia,eigenumber,eigevec,1)
            yodescl = yode*sclOmega
            LOD0 = 2*pi/(absv(omega0)*sclOmega)
            if(icetest) write(99,'(20e30.20)') x/yr, yode(1:2)/absv(omega0), LOD0*(absv(omega0)/yode(3)-1.),&
               real(ur2(1,0)-ur20), real(ur2(1,1)), imag(ur2(1,1)), real((rmaxsq)*phi(2,1,yodescl)), imag((rmaxsq)*phi(2,1,yodescl)),&
               real(Vsurf(1)), imag(Vsurf(1)), real((rmaxsq)*(phi(2,0,yodescl)-phi(2,0,(/0.,0.,1./)*sclOmega))),&
               imag((rmaxsq)*(phi(2,0,yodescl)-phi(2,0,(/0.,0.,1./)*sclOmega))), real(Vsurf(0)-Vsurf0), imag(Vsurf(0)-Vsurf0), eigevec(1:2,3)
            
            if(x==0.) print *,'time [yr], dt [yr], colatitude, longitude, w_3'
            !colatitude = 90.0 - atan(yode(3)/sqrt(yode(1)**2+yode(2)**2))/pi*180
            colatitude = acos(yode(3)/absv(yode))/pi*180
            longitude = 0.
            if(abs(yode(1))>0.) longitude = (atan(yode(2)/yode(1))/pi*180)
            print '(7e20.10)', x/yr, h/yr, colatitude, longitude, yode(3)/absv(omega0)
        CASE(3)
            print '(6e26.17)',x/yr,h/yr,abs(urder),abs((real(ur2(1,0))-urold)/h),1.-(urder*h/(real(ur2(1,0))-urold)),urold  
        CASE(4)
            open(53,file='inistate.dat')
            write(53,'(1e30.20)') real(ylast(:,0)),real(Memlast(:,0)),dydx,yode,Inertia,tlast,real(selfg1(:,0)),real(selfg2(:,0)),&
                                  real(Vpres(0)),real(k2t(0)),real(Vsurf(0)),dJdtw0(3)
            close(53)
        CASE(5)
            r1=rmax+real(ur2(1,0)+extload(0))*sqrt(5/pi)/2;
            r2=rmax-real(ur2(1,0)+extload(0))*sqrt(5/pi)/4;    
            alfa0=merge(2*acos(r1/r2)/sqrt(r2*r2-r1*r1),2/r2,abs((r2-r1)/r2)>1.e-10)
            print *,'w/w0-1, H/H0-1  ',absv(yode)/absv(omega0)-1.,absH/absH0-1.                
            print *,'bilance/maxclen'
            print *,(SuTrEn + CMBTopoEn + CMBSgPrEn + CMBRoPrEn - (Eel-Eel0+Edisip) + BoFoEn + BoRoEn)&
                /max(abs(SuTrEn),abs(CMBTopoEn),abs(CMBSgPrEn),abs(CMBRoPrEn),abs(Edisip),abs(Eellost),abs(BoFoEn),abs(BoRoEn))
            print *,'Inertia diagonal: '
            print *,Inertia(3,3), Inertia(1,1), Inertia(2,2)
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
            do k1=1,N
                trakce_rad = a0*(a1*ylast(6*k1-3,mp)+a2*ylast(6*k1-2,mp)+a3*ylast(6*k1-1,mp)) &
                    + b0*(b1*ylast(6*k1-3,mp)+b3*ylast(6*k1-1,mp)+b4*ylast(6*k1,mp))
                trakce_tan =-b0*(a1*ylast(6*k1-3,mp)+a2*ylast(6*k1-2,mp)+a3*ylast(6*k1-1,mp)) &
                    + a0*(b1*ylast(6*k1-3,mp)+b3*ylast(6*k1-1,mp)+b4*ylast(6*k1,mp))
                write(35,'(10e30.17)') redge(k1),real(trakce_rad),real(trakce_tan)
                ! outputting [\rho_0] as it is computed in equation of motion
                if(k1>1) write(38,'(10e30.17)') rcent(k1),(rho0r(k1)-rho0r(k1-1)),g0(k1),&
                    real(PhiIn(k1,0)),real(PhiEx(k1,0)),rho0(k1)
            enddo
            close(35); close(38);
        CASE(10)
            if(x>0.) alongtrack = alongtrack + acos(dot_product(yode, omold) / (absv(yode)*absv(omold))) * 180./pi
            normaldir = acos(dot_product(yode, Venload) / absv(yode)) * 180./pi
            omold = yode
            write(24,'(3e26.16)') x/yr, normaldir, alongtrack 
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
        iblb=1;
        
        DO k1=2,N
            if(abs(rho0r(k1)-rho0r(k1-1))>0.) then; sjmp(iblb)=k1; iblb=iblb+1; endif;        
        ENDDO
        
        ! Love numbers determined in relaxation runs for model Spada
        SELECT CASE(N)
        CASE(100)          
            !model s mekkou litosferou
            k2Le=-0.243988441390229; k2Te=0.303380573895915; k2Tf=0.975935482690866;
            !model s tuhou litosferou
            k2Le=-0.243988441390229; k2Te=0.303380573895915; k2Tf=0.961946158868269;
        CASE(300)
            !model s mekkou litosferou
            k2Le=-0.244424717369092; k2Te=0.303356238198566; k2Tf=0.973337019158345;           
            !model s tuhou litosferou
            k2Le=-0.244424717369092; k2Te=0.303356238198564; k2Tf=0.959426119913619;
        CASE default                
            ! Spada2011
            k2Le=-0.24398316;        k2Te=0.30346466;        k2Tf=0.96672389;
        END SELECT          
      END SUBROUTINE layered_model
      
      SUBROUTINE load_with_ice(colat,long,h,rhoi,alpha,load)
      IMPLICIT NONE
      REAL, INTENT(IN) :: colat, long, h, rhoi, alpha
      COMPLEX, INTENT(OUT) :: load(0:jmax)
        ! degree, order, latitude, longitude
        load(0)=icecap(2,0,colat,long,h,rhoi,alpha);
        load(1)=icecap(2,1,colat,long,h,rhoi,alpha);
        load(2)=icecap(2,2,colat,long,h,rhoi,alpha);
        ! load is treated as a surface topography loading defined in meters 
        load = load/rho0(1)
      END SUBROUTINE load_with_ice

      FUNCTION IfromC(C20, C21, C22, S21, S22)
      REAL, INTENT(IN) :: C20, C21, C22, S21, S22
      REAL :: IfromC(3,3)
      REAL :: VtoJcst, PhiR20, PhiR21, PhiR22, PhiI21, PhiI22
         VtoJcst=(5.*rmax**3.)/(2*pi*G)
         PhiR20 = C20*G*Mearth/rmax/ sqrt((2.*2+1.)*fac(2-0)/fac(2+0)/4.*pi)
         PhiR21 = C21*G*Mearth/rmax/ sqrt((2.*2+1.)*fac(2-1)/fac(2+1)/4.*pi) /2.0
         PhiR22 = C22*G*Mearth/rmax/ sqrt((2.*2+1.)*fac(2-2)/fac(2+2)/4.*pi) /2.0
         PhiI21 = S21*G*Mearth/rmax/ sqrt((2.*2+1.)*fac(2-1)/fac(2+1)/4.*pi) /2.0
         PhiI22 = S22*G*Mearth/rmax/ sqrt((2.*2+1.)*fac(2-2)/fac(2+2)/4.*pi) /2.0
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
      REAL :: axis(3), eigen(3), eigev(3,3), theta, zunit(3)=(/0.,0.,1./)    
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
            !All mass anomalies
            Iload = IfromC(-77.e-6,-56.e-6,-16.5e-6,-15.6e-6,0.3e-6)            
            !Fossil buldge
            Ifoss = IfromC(-126.2e-6,56.e-6,38.9e-6,15.6e-6,-0.3e-6)

            call print_eige(Ifoss,'FOSSILBEFORE',.true.)
            call print_eige(Iload,'LOADBEFORE',.true.)
            call print_eige(Iload+Ifoss,'TOTALBEFORE',.true.)
            call eige_problem(Ifoss, eigen, eigev, 1)
            indice = MAXLOC((eigen), DIM=1)              
            theta = acos(dot_product(eigev(:,indice),zunit)) /pi*180.
            print *,'Rotating the fossil bulge and impact basins by ', theta,' degrees'
            axis = cross_product( eigev(:,indice), zunit )
            axis = axis / absv(axis)
            open(94,file='run/principal.dat',Access='append')
                write(94,'(7e25.15,a,a)') 0.,0.,0.,axis(1),axis(2),axis(3),theta, '  AXIS'
            close(94)
            call rotIbyVec(axis, theta, Ifoss)
            call rotIbyVec(axis, theta, Iload)            
            call print_eige(Ifoss,'FOSSILROT',.true.)
            call print_eige(Iload,'LOADROT',.true.)
            call print_eige(Iload+Ifoss,'TOTALROT',.true.)

            !Ifoss = IfromC(-156.1e-6,0.,38.8e-6,0.,0.)
            !call print_eige(Ifoss,'DIAGONAL',.false.)

            Iload = (Iload + Ifoss)*SPAf
            !stop 'Iload called'
        CASE(3)
            print *,'Icetest loading will be used, but isostatic compensation is disabled (zero surface traction form load)'
            print *,'Iload is not hard-wired: it is a surface topography that generates self-gravity, check also ice00'
            IF(.not.icetest) print *,'WARNING: No loading'
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
        
        call eige_problem(Irot,eigen,eigev,1)
        diag=0.        
        diag(1,1) = eigen(1)
        diag(2,2) = eigen(2)
        diag(3,3) = eigen(3)        
        eigev(:,1)=matmul(Rot,eigev(:,1))
        eigev(:,2)=matmul(Rot,eigev(:,2))
        eigev(:,3)=matmul(Rot,eigev(:,3))

        call inverseI(eigev, eigevinv)
        Irot=matmul(matmul(eigev, diag), eigevinv)
        !print *,'Inverse: '
        !  print *,eigevinv(:,1)
        !  print *,eigevinv(:,2)
        !  print *,eigevinv(:,3)
        !print *,'Direct: '
        !  print *,eigev(:,1)
        !  print *,eigev(:,2)
        !  print *,eigev(:,3)
        !print *,'matmul(eigev, diag)', matmul(eigev, diag)
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
        print *,'    LIOUSHELL alive      '
        print *,'<------------/-----------'        
        print *,'<--------xxx/------------'
        print *,'<------xoooVox-----------'
        print *,'<----xooo~/~ooox-----_---'
        print *,'<------xo/ooox------/ \--'
        print *,'<-------/xxx--------\o/--'
        print *,'<------/-------------^---'
        print *
      END SUBROUTINE LIve

      END PROGRAM modul_verze
