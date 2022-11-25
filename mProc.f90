      MODULE mProc
      USE mShared
      USE mConst
      USE nr
      IMPLICIT NONE
      
    CONTAINS
    
      SUBROUTINE derivs(t,omegain,dydx)
      REAL :: t, omegain(3)
      REAL, INTENT(OUT) :: dydx(3)
      
        call fcn(t,omegain,dydx) 
      
      END SUBROUTINE
      
      SUBROUTINE fcn(t,omegain,dydx)
      IMPLICIT NONE
      REAL :: t, omegain(3)       !INTENT not stated because of wMIA approximation possibility
      REAL, INTENT(OUT) :: dydx(3)
      REAL :: Itp(3,3),rhs(3),ewJw(3),Jw(3),dd,dt,dJdtw(3),omega(3),Imoment(3),Iaxis(3,3)
      REAL :: dmkl(3),tau(3),th3,th4,th5,th6,th7,rat(0:jmax),tjpo,Qbulge(3,3),loadbg(3)
      REAL, SAVE :: Inertiaprev(3,3,2), dtfakt=-1.
      COMPLEX, DIMENSION(0:jmax) :: dVget,dVgit,Vextload,rotidpot
      COMPLEX ::y(NN,0:jmax), Mem(3*N,0:jmax), dVge(2:N,0:jmax), dVgi(2:N,0:jmax), ko(Nvinp,0:jmax)
      INTEGER :: indxx(3),k1,k2,k3,iterace,m,qfi,qfitot,ip
      LOGICAL :: wesave,timetodraw,apply_wMIA
      
      dt=t-tlast;         
      IF(t>=rozjezd) THEN; nrhse=nrhse+1; IF(nrhse>1) trialdtmax=max(dt,trialdtmax); ENDIF;
      IF ( (memschema==1 .or. memschema==3).and.(dt/=dtfakt) ) THEN
          call faktorizace(dt,0,jmax);
          dtfakt=dt
      ENDIF    
      call cpu_time(th3)
      wesave = (t==trel).and.(savestate==1).and.(imprnt==1)       
      
      IF(selfg==1 .or. memschema==2 .or. memschema==4) THEN
        iterace = sgit + merge(30,0,dt==0.)
        IF(sgit<1) stop 'sgit must be non-zero: selfgravity term requires iterations'
        IF(memschema==4 .and. sgit<2) stop 'sgit must be > 2: Simpsons rule requires iterations'
      ELSE
        iterace = 0
      ENDIF
      IF (sgres) THEN
        selfg1=0.; selfg2=0.; Vpres=0.; 
      ENDIF
      ! Apply_wMIA changes omegain, no additional ODE solver that would automatically advance yode should thus be used
      apply_wMIA = ((pr==1).or.(couple==1)).and.(wMIA>0).and.(imprnt==1)
      ! Keep omegain constant for the first step in order to determine absH0 that is conserved in the next steps
      IF(pr==0.and.t==0.0) apply_wMIA = .false.
      qfitot = 1
      IF(apply_wMIA) qfitot = qfit * merge(5,1,dt==0.0)
      ! Slow-rotator case: you are not iteratively converging to w_(n+1), but incrementing it per each timestep (dt involved)
      IF(wMIA==4 .AND. Hu_mode==1) qfitot = 1
      omega = omegain

      DO qfi=1,qfitot
          j = jmax    
          tjpo = 2.*j + 1.
!$OMP PARALLEL DO
          DO m=0,j  ! We iterate for each m=0,1,2 separately
            DO k3=1,1+iterace                
                IF(selfg==0) THEN; selfg1(:,m)=0.; selfg2(:,m)=0.; Vpres(m)=0.; ENDIF;
                rotidpot(m) = phi(j,m,omega,'rot') + phi(j,m,tidax,'tid')
                IF(selfg==1.and.j/=2) stop 'selfgravity calculations valid only for j==2'
                DO k1=1,3
                    ! NOTE: The original equation is divided by the shear modulus 
                    SELECT CASE(memschema) 
                    CASE(0)
                        Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*ylast((3+k1):NN:6,m)/visko
                    CASE(1)
                        Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*ylast((3+k1):NN:6,m)/visko/2.
                    CASE(2)                   
                        IF(k3==1) THEN
                            Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*ylast((3+k1):NN:6,m)/visko
                        ELSE
                            Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*(ylast((3+k1):NN:6,m)+y((3+k1):NN:6,m))/visko/2. 
                        ENDIF 
                    CASE(3)
                        Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m)
                    CASE(4)
                        IF(k3==1) THEN
                            Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt/2.*ylast((3+k1):NN:6,m)/visko
                        ELSEIF(k3 <  1+iterace/2) THEN
                            Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt/2.*y((3+k1):NN:6,m)/visko                    
                        ELSE                                                
                            Mem(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + &
                            dt/6.*(ylast((3+k1):NN:6,m) + 4.*ymidp((3+k1):NN:6,m) + y((3+k1):NN:6,m))/visko 
                        ENDIF            
                    END SELECT                
                ENDDO  
                y(:,m) = 0.; 
                ! Surface loading in case isotest or icetest
                y(1,m) = -sqrt(j/tjpo)*rho0(1)*g0(1)*extload(m)*sclBC                        
                y(2,m) = sqrt((j+1)/tjpo)*rho0(1)*g0(1)*extload(m)*sclBC
                ! Neglecting surface traction of the load for fixload==3
                IF(fixload==3) y(:,m) = 0.
                y(4,m) = Mem(1,m)*sclREO; 
                y(5,m) = Mem(2,m)*sclREO; 
                y(6,m) = Mem(3,m)*sclREO;
                DO k1=2,N               
                    ! centrifugal force 
                    IF(j==2) y(6*k1-5,m) = -5*sqrt(2./5)*rcent(k1)*rotidpot(m)*0.5*(rho0r(k1)+rho0r(k1-1))                                
                    ! selfgravitation
                    y(6*k1-5,m) = y(6*k1-5,m) + selfg1(k1,m); 
                    y(6*k1-4,m) = y(6*k1-4,m) + selfg2(k1,m);
                    ! the zero degree component of centrifugal force
                    IF(j==0) y(6*k1-4,m) = 2*rcent(k1)*rotidpot(m)*rho0(k1)
                    ! memory term in the rheological equation
                    y(6*k1-2,m) = Mem(3*k1-2,m)*sclREO; 
                    y(6*k1-1,m) = Mem(3*k1-1,m)*sclREO; 
                    y(6*k1,m) = Mem(3*k1,m)*sclREO;
                ENDDO             
                ! Term drho*varphi (neglected in the standard formulation)
                IF (model==5.and.j==2.and.k3>1.and.extended) THEN
                    DO k1=1,Nm5l-2
                        ip = sjmp(k1) 
                        y(6*ip-5,m) = y(6*ip-5,m)-5.*sqrt(2./5.)*rcent(ip)*rotidpot(m)*(rho_inp(k1+1)-rho_inp(k1))*ur(ip,m)/dr
                    ENDDO
                ENDIF
                ! Core pressure
                y(6*N+1,m) = -sqrt(j/tjpo)*(rho0(N+1)+rho0(N))*((redge(N)**2)*rotidpot(m)+Vpres(m))*sclBC
                y(6*N+2,m) = sqrt((j+1)/tjpo)*(rho0(N+1)+rho0(N))*((redge(N)**2)*rotidpot(m)+Vpres(m))*sclBC  
            
                IF(check_solver==1) ycheck = y(:,m)      
                yhelp(:,1,m) = real(y(:,m)) 
                yhelp(:,2,m) = imag(y(:,m)) 
                
                IF(m==0) call cpu_time(th6)     ! assuming a parallel run, measuring time for one core
                ! SOLVING THE BIG LINEAR SYSTEM AA.x = y
                SELECT CASE (ilib)
                    CASE(1)
                        ! The result will be stored in vector y
                        call lubksb(AA,6*N+2,6*N+2,indx,yhelp(:,1,m))     
                        call lubksb(AA,6*N+2,6*N+2,indx,yhelp(:,2,m))
                    CASE(2)  
                        stop 'disabled: MKL needed'
                        !call dgbtrs('N',NN,kl,ku,2,ab,ldab,ipiv,yhelp(:,:,m),NN,info);   
                        !IF (info/=0) print *,info,j,m ;     
                    CASE(3)
                        call dsvbksb(AA,w,v,NN,NN,NN,NN,yhelp(:,1,m),rx)
                        call dsvbksb(AA,w,v,NN,NN,NN,NN,yhelp(:,2,m),ix)
                        yhelp(:,1,m)=rx; yhelp(:,2,m)=ix; 
                    CASE(4)
                        call banbks(ab,NN,kl,ku,NN,kl+ku+1,al,kl,indx,yhelp(:,1,m))
                        call banbks(ab,NN,kl,ku,NN,kl+ku+1,al,kl,indx,yhelp(:,2,m))
                END SELECT    
                y(:,m)=cmplx(yhelp(:,1,m),yhelp(:,2,m)); 
                IF(m==0) call cpu_time(th7)     ! assuming a parallel run, measuring time for one core
                
                IF (check_solver==1 .and. m==0 .and. t==trel .and. k3==1) THEN;
                    call test_presnosti(y(:,m))   
                    call zpetna_kontrola(y(:,m))      
                    call printAA(y(:,m))
                ENDIF;
                ! ur is defined at centres of layers, ur2 and ut2 are defined at layer edges
                ur(:,m) = sqrt(j/tjpo)*y(1::6,m) - sqrt((j+1)/tjpo)*y(2::6,m)           
                ur2(:,m) = 0.5*( ur(1:N,m) + ur(2:N+1,m) )
                ut2(:,m) = sqrt((j+1)/tjpo)*(y(1:6*N-5:6,m)+y(7::6,m))/2. + sqrt(j/tjpo)*(y(2:6*N-4:6,m)+y(8::6,m))/2. 
                IF(memschema==4.and.k3==iterace/2) ymidp(:,m)=y(:,m);
                
                ! SELFGRAVITY
                dVgi(:,m) = 0.; dVge(:,m) = 0.;    
                dVget(m) = 4.*pi*Gconst/tjpo*rho0(N+1)*(redge(N)**(j+2))*ur2(N,m)
                dVgit(m) = 4.*pi*Gconst/tjpo*rho0(1)*(redge(1)**(1-j))*(ur2(1,m)+extload(m))
            
                DO k1=2,N                
                    IF(model==5.and.layered) THEN
                        DO k2 = 1,Nm5l-2
                            ko(k2,m) = 4.*pi*Gconst/tjpo*(rho_inp(k2+1)-rho_inp(k2))*ur(sjmp(k2),m)
                        ENDDO
                        DO k2 = 2,Nm5l-1
                            IF (rcent(k1)>rcent(sjmp(k2-1))) THEN
                                dVge(k1,m) = dVge(k1,m) + ko(k2-1,m)*rcent(sjmp(k2-1))**(j+2)
                            ELSEIF (rcent(k1)>rmin) THEN
                                dVgi(k1,m) = dVgi(k1,m) + ko(k2-1,m)*rcent(sjmp(k2-1))**(1-j)
                            ENDIF
                        ENDDO

                    ! Model 1 is homogeneous -> density jumps are only at the surface and at CMB
                    ELSEIF (model/=1) THEN 
                      DO k2 = k1,N-1
                        rat(m) = rcent(k2+1)/rcent(k2)                    
                        dVge(k1,m)=dVge(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/(dr*(j+3.))*(rcent(k2)**(j+3.)*(1.-rat(m)**(j+3.)))&
                        &-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr*(j+4.))*(rcent(k2)**(j+4.)*(1.-rat(m)**(j+4.))))*4.*pi*Gconst/tjpo                    
                        !computationally faster option using Taylor expansion (could be very accurate if rat**(j+3) differs little from 1)
                        !dVge(k1,m)=dVge(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/(dr*(j+3))*(rcent(k2)**(j+3)*(j+3)*(1.-rat))&
                        !&-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr*(j+4))*(rcent(k2)**(j+4)*(j+4)*(1.-rat)))*4*pi*Gconst/tjpo
                        !print *,rcent(k2)**(j+4)*(1.-rat**(j+4)),rcent(k2)**(j+4)*(j+4)*(1.-rat),(rcent(k2)**(j+4)-rcent(k2+1)**(j+4))
                      ENDDO
                      DO k2 = 2,k1-1            ! Shape for j=2
                        dVgi(k1,m) = dVgi(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/dr*log(rcent(k2)/rcent(k2+1))&
                        &-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr)*(rcent(k2)-rcent(k2+1)))*4.*pi*Gconst/tjpo
                      ENDDO
                        dVgi(k1,m) = dVgi(k1,m)+((ur(1,m)+(ur(2,m)-ur(1,m))*rcent(1)/dr)*2*(rho0(2)-rho0(1))/dr*log(redge(1)/rcent(2))&
                        &-(ur(2,m)-ur(1,m))*2*(rho0(2)-rho0(1))/(dr*dr)*(redge(1)-rcent(2)))*4.*pi*Gconst/tjpo
                    ENDIF
                    
                    PhiIn(:,m)=dVgi(:,m)*rcent(2:N)**j
                    PhiEx(:,m)=dVge(:,m)/rcent(2:N)**(j+1)
                    selfg1(k1,m)=-0.5*(rho0r(k1)+rho0r(k1-1))*sqrt(tjpo*j)*rcent(k1)**(j-1)*(dVgi(k1,m)+Dvgit(m))
                    selfg2(k1,m)=-0.5*(rho0r(k1)+rho0r(k1-1))*sqrt(tjpo*(j+1))/(rcent(k1)**(2+j))*(dVge(k1,m)+Dvget(m))
                ENDDO
                ! Term drho*Phi (neglected in the standard formulation)
                IF (model==5.and.j==2.and.k3>1.and.extended) THEN
                    DO k2 = 1,Nm5l-2
                        selfg1(sjmp(k2),m) = selfg1(sjmp(k2),m) - &
                         (rho_inp(k2+1)-rho_inp(k2))/dr*ur(sjmp(k2),m)*sqrt(tjpo*j)*rcent(sjmp(k2))**(j-1)*(dVgi(sjmp(k2),m)+Dvgit(m))
                        selfg2(sjmp(k2),m) = selfg2(sjmp(k2),m) - &
                         (rho_inp(k2+1)-rho_inp(k2))/dr*ur(sjmp(k2),m)*sqrt(tjpo*(j+1))/(rcent(sjmp(k2))**(2+j))*(dVge(sjmp(k2),m)+Dvget(m))                        
                    ENDDO
                ENDIF
                
                Vpres(m) = dVget(m)/redge(N)**(j+1) + (dVgi(N,m)+dVgit(m))*redge(N)**j
                ! Self-gravity at the outer surface
                Vsurf(m) = (dVget(m)+dVge(2,m))/redge(1)**(j+1) + dVgit(m)*redge(1)**j
                selfg1Tr(m,2) = rho0(N)*sqrt(tjpo*j) * (redge(N)**(j-1.)) * (dVgi(N,m)+Dvgit(m));
                selfg2Tr(m,2) = rho0(N)*sqrt(tjpo*(j+1)) / (redge(N)**(2.+j)) * Dvget(m);
                selfg1Tr(m,1) = rho0(1)*sqrt(tjpo*j) * (redge(1)**(j-1.)) * Dvgit(m); 
                selfg2Tr(m,1) = rho0(1)*sqrt(tjpo*(j+1)) / (redge(1)**(2.+j)) * (dVge(2,m)+Dvget(m));
            
                IF(m==0) taasolve=taasolve+th7-th6      ! assuming a parallel run, measuring time for one core
            ENDDO   ! end of selfg iterations                
          ENDDO     ! end of cycle over m=0,1,2
!$OMP END PARALLEL DO
                 
            ! INERTIA TENSOR
            Itp = Ball + Inertia_Tensor(Vsurf)            
            IF (pr==1) THEN
                IF ((t<1.5*tgrowth).and.sinsmooth) THEN
                    Itp = Itp + Iload * merge(1., sin(t/tgrowth*pi/2.), t>=tgrowth)
                ELSEIF (t<1.5*tgrowth) THEN
                    Itp = Itp + Iload * merge(1., t/tgrowth, t>=tgrowth)
                ELSE
                    Itp = Itp + Iload
                ENDIF
                IF(fosslit.and.add_Ifb) Itp = Itp + (J0nolit-J0)*sclJ
            ENDIF
            Itp = Itp / sclJ        

            ! Equilibrium position balances inertia contribution of the load with unrelaxed part of rotational bulge
            IF (apply_wMIA) THEN
                omegain = new_omg(Itp, Inertiaprev, omega, dt, t, qfi)
                IF ((absv(omegain-omega)/absv(omega) < eps0).and.(qfit_flexible)) THEN
                    omega = omegain
                    meanqfit(1) = meanqfit(1) + qfi
                    meanqfit(2) = meanqfit(2) + 1.0
                    exit;
                ELSEIF (check_solver==1) THEN
                    print *,'wMIA |w_(n+1)|/|w_n|: ', absv(omegain-omega)/absv(omega), angledeg(omegain,omega)
                ENDIF
                omega = omegain
            ENDIF
      ENDDO                     ! end of quasi-fluid iterations
      qfit_used = min(qfi,qfitot)

        IF (imprnt==1) THEN
            DO m=0,j
                DO k1=1,3;
                    SELECT CASE(memschema)
                    CASE(0)
                        Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*ylast((3+k1):NN:6,m)/visko
                    CASE(1,2)
                        Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*(ylast((3+k1):NN:6,m) + y((3+k1):NN:6,m))/visko/2.
                    CASE(3)
                        Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*y((3+k1):NN:6,m)/visko
                    CASE(4)                        
                        Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + &
                        dt/6.*(ylast((3+k1):NN:6,m) + 4.*ymidp((3+k1):NN:6,m) + y((3+k1):NN:6,m))/visko 
                    END SELECT                   
                ENDDO;
                ylast(:,m) = y(:,m) 
                tlast = t 
                timetodraw = ((t>=trel*2.**real(npicsav-NpicSavMax+1).or.t==0.).and.pr==0)
                IF(loadstate==1) timetodraw = (t>=npicsav*tmax/NpicSavMax)
                IF(modd/=0.and.timetodraw) call rez_radialni_slozky_v(y(:,m),Ntheta,Nphi,N,j,m,modd)    
                Vextload(m) = 4.*pi*Gconst/tjpo*rho0(1)*rmax*extload(m)
                IF(.not.isotest) k2t(m) = real(Vsurf(m)/(rotidpot(m)*redge(1)**2.)) 
                IF(isotest) k2t(m) = real((Vsurf(m)-Vextload(m)) / Vextload(m))
            ENDDO
        ENDIF

        ! Inertia array is shared among subroutines outside of derivs and recorded for the evaluation of dJdt
        Inertia = Itp
        IF (imprnt==1) THEN
            Erot = 0.5*dot_product( matmul(Itp, omega), omega ) *sclJ * sclOmega**2.
            Etid = 0.5*dot_product( matmul(Itp, tidax), tidax ) *sclJ * eas**2. ! - 2.*pi*Isph*dot_product(tidax,tidax)
            call dissipation(omega,tidax,dt,t)
        ENDIF

        IF (imprnt==1 .and. modd/=0 .and. timetodraw ) THEN
          call zapis3D(N,Ntheta,Nphi,npicsav); npicsav=npicsav+1;
        ENDIF;
        call cpu_time(th4)

        ! Vsurf and phi in the code both have unconvetional signs, but in the Love numbers (ratios) it cancels out
        IF (isotest) THEN
            IF(pr==0.and.t==0.) print *,'loading k2 at time=0 (elastic limit): ', k2t(0)
            IF(pr==0.and.t==trel) print *,'loading k2 at t=trel (fluid limit): ', k2t(0)
        ELSE
            IF(pr==0.and.t==trel) THEN; print *,'tidal k2 at t=trel (fluid limit): ',k2t(0); k2Tf=k2t(0); ENDIF;
            IF(pr==0.and.t==0.) THEN; print *,'tidal k2 at time=0 (elastic limit): ',k2t(0); k2Te=k2t(0); ENDIF;
        ENDIF
        
        Jw = matmul(Itp,omega)
        IF (tcross==1) THEN
            call eige_problem(Itp,Imoment,Iaxis)
            ewJw = matmul(Imoment*matmul(omega,Iaxis),tcross_product(omega,Iaxis))*sclOmega
        ELSE
            ewJw = cross_product(omega,Jw)*sclOmega        
        ENDIF
        
        dJdtw=matmul(dJdt(Itp,Inertiaprev,dt,t),omega)

        IF (dt/=0.) THEN                
            rhs=-(dJdtw+ewJw)      !multiplied by -1 in order to use directly Itp in the solver
            IF(imprnt==1) dJdtw0=dJdtw
        ELSE
            rhs=-(dJdtw0+ewJw)     !initial condition dJdtw0 is either zero or loaded from inistate.dat
        ENDIF
        IF (imprnt==1) THEN
            Inertiaprev(:,:,1) = Inertiaprev(:,:,2)
            Inertiaprev(:,:,2) = Inertia
        ENDIF
           
        call ludcmp(Itp,3,3,indxx,dd)
        IF (crash==1) return
        call lubksb(Itp,3,3,indxx,rhs)
        dydx=rhs         
        
        call cpu_time(th5)
        tnapln = tnapln + th4 - th3;
        tliou = tliou + th5 - th4;
        IF(t==trel.and.pr==0) nstp=1;
      
      END SUBROUTINE fcn

      FUNCTION Inertia_Tensor(Vpot)
      COMPLEX, INTENT(IN) :: Vpot(0:jmax)
      REAL :: Inertia_Tensor(3,3), cst
        cst=(5.*rmax**3.)/(2*pi*Gconst)
        Inertia_Tensor(1,1) = (sqrt(pi/5)*cst/3*real(Vpot(0))-cst*sqrt(2*pi/15)*real(Vpot(2)))
        Inertia_Tensor(2,2) = (sqrt(pi/5)*cst/3*real(Vpot(0))+cst*sqrt(2*pi/15)*real(Vpot(2)))
        Inertia_Tensor(3,3) = (-2*sqrt(pi/5)*cst/3*real(Vpot(0)))
        Inertia_Tensor(1,3) = cst*sqrt(2*pi/15)*real(Vpot(1))
        Inertia_Tensor(2,3) = -cst*sqrt(2*pi/15)*imag(Vpot(1))
        Inertia_Tensor(1,2) = cst*sqrt(2*pi/15)*imag(Vpot(2))  
        Inertia_Tensor(2,1) = Inertia_Tensor(1,2);
        Inertia_Tensor(3,1) = Inertia_Tensor(1,3);
        Inertia_Tensor(3,2) = Inertia_Tensor(2,3);
      END FUNCTION Inertia_Tensor

      SUBROUTINE faktorizace(dt,inicializace,j)
      REAL :: dt
      INTEGER :: inicializace,j,i1,i2,k1,k2
        ! CHECK THIS SUBROUTINE OVER NOTES: Why is there dt/visko in A(i,i)? For initialization dt=0, and when dt/=0 that part is never called!
        call cpu_time(th1)
        IF (inicializace==1) THEN
            call napln_A(j,A)       ! Thanks to dividing rheo Eq. by viscosity, only the diagonal elements of A are r-dependent.
            call napln_B(j,B)       ! The diagonal elements are filled only later in only to avoid too frequent matrix assembly.
            
            IF(ilib/=2.or.opt/=0) THEN
                AA=0.
                IF(bot_noslip) THEN
                    call BC_quasifreesurf(j,'outer'); call BC_noslip('inner');
                ELSE
                    call BC_quasifreesurf(j,'outer'); call BC_quasifreesurf(j,'inner');
                ENDIF
            ELSE
                abini=0.
                IF(bot_noslip) THEN
                    call BCab_quasifreesurf(j,'outer'); call BCab_noslip('inner');
                ELSE
                    call BCab_quasifreesurf(j,'outer'); call BCab_quasifreesurf(j,'inner');
                ENDIF
            ENDIF
            DO k1=1,N       
                DO k2=4,6   ! Completing the diagonal elements of A + including the shear modulus, by which it is divided in reo
                    IF (memschema==1) THEN      
                        ! Should be -r/shear, but in rheology (discretization in r) we omit the division by r
                        A(k2,k2)=-sclREO*(1./shear(k1) + dt/visko(k1)/2.)        
                    ELSEIF (memschema==3) THEN
                        A(k2,k2)=-sclREO*(1./shear(k1) + dt/visko(k1))        
                    ELSE                    
                        A(k2,k2)=-sclREO/shear(k1)
                    ENDIF    
                ENDDO                   
                IF(ilib/=2.or.opt/=0) THEN
                    ! In this order we fill AA: first BC, then R.K, then Reo, then P.R.
                    call rovkon(k1)     
                    call reovztah(k1)
                ELSE
                    call rovkonab(k1)
                    call reovztahab(k1)
                ENDIF
                
                IF (k1/=1) THEN
                    ! Adding -(u.grad (rho0))g0: again we neglect dividing A by r (it is done later in EOM) 
                    A(1,1)=-g0(k1)*j/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr
                    A(1,2)=g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr
                    A(2,1)=g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr
                    A(2,2)=-g0(k1)*(j+1)/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr
                    IF(ilib/=2.or.opt/=0) call pohrov(k1-1)     
                    IF(ilib==2.and.opt==0) call pohrovab(k1-1)      
                ENDIF   
            ENDDO
            SELECT CASE(ilib)
                CASE(2)
                    if(opt/=0) then
                        forall (i1=1:NN,i2=1:NN,((max(1,i2-ku) <= i1 ).and.(i1 <= min(NN,i2+kl)))) ; abini(kl+ku+1+i1-i2,i2)=AA(i1,i2); 
                    end forall; endif;    
                CASE(4)
                    forall (i2=1:NN,i1=1:(kl+ku+1),((i2 >= kl+2-i1 ).and.(i2 <= kl+NN+1-i1 ))) ; abini(i2,i1)=AA(i2,i2-kl+i1-1); end forall
            END SELECT
        ENDIF
           
        IF(check_solver==1) AAcheck=AA
        SELECT CASE (ilib)
            CASE(1)                    
                call ludcmp(AA,6*N+2,6*N+2,indx,d)              !solution of the "big" linear problem AAx=y
            CASE(2)
                stop 'disabled: MKL needed'
                !forall (i1=1:N,i2=0:2) ; abini(kl+ku+1,6*i1-i2)=-sclREO/shear(i1)*(1+shear(i1)*dt/visko(i1)/2); END forall
                !ab=abini;
                !call dgbtrf(NN,NN,kl,ku,ab,ldab,ipiv,info);
            CASE(3)
                call dsvdcmp(AA,NN,NN,NN,NN,w,v)
                wmax=0.             !Will be the maximum singular value obtained.
                DO k1=1,NN
                    IF(w(k1).gt.wmax)wmax=w(k1)
                ENDDO 
                wmin=wmax*1.0e-12   !This is where we set the threshold for singular values allowed to be nonzero.
                DO k1=1,NN
                    IF (w(k1).lt.wmin) THEN; w(k1)=0. ; ENDIF;
                ENDDO 
            CASE(4)      
                forall (i1=1:N,i2=0:2) ; abini(6*i1-i2,kl+1)=-sclREO/shear(i1)*(1+shear(i1)*dt/visko(i1)/2); END forall
                ab=abini;
                call bandec(ab,NN,kl,ku,NN,kl+ku+1,al,kl,indx,nrd)
        END SELECT    
        call cpu_time(th2)
        tfakt=tfakt+th2-th1
      END SUBROUTINE faktorizace
            
      SUBROUTINE fcnj(n,x,y,dfdy)
      INTEGER :: n,i
      REAL :: x,y(n),dfdy(n,n),f1(n),f3(n,n)
      REAL :: deltay(n,n)
        imprnt=0;
        deltay=0.;
        DO i=1,n; 
            deltay(i,i)=maxval(abs(y))/1000; 
        ENDDO;
        call fcn(x,y,f1);
        DO i=1,n; 
            call fcn(x,y+deltay(:,i),f3(:,i)); 
            dfdy(:,i)=(f3(:,i)-f1)/deltay(i,i);   
        ENDDO;      
        imprnt=1;
      END SUBROUTINE      
        
      SUBROUTINE jac(x,y,dfdx,dfdy,n)
      INTEGER :: n,i,tmpimp
      REAL :: x,y(n),dfdx(n),dfdy(n,n),f2(n),f1(n),f3(n,n)
      REAL :: deltax,deltay(n,n)
        deltax = yr/1000
        deltay = 0.0
        DO i=1,n
            deltay(i,i) = maxval(abs(y))/1000
        ENDDO
        tmpimp = imprnt 
        imprnt = 0
        call derivs(x,y,f1)
        call derivs(x+deltax,y,f2)
        DO i=1,n
            call derivs(x,y+deltay(:,i),f3(:,i))
            dfdy(:,i) = (f3(:,i)-f1)/deltay(i,i)
        ENDDO
        imprnt = tmpimp
        dfdx = (f2-f1)/deltax
      END SUBROUTINE

      ! Components of the rotational (or tidal) potential
      FUNCTION phi(j,m,w,mode)
      INTEGER, INTENT(IN) :: j,m
      REAL, INTENT(IN) :: w(3)
      CHARACTER(3), INTENT(IN) :: mode
      REAL :: Amp
      COMPLEX :: phi
        IF (mode=='rot') THEN
            Amp = sclOmega**2.
        ELSE ! tidal potential
            Amp = -eas**2.
        ENDIF
        IF (j==2) THEN
            SELECT CASE(m)
            CASE(0)
                phi = sqrt(pi/5.) * (w(1)*w(1)+w(2)*w(2)-2*w(3)*w(3)) * Amp/3.
            CASE(1)
                phi = sqrt(6.*pi/5.) * (w(1)*w(3)-(0.,1.)*w(2)*w(3)) * Amp/3.
            CASE(2)
                phi = sqrt(3.*pi/10.) * (w(2)*w(2)-w(1)*w(1)+2*w(1)*w(2)*(0.,1.)) * Amp/3.
            END SELECT  
        ELSEIF (j==0) THEN
            phi = 2.*sqrt(pi)*dot_product(w,w) * Amp/3.
        ENDIF
      END FUNCTION phi
      
      FUNCTION  tcross_product(a,b)
      INTEGER :: i
      REAL,INTENT(IN) :: a(3),b(3,3)
      REAL :: tcross_product(3,3)
        DO i=1,3
            tcross_product(i,1)=a(2)*b(3,i)-a(3)*b(2,i)
            tcross_product(i,2)=a(3)*b(1,i)-a(1)*b(3,i)
            tcross_product(i,3)=a(1)*b(2,i)-a(2)*b(1,i)      
        ENDDO 
      END FUNCTION
      
      FUNCTION  cross_product(a,b)
      INTEGER :: i
      REAL,INTENT(IN) :: a(3),b(3)
      REAL :: cross_product(3)
        cross_product(1)=a(2)*b(3)-a(3)*b(2)
        cross_product(2)=a(3)*b(1)-a(1)*b(3)
        cross_product(3)=a(1)*b(2)-a(2)*b(1)      
      END FUNCTION
      
      FUNCTION EnG(mode)
      REAL :: EnG,ura,HomoEnG
      INTEGER :: NR=130,Ntheta=180,Nphi=360,k1,k2,k3,m,mode,faktor
      REAL :: phi,theta,r,GE,r_topo,r_topo_CMB,r_low,r_top,drho,dz,dV,dS,dS_CMB
      REAL :: GEcontrol,Vcontrol,Vcontrol2,GEcontrol2,mvalec,mvalec_CMB,chyba,Vsg,summ1,summ2
      REAL,SAVE :: GE0=0.

      IF (mode==-1 .or. mode==0) THEN
          GE=0.; chyba=0.; SEng=0.; Vsg=0.; Ssg=0.; SEng_CMB=0.; Ssg_CMB=0.; EdrhoV=0.;
          Vcontrol=0.; GEcontrol=0.; GEcontrol2=0.; Vcontrol2=0.;
          DO k2=0,Ntheta-1              
            theta=k2*pi/Ntheta            
            DO k3=0,Nphi-1 
                phi=2.*k3*pi/Nphi                
                !Summing over the harmonic components m=-1, m=-2 is done in function SphHarm
                summ1=0.; summ2=0.;
                DO m=0,jmax
                    IF(order_analyse .and. m/=mp) cycle;
                    summ1 = summ1 + real( (ur2(1,m)+extload(m)) * SphHarm(theta,phi,j,m) )
                    summ2 = summ2 + real( ur2(N,m) * SphHarm(theta,phi,j,m) )
                ENDDO     
                r_topo=rmax + summ1
                r_topo_CMB=rmin + summ2
                r_low=min(rmax,r_topo); r_top=max(rmax,r_topo)
                drho=sign(rho0(1),r_topo-rmax)                
                dz=(r_top-r_low)/NR                               
                DO k1=0,(NR-1)*merge(1,-1,mode==0)
                    !pouze horni cast hranice - pocital jsem pro srovnani s plosnymi integraly                    
                    r= r_low+k1*dz
                    dV= abs(r**2*sin(theta)*dz*(2*pi/Nphi)*(pi/Ntheta))
                    GE= GE - drho*G_pot(r)*dV 
                    GEcontrol= GEcontrol + drho*dV
                    GEcontrol2= GEcontrol2 + drho*g0(1)*(r-rmax)*dV
                    chyba= chyba + drho*(-G_pot(r)-g0(1)*(r-rmax))*dV
                    ! here the radial dependence of selfg is neglected
                    Vsg= Vsg + dV*drho *real( Vsurf(0) * SphHarm(theta,phi,j,0) ) 
                ENDDO                
                mvalec=g0(1)*rho0(1)*(r_topo-rmax)**2./2.
                mvalec_CMB=g0(N+1)*rho0(N+1)*(r_topo_CMB-rmin)**2./2.
                dS=rmax**2.*sin(theta)*(2.*pi/Nphi)*(pi/Ntheta)
                dS_CMB=rmin**2.*sin(theta)*(2.*pi/Nphi)*(pi/Ntheta)
                Vcontrol = Vcontrol + dS*abs((r_topo-rmax))
                Vcontrol2 = Vcontrol2 + dS*(r_topo-rmax)
                SEng = SEng + dS*mvalec
                SEng_CMB = SEnG_CMB + dS_CMB*mvalec_CMB
                summ1=0.; summ2=0.;
                DO m=0,jmax
                    IF(order_analyse .and. m/=mp) cycle;
                    summ1 = summ1 + real( Vsurf(m) * SphHarm(theta,phi,j,m) )
                    summ2 = summ2 + real( Vpres(m) * SphHarm(theta,phi,j,m) )
                ENDDO                
                Ssg= Ssg + dS*rho0(1)*summ1*(r_topo-rmax)
                Ssg_CMB= Ssg_CMB + dS_CMB*rho0(N+1)*summ2*(r_topo_CMB-rmin)                
            ENDDO
          ENDDO
          
          DO k2=0,Ntheta-1              
            theta=k2*pi/Ntheta            
            DO k3=0,Nphi-1 
                phi=2.*k3*pi/Nphi
                DO k1=2,N                    
                    summ1=0.;
                    DO m=0,jmax
                        summ1 = summ1 + real( ur(k1,m) * SphHarm(theta,phi,j,m) )
                    ENDDO     
                    r_topo=rcent(k1) + summ1
                    dS=rcent(k1)**2.*sin(theta)*(2.*pi/Nphi)*(pi/Ntheta)
                    mvalec=g0(k1)*(rho0r(k1)-rho0r(k1-1))*(r_topo-rcent(k1))**2./2.
                    EdrhoV=EdrhoV + dS*mvalec
                ENDDO                
            ENDDO
          ENDDO
          
          IF(mode==0) THEN
            print *,'GE, chyba proti GE pocitano objemove s konst. g0, chyba proti pouzitemu povrchovemu vypoctu'
            print *,GE, chyba, GE-SEng
            print *,'GE pomoci konstantniho g0, GE pocitano povrchove'
            print *,GEcontrol2,SEng
            print *,'kontrola zmeny celkoveho objemu telesa pri povrchovem vypoctu'
            print *,Vcontrol2/Vcontrol
            print *,'Chyba celkoveho objemu pri objemovem vypoctu, G_pot(rmax)'
            print *,GEcontrol,G_pot(rmax)
            print *,'Srovnani objemoveho a povrchoveho vypoctu Ssg (no core, m=0)'
            print *,Vsg,Ssg
            stop  
          ENDIF
          
          Ssg = Ssg*selfg;  Ssg_CMB = Ssg_CMB*selfg;  EdrhoPhi = EdrhoPhi*selfg;        
          ! without selfgravity only the energy of mass in reference potential should be computed
          ! ve smyslu rovnice z clanku se jedna o (v tomto poradi): [rho]V, (drho/dr)V, [rho]Phi, (drho/dr)Phi 
          EnG = (SEng + SEng_CMB) + EdrhoV - (Ssg + Ssg_CMB)/2. + EdrhoPhi
          ! V kodu je znamenko sg potencialu obracene (sila = grad(\Phi))- proto zde musi byt minus
          write(33,'(20e27.15)') tlast/yr,SEng,real(ur2(1,0)*redge(1))**2*g0(1)*rho0(1)/2,real(ur2(1,:))
          
      ! EnG(1), the default mode    
      ELSE  
          summ1=0.; summ2=0.; EdrhoV=0.;
          DO m=0,jmax
            faktor=merge(1.,2.,m==0)  
            IF(order_analyse .and. m/=mp) cycle;
            IF(fixload/=3) THEN
                summ1 = summ1 + faktor*real( (ur2(1,m)+extload(m)) * conjg(ur2(1,m)+extload(m)) )*redge(1)**2*g0(1)*rho0(1)/2
            ELSE
                ! Surface traction of the load is neglected and thus also its EnG in reference g0 is neglected
                summ1 = summ1 + faktor*real( ur2(1,m) * conjg(ur2(1,m)) )*redge(1)**2*g0(1)*rho0(1)/2
            ENDIF
            summ2 = summ2 + faktor*real( ur2(N,m) * conjg(ur2(N,m)) )*redge(N)**2*g0(N)*rho0(N+1)/2
          ENDDO
          SEng= summ1; SEng_CMB= summ2;
          
          DO k1=2,N                    
            summ1=0.;
            DO m=0,jmax
                faktor=merge(1.,2.,m==0)
                summ1 = summ1 + faktor*real( ur(k1,m)*conjg(ur(k1,m)) )*rcent(k1)**2*g0(k1)*(rho0r(k1)-rho0r(k1-1))/2.
            ENDDO     
            EdrhoV=EdrhoV + summ1
          ENDDO
          
          !prvni tri cleny jsou energie vnejsi, CMB a "vnitrni" topografie v referencnim potencialu - rostou s rostouci deformaci (V0 se brani deformaci)
          !posledni clen zahrnuje energie vnejsi, CMB a "vnitrni" topografie v SG potencialu - klesa s rostouci deformaci, protoze SG pomaha deformaci
          !"vnitrni" topografie jsou dusledkem toho, ze na nehomogenni model se divame jako na mnohovrstevnaty model s infinitesimalnimi skoky v rho
          EnG = SEng + SEng_CMB + EdrhoV + Esgrav*selfg
      ENDIF      
      END FUNCTION EnG
      
      FUNCTION G_pot(r)
      REAL :: G_pot,r,Mbig,Msmall
        MBig=4.*pi*rho0(1)/3.*rmax**3.
        Msmall=4.*pi*(rhocore-rho0(1))/3.*rmin**3
        IF(r>=rmax) THEN
            G_pot=-Gconst*(Mbig+Msmall)/r    
        ELSEIF(r>=rmin) THEN
            G_pot=-Gconst*Mbig/(2.*rmax**3.)*(3.*rmax**2.-r**2.)-Gconst*Msmall/r
        ELSE
            G_pot=-Gconst*Mbig/(2.*rmax**3.)*(3.*rmax**2.-r**2.)-Gconst*Msmall/(2.*rmin**3.)*(3.*rmin**2.-r**2.)
        ENDIF         
      END FUNCTION G_pot
      
      SUBROUTINE eige_problem(J,eigen,eigev)
      REAL, INTENT (IN) :: J(3,3)
      REAL :: Itp(3,3),eigen(3),eigev(3,3)
      INTEGER :: nrot
        ! Subroutine jacobi() destroys the analyzed tensor, and so a temporary array must be used
        Itp = J
        call jacobi(Itp,3,3,eigen,eigev,nrot)
      END SUBROUTINE eige_problem
      
      SUBROUTINE dissipation(omega,tidax,dt,t)
      IMPLICIT NONE
      REAL, INTENT(IN) :: omega(3),tidax(3),dt,t
      REAL :: faktor,ymin,ymax,velicina(2:N-1),velicinas(2:N),rmax2,rmin2
      REAL, SAVE :: BoFoPo(nintp),BoRoPo(nintp),BoTiPo(nintp),SuTrPo(nintp),CMBTopoPo(nintp),BougrPo(nintp)
      REAL, SAVE :: CMBSgPrPo(nintp),CMBRoPrPo(nintp),CMBTiPrPo(nintp),DD(nintp),DDdot(nintp)
      REAL, SAVE :: rot_src(nintp),tid_src(nintp),taxprev(3),Io(3)
      COMPLEX :: v(N+1,0:jmax,2),v2(N,0:jmax,2),Ddotvek(1:N),rot_pom,vrad(N,0:jmax)
      COMPLEX :: ugr(2:N,2),PhiExt(2:N),PhiInt(2:N),trakce_rad,pcore,rotpot(0:jmax),tidpot(0:jmax)
      COMPLEX, SAVE :: ylastprev(2,NN,0:jmax)
      INTEGER :: m,q,ia,i,ip,k1
      LOGICAL, SAVE :: emplacement, virgin=.true.
      
        rmax2=rmax**2.; rmin2=rmin**2.;
        DD(nintp)=0.; DDdot(nintp)=0.; Esgrav=0.; Edef=0.; Eel=0.; EdrhoPhi=0.;
        SuTrPo(nintp)=0.; CMBTopoPo(nintp)=0.; BoFoPo(nintp)=0.; CMBSgPrPo(nintp)=0.; rot_src(nintp)=0.; tid_src(nintp)=0.
        BoRoPo(nintp)=0.; BoTiPo(nintp)=0.; CMBRoPrPo(nintp)=0.; CMBTiPrPo(nintp)=0.; BougrPo(nintp)=0.
        ! The integral quantities are not smooth over the pr=0 -> pr=1 transition (everything is restarted)
        IF(t==0.) ylastprev=(0.,0.)
        
        DO m=0,jmax
            IF(order_analyse .and. m/=mp) cycle;
            ip=1;
            rotpot(m) = phi(j,m,omega,'rot')
            tidpot(m) = phi(j,m,tidax,'tid')
            DO i=1,NN,6
                ! Using higher-order derivative when calculating velocity is needed for accurate integrals during free oscillations
                v(ip,m,1) = slope(ylast(i,m), ylastprev(:,i,m), dt, t)
                v(ip,m,2) = slope(ylast(i+1,m), ylastprev(:,i+1,m), dt, t)
                ip=ip+1    
            ENDDO            
            v2(:,m,1) = (v(1:N,m,1)+v(2:N+1,m,1))/2.
            v2(:,m,2) = (v(1:N,m,2)+v(2:N+1,m,2))/2.            
            vrad(:,m) = sqrt(j/(2.*j+1))*v2(:,m,1) - sqrt((j+1)/(2.*j+1))*v2(:,m,2)
            faktor=merge(1.,2.,m==0)
            
            ! Summing over our three components of deviatoric stress
            DO q=0,2
                ip=1;
                DO i=6-q,6*N-q,6        
                  Ddotvek(ip) = slope(ylast(i,m), ylastprev(:,i,m), dt, t)
                  ip=ip+1
                ENDDO
                ymin = real( ylast(6*N-q,m) * conjg(Ddotvek(N)) ) / (2.*shear(N))
                ymax = real( ylast(6-q,m) * conjg(Ddotvek(1)) ) / (2.*shear(1))
                velicina = real( ylast(6*2-q:6*(N-1)-q:6,m) * conjg(Ddotvek(2:N-1)) ) / (2.*shear(2:N-1))
                DDdot(nintp) = DDdot(nintp) + faktor*integrate(redge(2:N-1),rmin,rmax,ymin,ymax,velicina,2)
                
                ymin = real( ylast(6*N-q,m) * conjg(ylast(6*N-q,m)) ) / (2.*visko(N))
                ymax = real( ylast(6-q,m) * conjg(ylast(6-q,m)) ) / (2.*visko(1))
                velicina = real( ylast(6*2-q:6*(N-1)-q:6,m) * conjg(ylast(6*2-q:6*(N-1)-q:6,m)) ) / (2.*visko(2:N-1))
                DD(nintp) = DD(nintp) + faktor*integrate(redge(2:N-1),rmin,rmax,ymin,ymax,velicina,2)
                
                ymin = real( ylast(6*N-q,m) * conjg(ylast(6*N-q,m)) ) / (4.*shear(N))
                ymax = real( ylast(6-q,m) * conjg(ylast(6-q,m)) ) / (4.*shear(1))
                velicina = real( ylast(6*2-q:6*(N-1)-q:6,m) * conjg(ylast(6*2-q:6*(N-1)-q:6,m)) ) / (4.*shear(2:N-1))
                Eel = Eel + faktor*integrate(redge(2:N-1),rmin,rmax,ymin,ymax,velicina,2)
            ENDDO            
            ymin = rho0(N)*real( v2(N,m,1)*conjg(v2(N,m,1)) + v2(N,m,2)*conjg(v2(N,m,2)) )/2.
            ymax = rho0(1)*real( v2(1,m,1)*conjg(v2(1,m,1)) + v2(1,m,2)*conjg(v2(1,m,2)) )/2.
            velicina = rho0(2:N-1)*(real( v2(2:N-1,m,1) * conjg(v2(2:N-1,m,1)) + v2(2:N-1,m,2) * conjg(v2(2:N-1,m,2)) ))/2.
            IF(t/=0.) Edef = Edef + faktor*integrate(redge(2:N-1),rmin,rmax,ymin,ymax,velicina,2)
            
            ymin = -(real( selfg1Tr(m,2) * conjg(v2(N,m,1)) + selfg2Tr(m,2) * conjg(v2(N,m,2)) ))
            ymax = -(real( selfg1Tr(m,1) * conjg(v2(1,m,1)) + selfg2Tr(m,1) * conjg(v2(1,m,2)) ))
            velicinas = -(real( conjg(selfg1(:,m)) * v(2:N,m,1) + conjg(selfg2(:,m)) * v(2:N,m,2) ))
            BoFoPo(nintp) = BoFoPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,2)
            
            rot_pom = -5.*sqrt(2./5)*rotpot(m)   ! *rho0(k1)*rcent(k1) gives the force
            ymin = -real( rot_pom*rho0(N) * conjg(v2(N,m,1)) )
            ymax = -real( rot_pom*rho0(1) * conjg(v2(1,m,1)) )
            velicinas = -real( rot_pom*rho0(2:N) * conjg(v(2:N,m,1)) )
            BoRoPo(nintp) = BoRoPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,3)

            rot_pom = -5.*sqrt(2./5)*tidpot(m)
            ymin = -real( rot_pom*rho0(N) * conjg(v2(N,m,1)) )
            ymax = -real( rot_pom*rho0(1) * conjg(v2(1,m,1)) )
            velicinas = -real( rot_pom*rho0(2:N) * conjg(v(2:N,m,1)) )
            BoTiPo(nintp) = BoTiPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,3)
            
            trakce_rad = a0*(a1*ylast(3,m)+a2*ylast(4,m)+a3*ylast(5,m)) + b0*(b1*ylast(3,m)+b3*ylast(5,m)+b4*ylast(6,m))
            SuTrPo(nintp) = SuTrPo(nintp) + real( trakce_rad * conjg(vrad(1,m)) )*rmax2*faktor
            pcore = rhocore*((rotpot(m)+tidpot(m))*(rmin2) + Vpres(m))
            trakce_rad = a0*(a1*ylast(6*N-3,m)+a2*ylast(6*N-2,m)+a3*ylast(6*N-1,m))+b0*(b1*ylast(6*N-3,m)+b3*ylast(6*N-1,m)+b4*ylast(6*N,m))          
            CMBTopoPo(nintp) = CMBTopoPo(nintp) + real( conjg(-trakce_rad-pcore)*vrad(N,m) )*rmin2*faktor
            
            CMBSgPrPo(nintp) = CMBSgPrPo(nintp) + real( (rho0(N)+rho0(N+1))*conjg(Vpres(m)) * vrad(N,m) )*rmin2*faktor
            CMBRoPrPo(nintp) = CMBRoPrPo(nintp) + real( (rho0(N)+rho0(N+1))*rmin2*conjg(rotpot(m)) * vrad(N,m) )*rmin2*faktor
            CMBTiPrPo(nintp) = CMBTiPrPo(nintp) + real( (rho0(N)+rho0(N+1))*rmin2*conjg(tidpot(m)) * vrad(N,m) )*rmin2*faktor
            
            Esgrav = Esgrav - 0.5*faktor*(real( (ur2(1,m)+extload(m))*rho0(1) * conjg(Vsurf(m))*redge(1)**2. ))
            Esgrav = Esgrav - 0.5*faktor*(real( ur2(N,m)*rho0(N+1) * conjg(Vpres(m))*redge(N)**2. ))
            
            DO k1=2,N                
                ugr(k1,1)= -g0(k1)*j/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-5,m)
                ugr(k1,1)= ugr(k1,1) + g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-4,m)
                ugr(k1,2)= g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-5,m)
                ugr(k1,2)= ugr(k1,2) - g0(k1)*(j+1)/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-4,m)                
            ENDDO          

            IF(model==5.and.layered) THEN 
                ! Iterating over density jumps in the interior only
                DO k1=1,Nm5l-2     
                    ip=sjmp(k1)
                    BougrPo(nintp) = BougrPo(nintp) - faktor*real(conjg(ugr(ip,1))*v(ip,m,1) + conjg(ugr(ip,2))*v(ip,m,2))*dr*rcent(ip)**2
                ENDDO
            ELSE
              velicinas = -(real( conjg(ugr(:,1)) * v(2:N,m,1) + conjg(ugr(:,2)) * v(2:N,m,2) ))
              ymin=0.; ymax=0.;     
              ! grad(rho) is considered 0 in the first and last half-layers
              ! Note that integrate does not work well for layered models with large discontinuities
              BougrPo(nintp) = BougrPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,2)    
            ENDIF
            
            DO k1=2,N                
                PhiExt(k1)= 4*pi*Gconst/(2*j+1)*rho0(N+1)*(redge(N)**(j+2))*ur2(N,m) /rcent(k1)**(j+1)
                PhiInt(k1)= 4*pi*Gconst/(2*j+1)*rho0(1)*(redge(1)**(1-j))*(ur2(1,m)+extload(m)) *rcent(k1)**j
                velicinas(k1)= 0.5*real ( conjg(ur(k1,m)*(rho0r(k1)-rho0r(k1-1))/dr) * (PhiInt(k1)+PhiExt(k1)+PhiIn(k1,m)+PhiEx(k1,m)) )                
                ! Representing internal density gradient as a series of surface mass densities. Esgrav is just the self-gravity part of EnG
                Esgrav = Esgrav - 0.5*faktor*(real( ur(k1,m)*(rho0r(k1)-rho0r(k1-1)) * conjg(PhiInt(k1)+PhiExt(k1)+PhiIn(k1,m)+PhiEx(k1,m))*rcent(k1)**2. )) 
            ENDDO
            ymin=0.; ymax=0.;
            ! the minus sign is because Phi has opposite sign than in usual notations (here force = grad Phi)
            EdrhoPhi = EdrhoPhi - faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,2)

        ENDDO         ! end of m=0,1,2 cycle
        Io = matmul(Inertia,tidax) * sclJ*eas
        rot_src(nintp) = dot_product(Io, cross_product(omega,tidax)) * eas*sclOmega
        tid_src(nintp) = dot_product(Io, (tidax-taxprev)/dt) * eas
        
        Edisip = advance(Edisip,DD,dt,t,rule)
        CMBTopoEn = advance(CMBTopoEn,CMBTopoPo,dt,t,rule)
        CMBRoPrEn = advance(CMBRoPrEn,CMBRoPrPo,dt,t,rule)
        CMBTiPrEn = advance(CMBTiPrEn,CMBTiPrPo,dt,t,rule)
        CMBSgPrEn = advance(CMBSgPrEn,CMBSgPrPo,dt,t,rule)*selfg
        SuTrEn = advance(SuTrEn,SuTrPo,dt,t,rule)
        BoFoEn = advance(BoFoEn,BoFoPo,dt,t,rule)*selfg
        BougrEn = advance(BougrEn,BougrPo,dt,t,rule)
        BoRoEn = advance(BoRoEn,BoRoPo,dt,t,rule)
        BoTiEn = advance(BoTiEn,BoTiPo,dt,t,rule)
        Eellost = advance(Eellost,DDdot,dt,t,rule)
        Esrc(1) = advance(Esrc(1),rot_src,dt,t,rule)
        Esrc(2) = advance(Esrc(2),tid_src,dt,t,rule)
        
        emplacement = .false.
        IF (pr==1.and.t>=tgrowth.and.virgin) THEN
            emplacement = .true.
            virgin = .false.
        ENDIF
        IF ((t==0.).or.emplacement) THEN
            print *,'Setting reference energies and integrated powers (dissipation, t= ', t/yr,')'
            IF(write_totalE) THEN
               ! E0 values are subtracted from the total values on output
               Esgrav0=0.; Eel0=0.; EnG0=0.; EdrhoPhi0=0.; EdrhoV0=0.;
               SEng0=0.;  SEng_CMB0=0.;  Ssg0=0.;  Ssg_CMB0=0.; Erot0=0.; Etid0=0.;
            ELSE
               Esgrav0=Esgrav; Eel0=Eel; EnG0=EnG(1); Erot0=Erot; EdrhoV0=EdrhoV; Etid0=Etid;
               SEng0=SEng; SEng_CMB0=SEng_CMB; Ssg0=Ssg; Ssg_CMB0=Ssg_CMB; EdrhoPhi0=EdrhoPhi;
               IF (emplacement) THEN
                   ! In order to get a matching term-by-term analysis after smooth load emplacement
                   CMBTopoEn=0.; CMBRoPrEn=0.; CMBTiPrEn=0.; CMBSgPrEn=0.; SuTrEn=0.; BoFoEn=0.; 
                   BougrEn=0.; BoRoEn=0.; BoTiEn=0.; Eellost=0.; Edisip=0.; Esrc=0.;
               ENDIF
            ENDIF
        ENDIF
        DO q=1,nintp-1
            DD(q) = DD(q+1); 
            DDdot(q) = DDdot(q+1);
            SuTrPo(q) = SuTrPo(q+1); 
            CMBTopoPo(q) = CMBTopoPo(q+1);
            CMBSgPrPo(q) = CMBSgPrPo(q+1);
            BoFoPo(q) = BoFoPo(q+1);
            BougrPo(q) = BougrPo(q+1);
            BoRoPo(q) = BoRoPo(q+1);
            BoTiPo(q) = BoTiPo(q+1);
            CMBRoPrPo(q) = CMBRoPrPo(q+1);
            CMBTiPrPo(q) = CMBTiPrPo(q+1);
            rot_src(q) = rot_src(q+1);
            tid_src(q) = tid_src(q+1);
        ENDDO      
        ylastprev(1,:,:) = ylastprev(2,:,:)
        ylastprev(2,:,:) = ylast
        taxprev = tidax
      END SUBROUTINE dissipation
      
      FUNCTION advance(velicina,vykon,dt,time,rule)
      REAL, INTENT(IN) :: velicina,vykon(nintp),dt,time
      REAL :: advance
      REAL, SAVE :: t(nintp)=0.
      INTEGER, SAVE :: step=1       
      INTEGER :: k1
      CHARACTER(*), INTENT(IN) :: rule
        IF(time==0.) THEN
            step=1;
            t(nintp)=0.
        ENDIF        
        IF(time/=t(nintp)) THEN            
            step=step+1                        
            DO k1=1,nintp-1
                t(k1)=t(k1+1)
            ENDDO            
            t(nintp) = time
        ENDIF            
        IF(rule(1:7)=='trapezo') THEN
            IF(step==1) THEN                           
                advance = velicina
            ELSEIF(step==2) THEN
                advance = velicina + vykon(nintp)*dt
            ELSE
                advance = velicina + (vykon(nintp-1)+vykon(nintp))/2.*dt
            ENDIF
        ELSEIF(rule(1:7)=='simpson') THEN
            IF(step==2) THEN
                advance = velicina + vykon(nintp)*dt            
            ELSEIF(mod(step,2)==0) THEN                    
                advance = velicina + Simpson(vykon(nintp-2),vykon(nintp-1),vykon(nintp),t(nintp-2),t(nintp-1),t(nintp))                
            ELSE
                advance = velicina
            ENDIF   
        ELSEIF(rule(1:7)=='simpold') THEN
            IF(step==2) THEN
                advance = velicina + vykon(nintp)*dt            
            ELSEIF(mod(step,2)==0) THEN                    
                advance = velicina + modSim(vykon(nintp-2),vykon(nintp-1),vykon(nintp),t(nintp-2),t(nintp-1),t(nintp))                
            ELSE
                advance = velicina
            ENDIF   
        ELSE
            stop 'unrecognized rule for time integration of quantities'
        ENDIF        
        IF(step/=1 .and. dt<=0.) stop 'dt should be zero only for t==0.'
      END FUNCTION advance
      
      FUNCTION slope(vel,velprev,dt,time)
      COMPLEX, INTENT(IN) :: vel,velprev(2)
      REAL, INTENT(IN) :: dt,time
      COMPLEX :: slope, f0,f1,f2
      REAL :: dx1,dx2,dx
      REAL, SAVE :: t(nintp)=0.
      INTEGER, SAVE :: step=1       
      INTEGER :: k1      
        IF(time==0.) THEN
            step=1;
            t(nintp)=0.
        ENDIF        
        IF(time/=t(nintp)) THEN            
            step=step+1                        
            DO k1=1,nintp-1
                t(k1)=t(k1+1)
            ENDDO            
            t(nintp) = time
        ENDIF      
        
        f0=velprev(1);   f1=velprev(2);   f2=vel;
        dx=t(nintp)-t(nintp-2);
        dx1=t(nintp-1)-t(nintp-2);
        dx2=t(nintp)-t(nintp-1);                
        IF(step==1) THEN                           
            slope = 1.
        ELSEIF(step==2) THEN
            slope = (f2 - f1) / dt
        ELSE            
            slope = f0*dx2/(dx1*dx) - f1*dx/(dx1*dx2) + f2*(dx2+dx)/(dx*dx2)
        ENDIF
        IF(step/=1 .and. dt<=0.) stop 'dt should be zero only for t==0.'        
      END FUNCTION slope
      
      FUNCTION dJdt(J,Jprev,dt,time)      
      REAL, INTENT(IN) :: J(3,3),Jprev(3,3,2),dt,time
      REAL, DIMENSION(3,3) :: dJdt, f0,f1,f2
      REAL :: dx1,dx2,dx
      REAL, SAVE :: t(nintp)=0.
      INTEGER, SAVE :: step=1       
      INTEGER :: k1
        IF(time==0.) THEN
            step=1;
            t(nintp)=0.
        ENDIF        
        IF(time/=t(nintp).and.imprnt==1) THEN            
            step=step+1                        
            DO k1=1,nintp-1
                t(k1)=t(k1+1)
            ENDDO            
            t(nintp) = time
        ENDIF            
        f0=Jprev(:,:,1);   f1=Jprev(:,:,2);   f2=J;
        IF(imprnt==1) THEN
            dx=t(nintp)-t(nintp-2);
            dx1=t(nintp-1)-t(nintp-2);
            dx2=t(nintp)-t(nintp-1);
        ELSE
            dx=time-t(nintp-1);
            dx1=t(nintp)-t(nintp-1);
            dx2=time-t(nintp);
        ENDIF
        
        IF(step<=2) THEN
            IF (dt>0.) THEN
                dJdt = (f2 - f1) / dt
            ELSE
                dJdt = 0.0
            ENDIF
        ELSE            
            dJdt = f0*dx2/(dx1*dx) - f1*dx/(dx1*dx2) + f2*(dx2+dx)/(dx*dx2)
        ENDIF
      END FUNCTION dJdt
      
      SUBROUTINE ab5(yn,fn,t,h,yerr)
      REAL,INTENT(INOUT) :: yn(nvar),t
      REAL,INTENT(IN) :: yerr(nvar)
      REAL :: fn(nvar),h
      REAL, DIMENSION(nvar), SAVE :: ynp1,ynp2,ynp3,ynp4,ynp5,fn0,fn1,fn2,fn3,fn4,fnini(nvar,4)
      INTEGER,SAVE :: k1,ini=0,rozjezd=2
      IF (ini==0) THEN
        h=h/5.;
        SELECT CASE(rozjezd)
        CASE(0)
        ! rozjezd pomoci AB formuli nizsich radu, opraven 21.6.2016 (fn=fn0; added)        
            imprnt=1
            DO k1=1,5
                ynp1=yn+h*fn;  t=t+h;
                call derivs(t,ynp1,fn1)
                ynp2=ynp1+h*(3./2.*fn1 - fn/2.);  t=t+h;
                call derivs(t,ynp2,fn2)
                ynp3=ynp2+h*(23./12.*fn2 - 4./3.*fn1 + 5./12.*fn);  t=t+h;
                call derivs(t,ynp3,fn3)
                ynp4=ynp3+h*(55./24.*fn3 - 59./24.*fn2 + 37./24.*fn1 - 3./8.*fn);  t=t+h;
                call derivs(t,ynp4,fn4)
                ynp5=ynp4+h*(1901./720.*fn4 - 1387./360.*fn3 + 109./30.*fn2 - 637./360.*fn1 + 251./720.*fn)
                yn=ynp5;  t=t+h;
                IF (k1/=5) THEN
                    call derivs(t,yn,fn0); 
                    fnini(:,k1)=fn0; fn=fn0;
                ENDIF;
            ENDDO
            imprnt=0
        CASE(1)
        ! rozjezd, kdy AB formule nizsich radu pouziji jen na prvni ctyri kratsi kroky
        ! problem je, ze ten jeden krok ynp1=yn+h*fn jiz zpusobi onen zkoumany skok v Erot
            imprnt=1
            ynp1=yn+h*fn;  t=t+h;
            call derivs(t,ynp1,fn1)
            ynp2=ynp1+h*(3./2.*fn1 - fn/2.);  t=t+h;
            call derivs(t,ynp2,fn2)
            ynp3=ynp2+h*(23./12.*fn2 - 4./3.*fn1 + 5./12.*fn);  t=t+h;
            call derivs(t,ynp3,fn3)
            ynp4=ynp3+h*(55./24.*fn3 - 59./24.*fn2 + 37./24.*fn1 - 3./8.*fn);  t=t+h;
            call derivs(t,ynp4,fn4)
            ynp5=ynp4+h*(1901./720.*fn4 - 1387./360.*fn3 + 109./30.*fn2 - 637./360.*fn1 + 251./720.*fn)
            yn=ynp5;  t=t+h;
            call derivs(t,yn,fn)
            fnini(:,1)=fn 
            DO k1=1,20
                ynp4=yn; fn0=fn1; fn1=fn2; fn2=fn3; fn3=fn4; fn4=fn; 
                ynp5=ynp4+h*(1901./720.*fn4 - 1387./360.*fn3 + 109./30.*fn2 - 637./360.*fn1 + 251./720.*fn0)
                yn=ynp5; t=t+h;
                IF(k1/=20) call derivs(t,yn,fn);
                IF(k1==5) fnini(:,2)=fn
                IF(k1==10) fnini(:,3)=fn
                IF(k1==15) fnini(:,4)=fn
            ENDDO
            imprnt=0
        CASE(2)
        ! rozjezd pomoci Runge-Kutty na kratsich krocich, 5 je jeden standardni (stejne jako vyse)
            DO k1=1,25                
                call rkck(yn,fn,nvar,t,h,yn,yerr,derivs)
                t=t+h;
                IF(k1/=25) THEN; imprnt=1; call derivs(t,yn,fn); imprnt=0; ENDIF;
                IF(k1==5) fnini(:,1)=fn
                IF(k1==10) fnini(:,2)=fn
                IF(k1==15) fnini(:,3)=fn
                IF(k1==20) fnini(:,4)=fn
            ENDDO
        END SELECT        
        fn1=fnini(:,1); fn2=fnini(:,2); fn3=fnini(:,3); fn4=fnini(:,4);
        ini=1; 
        h=5.*h;     !h je pak jeste pouzito pro hdid (note: skutecne hdid inicializace je ovsem 25.*h) 
      ELSE          
          ynp4=yn; fn0=fn1; fn1=fn2; fn2=fn3; fn3=fn4; fn4=fn; 
          ynp5=ynp4+h*(1901./720.*fn4 - 1387./360.*fn3 + 109./30.*fn2 - 637./360.*fn1 + 251./720.*fn0)
          yn=ynp5; t=t+h;
          ini=2
      ENDIF
      IF (t>=tmax.or.crash==1) ini=0;
      END SUBROUTINE ab5
      
     ! For integrals over f(r)*r^exp
      FUNCTION integrate(rpolo,xmin,xmax,ymin,ymax,velicina,po)
      REAL,INTENT(IN) :: rpolo(:),velicina(:),xmin,xmax,ymin,ymax
      REAL :: integrate,f0,f2,f1,x0,x1,x2
      INTEGER :: po,n,i      
        integrate=0.;  n=size(rpolo);
        DO i=0,n,2
        ! NOTE: rpolo goes from xmax to xmin
            IF (i==0) THEN
                x0=rpolo(2);    x1=rpolo(1);    x2=xmax;
                f2=ymax*x2**po;    f1=velicina(1)*x1**po;     f0=velicina(2)*x0**po;                
            ELSEIF (i==n-1) THEN  
                x0=xmin;    x1=rpolo(n);    x2=rpolo(n-1);
                f0=ymin*x0**po;    f1=velicina(n)*x1**po;     f2=velicina(n-1)*x2**po
            ELSEIF (i==n) THEN            
                ! We cannot use three point rule for the last inerval because number of nodes is not odd    
                x0=xmin;    x1=rpolo(n);
                f0=ymin*x0**po;    f1=velicina(n)*x1**po;
                integrate = integrate + (x1-x0)*(f0+f1)/2.
                exit;
            ELSE
                x0=rpolo(i+2);  x1=rpolo(i-1);  x2=rpolo(i);
                f0=velicina(i+2)*x0**po;    f1=velicina(i-1)*x1**po;    f2=velicina(i)*x2**po;
            ENDIF
            integrate = integrate + Simpson(f0,f1,f2,x0,x1,x2)                        
        ENDDO              
      END FUNCTION
      
      FUNCTION Simpson(F0,F1,F2,x0,x1,x2)
      IMPLICIT NONE
      REAL, INTENT(IN) :: F0,F1,F2,x0,x1,x2
      REAL :: Simpson
      REAL :: h1,h2,dx
          h1=x1-x0;
          h2=x2-x1;
          dx=x2-x0;          
          Simpson=( (2.*h1-h2)*dx/h1*F0 + dx**3/(h1*h2)*F1 + dx*(2.*h2-h1)/h2*F2 )/6.          
      END FUNCTION Simpson
      
      FUNCTION modSim(F0,F1,F2,x0,x1,x2)
      IMPLICIT NONE
      REAL, INTENT(IN) :: F0,F1,F2,x0,x1,x2
      REAL :: modSim
      REAL :: dx1,dx2,dx      
      ! Formula I have derived in hand: Simpson's rule modified for uneven spacing
      ! Numerically not good for some reason, replaced by nicer formula in Simpson (viz notes)
        dx1=x1-x0;  dx2=x2-x1;  dx=x2-x0;
        modSim = F0/(dx1*dx)*intgr(x0,x2,x1,x2) - F1/(dx1*dx2)*intgr(x0,x2,x0,x2) + F2/(dx*dx2)*intgr(x0,x2,x0,x1)      
      END FUNCTION modSim
      
      FUNCTION intgr(A,B,x1,x2)
      IMPLICIT NONE
      REAL :: intgr
      REAL, INTENT(IN) :: A,B,x1,x2      
        intgr = (x1*x2*B + B**3./3. - B**2./2.*(x1+x2)) - (x1*x2*A + A**3./3. - A**2./2.*(x1+x2))        
      END FUNCTION intgr 
      
      FUNCTION icecap(j,m,colatitude,longitude,h,rhoi,alpha)
      INTEGER, INTENT(IN) :: j,m      
      REAL, INTENT(IN) :: colatitude,longitude,h,rhoi,alpha
      REAL :: theta, lambda, faktor
      COMPLEX :: ic=(0.,1.), icecap
      INTEGER :: n
        theta=colatitude*pi/180.; lambda=longitude*pi/180.
        faktor=merge(1.0,2.0,m==0);        
        icecap=4.*pi/(2.*j+1.)*legendre(j,h,rhoi,alpha)*conjg(SphHarm(theta,lambda,j,m))/faktor
        ! comparison with eq. 97 from Martinec2014
        ! if(j==2.and.m==1) print *,icecap, -sqrt(3.*pi/10.)*legendre(j,h,rhoi,alpha)*sin(2*theta)*exp(-ic*lambda)
      END FUNCTION icecap
      
      FUNCTION legendre(n, h, rhoi, alpha)
      REAL :: legendre
      REAL, INTENT(IN) :: h, rhoi, alpha
      REAL :: alfa, x
      INTEGER, INTENT(IN) :: n        
        alfa = alpha * pi/180.
        x = cos(alfa)
        SELECT CASE(load_shape)
        CASE(0) !cap
            legendre = -rhoi*h/(4.-4.*x) * ((cos((n+1)*alfa)-cos((n+2)*alfa))/(n+3./2.) - (cos((n-1)*alfa)-cos(n*alfa))/(n-1./2.))
        CASE(1) !disc
            legendre = rhoi*h/2.*(-(5./2.*x**3 - 3./2.*x) + (x))
        CASE(2) !point
            legendre = Mpoint/(4.*pi*rmax**2.)*(2.*n+1)
        END SELECT
      END FUNCTION legendre

      FUNCTION absv(vektor)
      REAL  :: absv
      REAL , INTENT(IN) :: vektor(:)
      absv=sqrt(dot_product(vektor,vektor))
      END FUNCTION
   
      ! Returns a new normalized rotation vector - must be multiplied by sclOmega to get the physical solution
      FUNCTION new_omg(Itp, Itpold, omega, dt, t, iter)
      IMPLICIT NONE
      REAL, INTENT(IN) :: Itp(3,3), omega(3), Itpold(3,3,2), dt, t
      INTEGER, INTENT(IN) :: iter
      REAL, DIMENSION(3,3) :: Jrate, Qmat, Umat, eigenvec
      REAL, DIMENSION(3) :: nomg, ntid, eigen, new_omg, new_tid, YQsys, PIAmax, PIAmin
      REAL :: dotprod, col, lon, absH, absT
      INTEGER :: idImax, idImin, k1, tmpimp
      CHARACTER(100) :: labiter
      LOGICAL :: test_conv=.false., fix_LOD=.false.
        nomg = omega / absv(omega)
        ntid = tidax / absv(tidax)
        call eige_problem(Itp,eigen,eigenvec)
        idImax = maxloc(eigen, 1)
        idImin = minloc(eigen, 1)
        PIAmax = eigenvec(:,idImax)
        PIAmin = eigenvec(:,idImin)
        ! Flipping the (arbitrary) sign of MAX/MIN eigenvectors: I take the one closer to the rot/tid vector
        dotprod = max(dot_product(eigenvec(:,idImax),nomg), -1.)
        IF(acos(dotprod)/pi*180. > 90.) PIAmax = -eigenvec(:,idImax)
        dotprod = max(dot_product(eigenvec(:,idImin),ntid), -1.)
        IF(acos(dotprod)/pi*180. > 90.) PIAmin = -eigenvec(:,idImin)

        SELECT CASE(wMIA)        
        CASE(1)
            ! setting |w|=1 at first, amplitude is later determined via the conservation of angular momentum
            new_omg(1:3) = PIAmax(1:3)
            ! Angular momentums, absH and absH0, are normalized by a time-constant value sclJ*sclOmega
            absH = absv(matmul(Itp,new_omg)) + Mbody/(1.+1./hostmassf)*adistsq / sclJ
            new_omg = new_omg * absH0 / absH
            IF (tides) THEN
                IF (locked) THEN
                    tidax(1:3) = PIAmin(1:3) * absv(new_omg)
                ELSE
                    absT = absv(matmul(Itp,PIAmin(1:3)))
                    tidax = PIAmin(1:3) * absT0 / absT
                    ! Updating the host-body distance
                    ! adistsq = ( Gconst*Mbody*(1.+hostmassf) / ((absv(tidax)*sclOrbit)**2.) )**(2./3.)
                ENDIF
            ENDIF
        CASE(2)     ! Change of amplitude accounted via m3
            new_omg(1:2) = PIAmax(1:2)
            new_omg(3)   = 1. - (Itp(3,3)-J0(3,3))/J0(3,3)
            new_omg = new_omg * omg0
        CASE(3)     ! Linearized Liouville equation (LLE)
            new_omg(1:2) = Itp(1:2,3) / (J0(3,3)-J0(1,1))            
            new_omg(3)   = 1. - (Itp(3,3)-J0(3,3)) / J0(3,3)
            new_omg = new_omg * omg0
        CASE(4)     ! Extended LLE from Hu et al., 2017a,b
            ! LLE is a perturbation from the hydrostatic state. As such, |Omega| and J0 are linked in the LLE eqs (Amp in Hu_axis is time-constant).
            ! Hu et al. mention updating |Omega| in their algorithm. Note, however, that the (extended) LLE give you the total perturbation at time t,
            ! with therein Omega being the angular speed for which the HYDROSTATIC state is obtained - it is not the time varying |omega(t)|. 
            ! Considering instead (m1,m2,m3) to be a change over dt, i.e. using new_omg = Hu_axis * absv(omega) or Amp = sclOmega * absv(omega) is WRONG.
            ! P.S. (m1,m2,m3) is small even for large angle TPW - the body is always close to hydrostatic. P.P.S. |Omega| plays little role in Eqs 11.
            tmpimp = imprnt;  imprnt = 0;
            Jrate = dJdt(Itp,Itpold,dt,t)
            imprnt = tmpimp
            IF (tides) THEN
                Umat = Umatrix(nomg,tidax/absv(tidax))
                ! tidax and omega amplitudes are initially linked (tid0=omg0), their physical scales differ (eas vs sclOmega)
                new_omg = Hu_axis(sclOmega*omg0,Umat,Itp,Jrate,.false.,dt,Hu_mode) * omg0
                new_tid = Hu_axis(eas*omg0,Umat,Itp,Jrate,.true., dt,Hu_mode) * omg0
                ntid = new_tid/absv(new_tid)
                YQsys = cross_product(new_omg, ntid)
                new_omg = cross_product(ntid, YQsys) / merge(absv(YQsys)/omg0, 1.0, fix_LOD)
                tidax = ntid * absv(new_omg)
                IF(.not.locked) stop 'Hu et al. algorithm is designed for tidally LOCKED bodies only'
            ELSE
                Qmat = Qmatrix(nomg)
                new_omg = Hu_axis(sclOmega*omg0,Qmat,Itp,Jrate,.false.,dt,Hu_mode) * omg0
            ENDIF
            IF(test_conv) THEN
                open(91,file='run/omega.dat',Access='append')
                open(92,file='run/tidax.dat',Access='append')
                write(91,'(4e24.14,a,i0)') t/yr, absv(new_omg), col_lon(new_omg),'  iteration  ',iter
                write(92,'(4e24.14,a,i0)') t/yr, absv(new_tid), col_lon(new_tid),'  iteration  ',iter
                close(91); close(92);
                write(labiter,"(a,i0)") "iteration  ",iter
                call print_eige(Itp,trim(labiter),2,t)
            ENDIF
            !print *,'Final new_omg - omega. time, dt: ',t,dt
            !print *, new_omg - omega
            !print *, absv(new_omg) / absv(omega), angledeg(new_omg,omega)
        CASE DEFAULT
            stop 'unrecognized wMIA option'
        END SELECT

        !col = acos(new_omg(3)/absv(new_omg))/pi*180
        !lon = (atan(new_omg(2)/new_omg(1))/pi*180)
        !print *,'Colatitude, Longitude  ', colat, longit                
      END FUNCTION new_omg

      FUNCTION Hu_axis(Amp,Q,Itp,Jrate,tidal_axis,dt,eqmode)
      REAL, INTENT(IN) :: Amp, Itp(3,3), Jrate(3,3), dt, Q(3,3)
      LOGICAL, INTENT(IN) :: tidal_axis
      INTEGER, INTENT(IN) :: eqmode
      REAL :: ItpQsys(3,3), JrateQsys(3,3), JLLEqs(3,3), CmA, CmB, denom, A2B, omgQsys(3), Hu_axis(3)
      REAL :: Atmp(3,3), A2tmp(3,3)
        JrateQsys = matmul(TRANSPOSE(Q), matmul(Jrate,Q))
        ! ItpQsys is the perturbation from what one gets when the hydrostatic state is transformed to the new coordinate system.
        ItpQsys = matmul(TRANSPOSE(Q), matmul(Itp,Q)) - J0
        ! A,B,C are transformed to X',Y,Z' and treated as anti-matter when LLE is applied to the tidal axis (H. Hu, personal comm.)
        IF (tidal_axis) THEN
            Atmp = ItpQsys; A2tmp = JrateQsys;
            JLLEqs = -1.0*matmul(TRANSPOSE(Smat), matmul(J0,Smat))
            ItpQsys = -1.0*matmul(TRANSPOSE(Smat), matmul(Atmp,Smat))
            JrateQsys = -1.0*matmul(TRANSPOSE(Smat), matmul(A2tmp,Smat))
        ELSE
            JLLEqs = J0
        ENDIF
        CmA = JLLEqs(3,3) - JLLEqs(1,1)
        CmB = JLLEqs(3,3) - JLLEqs(2,2)
        denom = Amp * CmA * CmB
        A2B = 2.*JLLEqs(1,1)*JLLEqs(2,2)
        SELECT CASE(eqmode)
        CASE(1)
            ! Equation 23 from Hu et al.,2017b - extension for SLOW-ROTATING bodies
            omgQsys(1)   = (CmB*ItpQsys(1,3)*(Amp*dt)**2 + 2.*JLLEqs(2,2)*Amp*dt*ItpQsys(2,3)) / A2B
            omgQsys(2)   = (CmA*ItpQsys(2,3)*(Amp*dt)**2 - 2.*JLLEqs(1,1)*Amp*dt*ItpQsys(1,3)) / A2B
            omgQsys(3)   = 1. - ItpQsys(3,3) / JLLEqs(3,3)           
        CASE(2)
            ! Equation 11 from Hu et al., 2017a - extension for FAST-ROTATING bodies
            omgQsys(1)   = ItpQsys(1,3) / CmA + JLLEqs(3,3)*JrateQsys(2,3)/denom
            omgQsys(2)   = ItpQsys(2,3) / CmB + JLLEqs(3,3)*JrateQsys(1,3)/denom
            omgQsys(3)   = 1. - ItpQsys(3,3) / JLLEqs(3,3)
        CASE(3)
            ! LLE
            omgQsys(1)   = ItpQsys(1,3) / CmA
            omgQsys(2)   = ItpQsys(2,3) / CmB
            omgQsys(3)   = 1. - ItpQsys(3,3) / JLLEqs(3,3)
        CASE DEFAULT
            ! LLE (just the PM part), used in the script for Triton's TPW (obtained from H. Hu)
            omgQsys(1)   = ItpQsys(1,3) / CmA
            omgQsys(2)   = ItpQsys(2,3) / CmB
            omgQsys(3)   = 1.
        END SELECT
        IF (tidal_axis) THEN
            Hu_axis = matmul(Q, matmul(Smat, omgQsys))            
        ELSE
            Hu_axis = matmul(Q, omgQsys)
        ENDIF
        !print *,'Hu_pert ', omgQsys(1), omgQsys(2), omgQsys(3)-1.0
      END FUNCTION Hu_axis

      FUNCTION Qmatrix(nomg)
      REAL, INTENT(IN) :: nomg(3)
      REAL :: Qmatrix(3,3)
        Qmatrix(1,1) = nomg(3) + nomg(2)**2/(1.+nomg(3))
        Qmatrix(2,2) = 1. - nomg(2)**2/(1.+nomg(3))
        Qmatrix(3,3) = nomg(3)
        Qmatrix(1,2) = - nomg(1)*nomg(2)/(1.+nomg(3))
        Qmatrix(1,3) = nomg(1)
        Qmatrix(2,3) = nomg(2)
        Qmatrix(2,1) = Qmatrix(1,2)  
        Qmatrix(3,1) = -Qmatrix(1,3)
        Qmatrix(3,2) = -Qmatrix(2,3)
      END FUNCTION Qmatrix

      FUNCTION Umatrix(nwr,nwt)
      REAL, INTENT(IN) :: nwr(3),nwt(3)
      REAL :: Umatrix(3,3)
      REAL :: cp(3)
        cp = cross_product(nwr,nwt)
        Umatrix(1,1) = nwt(1)
        Umatrix(2,1) = nwt(2)
        Umatrix(3,1) = nwt(3)
        Umatrix(1,2) = cp(1)
        Umatrix(2,2) = cp(2)
        Umatrix(3,2) = cp(3)
        Umatrix(1,3) = nwr(1)
        Umatrix(2,3) = nwr(2)
        Umatrix(3,3) = nwr(3)
      END FUNCTION Umatrix

      FUNCTION angledeg(vecA,vecB)
      REAL, INTENT(IN) :: vecA(3), vecB(3)
      REAL :: angledeg
        angledeg = acos(dot_product(vecA,vecB) / (absv(vecA)*absv(vecB))) / pi*180.
      END FUNCTION angledeg

      FUNCTION col_lon(vector)
      REAL, INTENT(IN) :: vector(3)
      REAL :: col_lon(2), colatitude, longitude
        colatitude = acos(vector(3)/absv(vector))/pi*180
        longitude = atan2(vector(2),vector(1))/pi*180
        col_lon(1) = colatitude
        col_lon(2) = longitude
      END FUNCTION col_lon

      SUBROUTINE print_eige(A,label,wo,time)
      REAL, INTENT(IN) :: A(3,3),time
      CHARACTER(*), INTENT(IN) :: label
      INTEGER, INTENT(IN) :: wo
      REAL :: eigen(3),eigenvec(3,3),zunit(3)=(/0.,0.,1./),valsort(3),vecsort(3,3),dotprod,Qbulge(3,3),loadbg(3)
      INTEGER :: idImax, idImin, idImed, k1
      LOGICAL :: verbal=.true.
        call eige_problem(A,eigen,eigenvec)
        idImax = maxloc(eigen, 1)
        idImin = minloc(eigen, 1)
        DO k1 = 1,3
            IF(k1 /= idImax .and. k1 /= idImin) idImed = k1
        ENDDO
        valsort(1) = eigen(idImax)
        valsort(2) = eigen(idImed)
        valsort(3) = eigen(idImin)
        vecsort(:,1) = eigenvec(:,idImax)
        vecsort(:,2) = eigenvec(:,idImed)
        vecsort(:,3) = eigenvec(:,idImin)
        ! Flipping MAX/MIN eigenvectors (sign of eigenvec is arbitrary): I take the ones closer to the rot/tid vector
        dotprod = max(dot_product(vecsort(:,1),yode/absv(yode)), -1.)
        IF(acos(dotprod)/pi*180. > 90.) vecsort(:,1) = -eigenvec(:,idImax)
        dotprod = max(dot_product(vecsort(:,3),tidax/absv(tidax)), -1.)
        IF(acos(dotprod)/pi*180. > 90.) vecsort(:,3) = -eigenvec(:,idImin)
        SELECT CASE(wo)
        CASE(0)
            IF (verbal) THEN
                print *,label
                print *,' MAX: colat, longit, eigenvalue ', col_lon(vecsort(:,1)), valsort(1)
                print *,' MIN: colat, longit, eigenvalue ', col_lon(vecsort(:,3)), valsort(3)
            ENDIF
            open(94,file='run/Jtensors.dat',Access='append')
            write(94,*) label
            write(94,'(3e24.14)') transpose(A)
            close(94)
            IF (label=='PURE LOAD') THEN
                lamp = valsort(1)
                loadax = vecsort(:,1)
                IF(abs(valsort(3)) > abs(valsort(1))) THEN
                    lamp = valsort(3)
                    loadax = vecsort(:,3)
                ENDIF
            ENDIF
            IF (label=='FLOAD + F.B.') THEN
                Qbulge = Umatrix(vecsort(:,1), vecsort(:,3))
                loadbg = matmul(TRANSPOSE(Qbulge), loadax)
                open(91,file='run/MIA.dat',Access='append')
                write(91,'(8e24.14,a)') time/yr,1.0,1.0,1.0,col_lon(loadbg),0.,0.,'  LOAD_in_B.F.  0'
                close(91)
            ENDIF
        CASE(1)
            open(94,file='run/principal.dat',Access='append')
            write(94,'(7e24.14,a,a)') 0.,0.,0.,eigenvec(:,1),eigen(1), '  ', label
            write(94,'(7e24.14,a,a)') 0.,0.,0.,eigenvec(:,2),eigen(2), '  ', label
            write(94,'(7e24.14,a,a)') 0.,0.,0.,eigenvec(:,3),eigen(3), '  ', label
            close(94)
        CASE(2)
            open(91,file='run/MIA.dat',Access='append')
            write(91,'(8e24.14E3,a,a)') time/yr,valsort,col_lon(vecsort(:,1)),col_lon(vecsort(:,3)),'  ',label
            close(91)
        END SELECT
      END SUBROUTINE print_eige
      
      INCLUDE 'check_solver.f90'
      INCLUDE '3Dvisualization.f90'
      INCLUDE 'equations.f90'
      
      END MODULE
