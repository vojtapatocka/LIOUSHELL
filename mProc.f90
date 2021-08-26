      MODULE mProc
      USE mShared
      USE mConst
      USE nr
      IMPLICIT NONE
      
    CONTAINS
    
      SUBROUTINE derivs(t,omegain,dydx)
      REAL :: t, omegain(3)
      REAL, INTENT(OUT) :: dydx(3)
      
        call fcn(nvar,t,omegain,dydx) 
      
      END SUBROUTINE
      
      SUBROUTINE fcn(neq,t,omegain,dydx)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: neq
      REAL :: t, omegain(3)       !INTENT not stated because of quasi_fluid approximation possibility
      REAL, INTENT(OUT) :: dydx(3)
      REAL :: Itp(3,3),rhs(3),ewJw(3),Jw(3),dd,dt,dJdtw(3),omega(3),VtoJcst
      REAL :: dmkl(3),e(2),tau(3),th3,th4,th5,th6,th7,rat(0:jmax),r(0:jmax),tjpo
      REAL, SAVE :: Inertiaprev(3,3,2), dtfakt=-1.
      COMPLEX, DIMENSION(0:jmax) :: dVget,dVgit,Vextload,ko1,ko2,ko3
      COMPLEX ::y(NN,0:jmax), Mem(3*N,0:jmax), dVge(2:N,0:jmax), dVgi(2:N,0:jmax)
      INTEGER :: indxx(3),k1,k2,k3,iterace,m,qfi,qfitot
      LOGICAL :: wesave,timetodraw,runquasifluid      
      
      dt=t-tlast;         
      IF(t>=rozjezd) THEN; nrhse=nrhse+1; IF(nrhse>1) trialdtmax=max(dt,trialdtmax); ENDIF;
      IF ( (memschema==1 .or. memschema==3).and.(dt/=dtfakt) ) THEN
          call faktorizace(dt,0,jmax);
          dtfakt=dt
      ENDIF    
      call cpu_time(th3)
      wesave = (t==tstart).and.(savestate==1).and.(imprnt==1)       
      
      IF(selfg==1 .or. memschema==2 .or. memschema==4) THEN
        iterace = sgit + merge(30,0,dt==0.)
        IF(sgit<1) stop 'sgit must be non-zero: selfgravity term requires iterations'
        IF(memschema==4 .and. sgit<2) stop 'sgit must be > 2: Simpsons rule requires iterations'
      ELSE
        iterace = 0
      ENDIF
      IF(sgres) THEN 
        selfg1=0.; selfg2=0.; Vpres=0.; 
      ENDIF         
      runquasifluid = ((pr==1).or.(couple==1)).and.(quasi_fluid>0)
      qfitot = 1
      IF(runquasifluid) qfitot = qfit
      IF(runquasifluid .AND. dt==0.0) qfitot = 5*qfit
      IF(quasi_fluid==4 .AND. Hu_Venus) qfitot = 1
      ! omegain is normalized, but the below cycle needs omega to be in the physical units
      omega = omegain*sclOmega  

      DO qfi=1,qfitot
          j = jmax    
          tjpo = 2.*j + 1.

          ! We iterate for each m=0,1,2 separately
!$OMP PARALLEL DO
          DO m=0,j        
            DO k3=1,1+iterace                
                IF(selfg==0) THEN; selfg1(:,m)=0.; selfg2(:,m)=0.; Vpres(m)=0.; ENDIF;    
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
                    IF(j==2) y(6*k1-5,m) = -5*sqrt(2./5)*rcent(k1)*phi(j,m,omega)*0.5*(rho0r(k1)+rho0r(k1-1))                                
                    ! selfgravitation
                    y(6*k1-5,m) = y(6*k1-5,m) + selfg1(k1,m); 
                    y(6*k1-4,m) = y(6*k1-4,m) + selfg2(k1,m);
                    ! the zero degree component of centrifugal force
                    IF(j==0) y(6*k1-4,m) = 2*rcent(k1)*phi(j,m,omega)*rho0(k1)
                    ! memory term in the rheological equation
                    y(6*k1-2,m) = Mem(3*k1-2,m)*sclREO; 
                    y(6*k1-1,m) = Mem(3*k1-1,m)*sclREO; 
                    y(6*k1,m) = Mem(3*k1,m)*sclREO;
                ENDDO             
                ! The term drho*varphi (neglected in standard formulation)
                IF (model==5.and.j==2.and.k3>1.and.extended) THEN                    
                    y(6*sjmp(3)-5,m) = y(6*sjmp(3)-5,m) - 5.*sqrt(2./5.)*rcent(sjmp(3))*phi(j,m,omega)*(rho_inp(4)-rho_inp(3))*ur(sjmp(3),m)/dr
                    y(6*sjmp(2)-5,m) = y(6*sjmp(2)-5,m) - 5.*sqrt(2./5.)*rcent(sjmp(2))*phi(j,m,omega)*(rho_inp(3)-rho_inp(2))*ur(sjmp(2),m)/dr
                    y(6*sjmp(1)-5,m) = y(6*sjmp(1)-5,m) - 5.*sqrt(2./5.)*rcent(sjmp(1))*phi(j,m,omega)*(rho_inp(2)-rho_inp(1))*ur(sjmp(1),m)/dr
                ENDIF
                ! Core pressure
                y(6*N+1,m) = -sqrt(j/tjpo)*(rho0(N+1)+rho0(N))*((redge(N)**2)*phi(j,m,omega)+Vpres(m))*sclBC
                y(6*N+2,m) = sqrt((j+1)/tjpo)*(rho0(N+1)+rho0(N))*((redge(N)**2)*phi(j,m,omega)+Vpres(m))*sclBC  
            
                IF(test_reseni==1) ycheck = y(:,m)      
                yhelp(:,1,m) = real(y(:,m)) 
                yhelp(:,2,m) = imag(y(:,m)) 
                
                call cpu_time(th6);
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
                y(:,m)=cmplx(yhelp(:,1,m),yhelp(:,2,m)); call cpu_time(th7);
                
                IF (test_reseni==1 .and. m==0 .and. t==tstart .and. k3==1) THEN;
                    call test_presnosti(y(:,m))   
                    call zpetna_kontrola(y(:,m))      
                    call printAA(y(:,m))
                ENDIF;              
                ur(:,m) = sqrt(j/tjpo)*y(1::6,m) - sqrt((j+1)/tjpo)*y(2::6,m)           ! centres of layers
                ur2(:,m) = 0.5*( ur(1:N,m) + ur(2:N+1,m) )                          ! edges of layers
                ut2(:,m) = sqrt((j+1)/tjpo)*(y(1:6*N-5:6,m)+y(7::6,m))/2. + sqrt(j/tjpo)*(y(2:6*N-4:6,m)+y(8::6,m))/2. 
                IF(memschema==4.and.k3==iterace/2) ymidp(:,m)=y(:,m);
                
                ! SELFGRAVITY
                dVgi(:,m) = 0.; dVge(:,m) = 0.;    
                dVget(m) = 4.*pi*G/tjpo*rho0(N+1)*(redge(N)**(j+2))*ur2(N,m)
                dVgit(m) = 4.*pi*G/tjpo*rho0(1)*(redge(1)**(1-j))*(ur2(1,m)+extload(m))
            
                DO k1=2,N                
                    IF(model==5.and.layered) THEN
                        r(m) = rcent(k1)
                        ko3(m) = 4.*pi*G/tjpo*(rho_inp(4)-rho_inp(3))*ur(sjmp(3),m)
                        ko2(m) = 4.*pi*G/tjpo*(rho_inp(3)-rho_inp(2))*ur(sjmp(2),m)
                        ko1(m) = 4.*pi*G/tjpo*(rho_inp(2)-rho_inp(1))*ur(sjmp(1),m)
                        r_vinp(2) = rcent(sjmp(1)); r_vinp(3) = rcent(sjmp(2)); r_vinp(4) = rcent(sjmp(3));
                        IF(r(m)>r_vinp(2))     THEN; dVge(k1,m)= ko3(m)*r_vinp(4)**(j+2) + ko2(m)*r_vinp(3)**(j+2) + ko1(m)*r_vinp(2)**(j+2)
                        ELSEIF(r(m)>r_vinp(3)) THEN; dVgi(k1,m)= ko1(m)*r_vinp(2)**(1-j)
                                                     dVge(k1,m)= ko3(m)*r_vinp(4)**(j+2) + ko2(m)*r_vinp(3)**(j+2) 
                        ELSEIF(r(m)>r_vinp(4)) THEN; dVgi(k1,m)= ko1(m)*r_vinp(2)**(1-j) + ko2(m)*r_vinp(3)**(1-j) 
                                                     dVge(k1,m)= ko3(m)*r_vinp(4)**(j+2)
                        ELSEIF(r(m)>rmin) THEN;      dVgi(k1,m)= ko1(m)*r_vinp(2)**(1-j) + ko2(m)*r_vinp(3)**(1-j) + ko3(m)*r_vinp(4)**(1-j);
                        ENDIF

                    ! Model 1 is homogeneous -> density jumps only at the surface and CMB
                    ELSEIF (model/=1) THEN 
                      DO k2 = k1,N-1
                        rat(m) = rcent(k2+1)/rcent(k2)                    
                        dVge(k1,m)=dVge(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/(dr*(j+3.))*(rcent(k2)**(j+3.)*(1.-rat(m)**(j+3.)))&
                        &-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr*(j+4.))*(rcent(k2)**(j+4.)*(1.-rat(m)**(j+4.))))*4.*pi*G/tjpo                    
                        !computationally faster option using Taylor expansion (could be very accurate if rat**(j+3) differs little from 1)
                        !dVge(k1,m)=dVge(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/(dr*(j+3))*(rcent(k2)**(j+3)*(j+3)*(1.-rat))&
                        !&-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr*(j+4))*(rcent(k2)**(j+4)*(j+4)*(1.-rat)))*4*pi*G/tjpo
                        !print *,rcent(k2)**(j+4)*(1.-rat**(j+4)),rcent(k2)**(j+4)*(j+4)*(1.-rat),(rcent(k2)**(j+4)-rcent(k2+1)**(j+4))
                      ENDDO
                      DO k2 = 2,k1-1            ! Shape for j=2
                        dVgi(k1,m) = dVgi(k1,m)+((ur(k2,m)+(ur(k2+1,m)-ur(k2,m))*rcent(k2)/dr)*(rho0(k2+1)-rho0(k2))/dr*log(rcent(k2)/rcent(k2+1))&
                        &-(ur(k2+1,m)-ur(k2,m))*(rho0(k2+1)-rho0(k2))/(dr*dr)*(rcent(k2)-rcent(k2+1)))*4.*pi*G/tjpo
                      ENDDO
                        dVgi(k1,m) = dVgi(k1,m)+((ur(1,m)+(ur(2,m)-ur(1,m))*rcent(1)/dr)*2*(rho0(2)-rho0(1))/dr*log(redge(1)/rcent(2))&
                        &-(ur(2,m)-ur(1,m))*2*(rho0(2)-rho0(1))/(dr*dr)*(redge(1)-rcent(2)))*4.*pi*G/tjpo
                    ENDIF
                    
                    PhiIn(:,m)=dVgi(:,m)*rcent(2:N)**j; PhiEx(:,m)=dVge(:,m)/rcent(2:N)**(j+1);                
                    selfg1(k1,m)=-0.5*(rho0r(k1)+rho0r(k1-1))*sqrt(tjpo*j)*rcent(k1)**(j-1)*(dVgi(k1,m)+Dvgit(m))
                    selfg2(k1,m)=-0.5*(rho0r(k1)+rho0r(k1-1))*sqrt(tjpo*(j+1))/(rcent(k1)**(2+j))*(dVge(k1,m)+Dvget(m))
                ENDDO
                ! Term drho*Phi (neglected in the standard formulation)
                IF (model==5.and.j==2.and.k3>1.and.extended) THEN                               
                    selfg1(sjmp(3),m) = selfg1(sjmp(3),m) - &
                        (rho_inp(4)-rho_inp(3))/dr*ur(sjmp(3),m)*sqrt(tjpo*j)*rcent(sjmp(3))**(j-1)*(dVgi(sjmp(3),m)+Dvgit(m))
                    selfg2(sjmp(3),m) = selfg2(sjmp(3),m) - &
                        (rho_inp(4)-rho_inp(3))/dr*ur(sjmp(3),m)*sqrt(tjpo*(j+1))/(rcent(sjmp(3))**(2+j))*(dVge(sjmp(3),m)+Dvget(m))    
                    selfg1(sjmp(2),m) = selfg1(sjmp(2),m) - &
                        (rho_inp(3)-rho_inp(2))/dr*ur(sjmp(2),m)*sqrt(tjpo*j)*rcent(sjmp(2))**(j-1)*(dVgi(sjmp(2),m)+Dvgit(m))
                    selfg2(sjmp(2),m) = selfg2(sjmp(2),m) - &
                        (rho_inp(3)-rho_inp(2))/dr*ur(sjmp(2),m)*sqrt(tjpo*(j+1))/(rcent(sjmp(2))**(2+j))*(dVge(sjmp(2),m)+Dvget(m))    
                    selfg1(sjmp(1),m) = selfg1(sjmp(1),m) - &
                        (rho_inp(2)-rho_inp(1))/dr*ur(sjmp(1),m)*sqrt(tjpo*j)*rcent(sjmp(1))**(j-1)*(dVgi(sjmp(1),m)+Dvgit(m))
                    selfg2(sjmp(1),m) = selfg2(sjmp(1),m) - &
                        (rho_inp(2)-rho_inp(1))/dr*ur(sjmp(1),m)*sqrt(tjpo*(j+1))/(rcent(sjmp(1))**(2+j))*(dVge(sjmp(1),m)+Dvget(m))    
                ENDIF
                
                Vpres(m)=dVget(m)/redge(N)**(j+1) + (dVgi(N,m)+dVgit(m))*redge(N)**j
                Vsurf(m)=(dVget(m)+dVge(2,m))/redge(1)**(j+1) + dVgit(m)*redge(1)**j
                ! Selfgravity potential
                Ek(m)=(dVget(m)+dVge(2,m))/(4*pi*G/tjpo) + rho0(1)*(redge(1)**(j+2))*ur2(1,m)
                
                selfg1Tr(m,2)=rho0(N)*sqrt(tjpo*j) * (redge(N)**(j-1.)) * (dVgi(N,m)+Dvgit(m));
                selfg2Tr(m,2)=rho0(N)*sqrt(tjpo*(j+1)) / (redge(N)**(2.+j)) * Dvget(m);
                selfg1Tr(m,1)=rho0(1)*sqrt(tjpo*j) * (redge(1)**(j-1.)) * Dvgit(m); 
                selfg2Tr(m,1)=rho0(1)*sqrt(tjpo*(j+1)) / (redge(1)**(2.+j)) * (dVge(2,m)+Dvget(m));
            
                taasolve=taasolve+th7-th6
            ENDDO   ! end of selfg iterations
                
                IF (imprnt==1 .and. qfi==qfitot) THEN
                    DO k1=1,3;
                        SELECT CASE(memschema)
                        CASE(0)
                            Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*ylast((3+k1):NN:6,m)/visko
                        CASE(1,2)
                            Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*(ylast((3+k1):NN:6,m)+y((3+k1):NN:6,m))/visko/2
                        CASE(3)
                            Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + dt*y((3+k1):NN:6,m)/visko
                        CASE(4)                        
                            Memlast(k1:3*N:3,m) = Memlast(k1:3*N:3,m) + &
                            dt/6.*(ylast((3+k1):NN:6,m) + 4.*ymidp((3+k1):NN:6,m) + y((3+k1):NN:6,m))/visko 
                        END SELECT                   
                    ENDDO;
                    ylast(:,m)=y(:,m); tlast=t; 
                    timetodraw = ((t>=tstart*2.**real(npicsav-NpicSavMax+1).or.t==0.).and.pr==0)
                    if(loadstate==1) timetodraw = (t>=npicsav*tmax/NpicSavMax)
                    IF (modd/=0.and.timetodraw) call rez_radialni_slozky_v(y(:,m),Ntheta,Nphi,N,j,m,modd)    
                    IF (.not.isotest) k2t(m)=real(Vsurf(m)/(phi(j,m,omega)*redge(1)**2.)) 
                    Vextload(m) = 4.*pi*G/tjpo*rho0(1)*rmax*extload(m)
                    IF (isotest) k2t(m)=real( (Vsurf(m)-Vextload(m)) / Vextload(m) )
                ENDIF
          ENDDO     ! end of cycle over m=0,1,2
!$OMP END PARALLEL DO
                 
            VtoJcst=(5.*redge(1)**3.)/(2*pi*G)
            Itp(1,1)=(sqrt(pi/5)*VtoJcst/3*real(Vsurf(0))-VtoJcst*sqrt(2*pi/15)*real(Vsurf(2)))
            Itp(2,2)=(sqrt(pi/5)*VtoJcst/3*real(Vsurf(0))+VtoJcst*sqrt(2*pi/15)*real(Vsurf(2)))
            Itp(3,3)=(-2*sqrt(pi/5)*VtoJcst/3*real(Vsurf(0)))
            Itp(1,3)=VtoJcst*sqrt(2*pi/15)*real(Vsurf(1))
            Itp(2,3)=-VtoJcst*sqrt(2*pi/15)*imag(Vsurf(1))
            Itp(1,2)=VtoJcst*sqrt(2*pi/15)*imag(Vsurf(2))  
            Itp(2,1)=Itp(1,2); Itp(3,1)=Itp(1,3); Itp(3,2)=Itp(2,3);
            !call print_eige(Itp,'ROTATIONAL BULGE',.false.); stop 'Rot Bulge';
            Itp = Itp + Ball
            
            IF(pr==1 .and. fixload>0) THEN
                IF(iceslow) THEN
                    Itp = Itp + Iload * merge(1., sin(t/tgrowth*pi/2.), t>=tgrowth)
                ELSE
                    Itp = Itp + Iload
                ENDIF
            ENDIF

            ! Replacing omegain. The equilibrium position of the rotation vector to which iterations converge
            ! balances inertia contribution of the load with that of the unrelaxed part of the rotational bulge
            IF(runquasifluid) omegain = fluid_axis(Itp,Inertiaprev,omega,dt,t)
            omega = omegain*sclOmega

      ! end of quasi-fluid iterations
      ENDDO       
        
        IF(imprnt==1)  THEN
            Inertia = Itp/sclJ              
            Erot = dot_product(matmul(Itp,omega),omega)/2.
            ! soucet dvou clenu nize nam dela "skok" v Erot pri rozjezdu pomoci AB (oproti Erot ktere dostanu kdyz pouzivam RK4)
            ! write(28,'(20e30.17)') t/yr,omega(1)*dot_product(Inertia(:,1),omega),omega(2)*dot_product(Inertia(:,2),omega),Erot            
            ! solved: zpusobeno explicitnim krokem s presnosti radu 1, ktery je pri rozjezdu vyuzit (Erot je extremne citliva na omegu)
        ENDIF

        IF (imprnt==1) call dissipation(omega,dt,t)
        IF (imprnt==1 .and. modd/=0 .and. timetodraw ) THEN
          call zapis3D(N,Ntheta,Nphi,npicsav); npicsav=npicsav+1;
        ENDIF;
        call cpu_time(th4)

        ! Vsurf and phi in the code both have unconvetional signs, but in the Love numbers (ratios) it cancels out
        IF(pr==0.and.t==tstart.and.(.not.isotest)) THEN; print *,'fluid tidal k2: ',k2t(0); k2Tf=k2t(0); ENDIF;
        IF(pr==0.and.t==0..and.isotest) THEN; print *,'loading k2 at time=0: ',k2t(0); k2Le=k2t(0); ENDIF;
        IF(pr==0.and.t==0..and.(.not.isotest)) THEN; print *,'tidal k2 at time=0: ',k2t(0); k2Te=k2t(0); ENDIF;
        
        omega=omega/sclOmega; 
        Itp=Itp/sclJ;
        
        IF (tcross==1) THEN; call eige_problem(Itp,eigenumber,eigevec,ilibb);  ENDIF;
        Jw=matmul(Itp,omega)
        IF (tcross/=1) ewJw=cross_product(omega,Jw)*sclOmega        
        IF (tcross==1) ewJw=matmul(eigenumber*matmul(omega,eigevec),tcross_product(omega,eigevec))*sclOmega
        
        !dJdtw=matmul((Itp-Inertia)/dt,omega)
        dJdtw=matmul(dJdt(Itp,Inertiaprev,dt,t),omega)

        IF (dt/=0.) THEN                
            rhs=-(dJdtw+ewJw)      !multiplied by -1 in order to use directly Itp in the solver
            IF(imprnt==1) dJdtw0=dJdtw
        ELSE
            rhs=-(dJdtw0+ewJw)     !initial condition dJdtw0 is either zero or loaded from inistate
        ENDIF
        IF (imprnt==1) THEN
            Inertiaprev(:,:,1) = Inertiaprev(:,:,2)
            Inertiaprev(:,:,2) = Inertia
        ENDIF
           
        SELECT CASE(ilibb)
        CASE(1)
            call ludcmp(Itp,3,3,indxx,dd)
            IF (crash==1) return
            call lubksb(Itp,3,3,indxx,rhs)
        CASE(2)
            stop 'disabled: MKL needed'
            !call dgetrf(3,3,Itp,3,indxx,info)
            !call dgetrs('N',3,1,Itp,3,indxx,rhs,3,info)
        CASE(3)
            stop 'disabled: MKL needed'
            !IF (lwork==-1) THEN
            !    IF (.not. allocated(work)) allocate(work(1));
            !    call dsytrf('U',3,Itp,3,indxx,work,lwork,info)
            !    lwork=int(work(1)); deallocate(work); allocate(work(lwork));
            !ENDIF;
            !call dsytrf('U',3,Itp,3,indxx,work,lwork,info)
            !call dsytrs('U', 3, 1, Itp, 3, indxx, rhs, 3, info )
        END SELECT
        dydx=rhs         
        
        call cpu_time(th5)
        tnapln=tnapln+th4-th3;      tliou=tliou+th5-th4;
        IF (t==tstart.and.pr==0) nstp=1;
      
      END SUBROUTINE fcn

      SUBROUTINE faktorizace(dt,inicializace,j)
      REAL :: dt
      INTEGER :: inicializace,j,i1,i2,k1,k2
      call cpu_time(th1)
      IF (inicializace==1) THEN
        call napln_A(j,A)       !diky deleni reovztahu viskozitou mame zavislost na r jen u diagonalnich prvku matice A, ty proto
        call napln_B(j,B)       !doplnime pozdeji, abychom nevolali plneni matice prilis casto. 
        
        IF(ilib/=2.or.opt/=0) THEN
            AA=0.;
            IF(model/=4) THEN; call BC_quasifreesurf(j,'outer'); call BC_quasifreesurf(j,'inner'); ENDIF; 
            IF(model==4) THEN; call BC_quasifreesurf(j,'outer'); call BC_noslip('inner'); ENDIF; 
        ELSE
            abini=0.;    
            IF(model/=4) THEN; call BCab_quasifreesurf(j,'outer'); call BCab_quasifreesurf(j,'inner'); ENDIF; 
            IF(model==4) THEN; call BCab_quasifreesurf(j,'outer'); call BCab_noslip('inner'); ENDIF; 
        ENDIF
        DO k1=1,N       
            DO k2=4,6   ! Completing the diagonal elements of A + inclusion of shear modulus, by which we divide in reo
                IF (memschema==1) THEN      
                    ! Should be -r/shear, but in rheology (discretization in r) we division by r skip
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
           
        IF (test_reseni==1) THEN;   AAcheck=AA;     ENDIF;      
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
        call fcn(nvar,x,y,f1);
        DO i=1,n; 
            call fcn(nvar,x,y+deltay(:,i),f3(:,i)); 
            dfdy(:,i)=(f3(:,i)-f1)/deltay(i,i);   
        ENDDO;      
        imprnt=1;
      END SUBROUTINE      
        
      SUBROUTINE jac(x,y,dfdx,dfdy,n)
      INTEGER :: n,i
      REAL :: x,y(n),dfdx(n),dfdy(n,n),f2(n),f1(n),f3(n,n)
      REAL :: deltax,deltay(n,n)
        imprnt=0;
        deltax=yr/1000; deltay=0.;
        DO i=1,n; deltay(i,i)=maxval(abs(y))/1000; ENDDO;
        call derivs(x,y,f1); call derivs(x+deltax,y,f2); 
        DO i=1,n; call derivs(x,y+deltay(:,i),f3(:,i)); dfdy(:,i)=(f3(:,i)-f1)/deltay(i,i);   ENDDO;
        dfdx=(f2-f1)/deltax;       
      END SUBROUTINE

      FUNCTION phi(j,m,w)
      INTEGER :: j,m
      REAL :: w(3)
      COMPLEX :: phi
      IF (j==2) THEN
        SELECT CASE(m)
        CASE(0) ; phi=sqrt(pi/5)*(w(1)*w(1)+w(2)*w(2)-2*w(3)*w(3))/3;
        CASE(1) ; phi=sqrt(2*pi/15)*(w(1)*w(3)-(0.,1.)*w(2)*w(3));
        CASE(2) ; phi=sqrt(pi/30)*(w(2)*w(2)-w(1)*w(1)+2*w(1)*w(2)*(0.,1.));
        END SELECT  
      ELSEIF (j==0) THEN
        phi=2.*sqrt(pi)*dot_product(w,w)/3.
      ENDIF
      END FUNCTION
      
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

      IF(selfg==0.and.model==1.and..false.) THEN
          IF(extload(1)/=0. .or. extload(2)/=0.) stop 'EnG not analyticaly computed for m/=0'
          ura = (ur2(1,0)+extload(0))*sqrt(5/pi)/4
          SEng = 2.*pi*rho0(1)*g0(1)*8./35.*(3.*ura**4+8.*ura**3*rmax/3.+7.*ura**2*rmax**2/2.)
          ura = ur2(N,0)*sqrt(5/pi)/4
          SEng_CMB = HomoEnG+2*pi*rho0(N+1)*g0(N+1)*8/35.*(3*ura**4+8*ura**3*rmin/3+7*ura**2*rmin**2/2)
          HomoEnG = SEng + SEng_CMB
          EnG = HomoEnG
      ELSEIF (mode==-1.or.mode==0) THEN
          GE=0.; chyba=0.; SEng=0.; Vsg=0.; Ssg=0.; SEng_CMB=0.; Ssg_CMB=0.; EdrhoV=0.;
          Vcontrol=0.; GEcontrol=0.; GEcontrol2=0.; Vcontrol2=0.;
          DO k2=0,Ntheta-1                  !bereme jen realnou cast harmonik
              
            theta=k2*pi/Ntheta            
            DO k3=0,Nphi-1 
                phi=2.*k3*pi/Nphi
                
                !POZOR: vyscitani pres harmoniky m=-1, m=-2 probiha ve funkci SphHarm
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
          EnG = SEng + SEng_CMB + EdrhoV + Egrav*selfg
      ENDIF
      
      END FUNCTION EnG
      
      FUNCTION G_pot(r)
      REAL :: G_pot,r,Mbig,Msmall
      
        MBig=4.*pi*rho0(1)/3.*rmax**3.
        Msmall=4.*pi*(rhocore-rho0(1))/3.*rmin**3
        IF(r>=rmax) THEN
            G_pot=-G*(Mbig+Msmall)/r    
        ELSEIF(r>=rmin) THEN
            G_pot=-G*Mbig/(2.*rmax**3.)*(3.*rmax**2.-r**2.)-G*Msmall/r
        ELSE
            G_pot=-G*Mbig/(2.*rmax**3.)*(3.*rmax**2.-r**2.)-G*Msmall/(2.*rmin**3.)*(3.*rmin**2.-r**2.)
        ENDIF 
        !G_pot=G_pot+G*(Mbig+Msmall)/rmax 
        
      END FUNCTION    
      
      SUBROUTINE eige_problem(J,eigen,eigev,eigelib)
      REAL, INTENT (IN) :: J(3,3)
      REAL :: Itp(3,3),eigen(3),eigev(3,3)
      INTEGER, INTENT(IN) :: eigelib
      REAL, SAVE :: vl,vu
      INTEGER, SAVE :: il,iu,mfound,isuppz(6)
      Itp=J
        SELECT CASE(eigelib)
        CASE(1)
            call jacobi(Itp,3,3,eigen,eigev,nrot)
        CASE(3)
            stop 'disabled: MKL needed'
            !IF (lworkj==-1) THEN
            !    allocate(workj(1)); allocate(iwork(1)); vl=2*dlamch('S');
            !    call dsyevr('V','A','U',3 ,Itp, 3, vl, vu, il, iu, vl, mfound, eigen, eigev, 3, isuppz, workj, lworkj, iwork, liwork, info)
            !    lworkj=int(workj(1)); deallocate(workj); allocate(workj(lworkj));
            !    liwork=iwork(1);  deallocate(iwork); allocate(iwork(liwork));
            !ENDIF;
            !call dsyevr('V','A','U',3 ,Itp, 3, vl, vu, il, iu, Tiny, mfound, eigen, eigev, 3, isuppz, workj, lworkj, iwork, liwork, info)
        END SELECT        
      END SUBROUTINE
      
      SUBROUTINE dissipation(omega,dt,t)
      IMPLICIT NONE
      REAL, INTENT(IN) :: omega(3),dt,t
      REAL :: faktor,ymin,ymax,velicina(2:N-1),velicinas(2:N),rmax2,rmin2
      REAL, SAVE :: BoFoPo(nintp),BoRoPo(nintp),SuTrPo(nintp),CMBTopoPo(nintp),BougrPo(nintp)
      REAL, SAVE :: CMBSgPrPo(nintp),CMBRoPrPo(nintp),DD(nintp),DDdot(nintp)
      COMPLEX :: v(N+1,0:jmax,2),v2(N,0:jmax,2),Ddotvek(1:N),rot_pom,vrad(N,0:jmax)
      COMPLEX :: ugr(2:N,2),PhiExt(2:N),PhiInt(2:N),trakce_rad,pcore
      COMPLEX, SAVE :: ylastprev(2,NN,0:jmax)
      INTEGER :: m,q,ia,i,ip,k1      
      REAL,SAVE :: Ellastfuck
      LOGICAL, SAVE :: nullit=.true.
      
      rmax2=rmax**2.; rmin2=rmin**2.;
      DD(nintp)=0.; DDdot(nintp)=0.; ErotHarm=0.; Egrav=0.; Edef=0.; Eel=0.; EdrhoPhi=0.;
      SuTrPo(nintp)=0.; CMBTopoPo(nintp)=0.; BoFoPo(nintp)=0.; CMBSgPrPo(nintp)=0.;
      BoRoPo(nintp)=0.; CMBRoPrPo(nintp)=0.; BougrPo(nintp)=0.
      ! The integral quantities are never smooth over the pr=0 -> pr=1 transition (everything is restarted)
      IF(t==0.) ylastprev=(0.,0.)
      
      DO m=0,jmax
          IF(order_analyse .and. m/=mp) cycle;
          ip=1;
          DO i=1,NN,6
              ! Not using more sophisticated calculation of derivative leads to wrong integration for wobbling process
              v(ip,m,1) = slope(ylast(i,m), ylastprev(:,i,m), dt, t)
              v(ip,m,2) = slope(ylast(i+1,m), ylastprev(:,i+1,m), dt, t)
              ip=ip+1    
          ENDDO
          
          v2(:,m,1) = (v(1:N,m,1)+v(2:N+1,m,1))/2.
          v2(:,m,2) = (v(1:N,m,2)+v(2:N+1,m,2))/2.
          
          vrad(:,m) = sqrt(j/(2.*j+1))*v2(:,m,1) - sqrt((j+1)/(2.*j+1))*v2(:,m,2)
          faktor=merge(1.,2.,m==0)  
          
          DO q=0,2
          ! Summing over our three components of deviatoric stress
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
          
          rot_pom = -5.*sqrt(2./5)*phi(j,m,omega)   !tvar pro silu ma jeste *rho0(k1)*rcent(k1)
          ymin = -real( rot_pom*rho0(N) * conjg(v2(N,m,1)) )
          ymax = -real( rot_pom*rho0(1) * conjg(v2(1,m,1)) )
          velicinas = -real( rot_pom*rho0(2:N) * conjg(v(2:N,m,1)) )
          BoRoPo(nintp) = BoRoPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,3)
          
          !SuTrPo(nintp) = SuTrPo(nintp) - real( (ur2(1,m)+extload(m))*rho0(1)*g0(1) * conjg(vrad(1,m)) )*rmax2*faktor
          !CMBTopoPo(nintp) = CMBTopoPo(nintp) - real( conjg(ur2(N,m)) * rho0(N+1)*g0(N+1)*vrad(N,m) )*rmin2*faktor
          trakce_rad = a0*(a1*ylast(3,m)+a2*ylast(4,m)+a3*ylast(5,m)) + b0*(b1*ylast(3,m)+b3*ylast(5,m)+b4*ylast(6,m))
          SuTrPo(nintp) = SuTrPo(nintp) + real( trakce_rad * conjg(vrad(1,m)) )*rmax2*faktor
          pcore = rhocore*(phi(j,m,omega)*(rmin2) + Vpres(m))     !scalar field, radial component of core pressure vector        
          trakce_rad = a0*(a1*ylast(6*N-3,m)+a2*ylast(6*N-2,m)+a3*ylast(6*N-1,m)) &
                    + b0*(b1*ylast(6*N-3,m)+b3*ylast(6*N-1,m)+b4*ylast(6*N,m))          
          CMBTopoPo(nintp) = CMBTopoPo(nintp) + real( conjg(-trakce_rad-pcore)*vrad(N,m) )*rmin2*faktor
          
          CMBSgPrPo(nintp) = CMBSgPrPo(nintp) + real( (rho0(N)+rho0(N+1))*conjg(Vpres(m)) * vrad(N,m) )*rmin2*faktor
          CMBRoPrPo(nintp) = CMBRoPrPo(nintp) + real( (rho0(N)+rho0(N+1))*rmin2*conjg(phi(j,m,omega)) * vrad(N,m) )*rmin2*faktor
          
          ErotHarm = ErotHarm + real( phi(j,m,omega) * conjg(Ek(m)) )*faktor
          Egrav = Egrav - 0.5*faktor*(real( (ur2(1,m)+extload(m))*rho0(1) * conjg(Vsurf(m))*redge(1)**2. ))
          Egrav = Egrav - 0.5*faktor*(real( ur2(N,m)*rho0(N+1) * conjg(Vpres(m))*redge(N)**2. ))  
          ! prvotni vypocet v poznamkach povazoval toto za jediny prispevek ke gravitacni energii, je jich ale vice
          
          DO k1=2,N                
                ugr(k1,1)= -g0(k1)*j/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-5,m)
                ugr(k1,1)= ugr(k1,1) + g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-4,m)
                ugr(k1,2)= g0(k1)*sqrt(j*(j+1.))/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-5,m)
                ugr(k1,2)= ugr(k1,2) - g0(k1)*(j+1)/(2.*j+1)*(rho0r(k1)-rho0r(k1-1))/dr * ylast(6*k1-4,m)                
          ENDDO          
          ! integrate does not work well for functions with discontinuities (for layered models)
          IF(model==5.and.layered) THEN 
              DO k1=1,3     ! number of internal density jumps
                ip=sjmp(k1)
                BougrPo(nintp) = BougrPo(nintp) - faktor*real(conjg(ugr(ip,1))*v(ip,m,1) + conjg(ugr(ip,2))*v(ip,m,2))*dr*rcent(ip)**2
              ENDDO
          ELSE
            velicinas = -(real( conjg(ugr(:,1)) * v(2:N,m,1) + conjg(ugr(:,2)) * v(2:N,m,2) ))
            ymin=0.; ymax=0.;     ! rho uvazujeme v prvni a posledni pulvrstve konstantni, tudiz grad(rho) je tam 0          
            BougrPo(nintp) = BougrPo(nintp) + faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,2)    
          ENDIF
          
          DO k1=2,N                
                PhiExt(k1)= 4*pi*G/(2*j+1)*rho0(N+1)*(redge(N)**(j+2))*ur2(N,m) /rcent(k1)**(j+1)
                PhiInt(k1)= 4*pi*G/(2*j+1)*rho0(1)*(redge(1)**(1-j))*(ur2(1,m)+extload(m)) *rcent(k1)**j
                velicinas(k1)= 0.5*real ( conjg(ur(k1,m)*(rho0r(k1)-rho0r(k1-1))/dr) * (PhiInt(k1)+PhiExt(k1)+PhiIn(k1,m)+PhiEx(k1,m)) )                
                Egrav = Egrav - 0.5*faktor*(real( ur(k1,m)*(rho0r(k1)-rho0r(k1-1)) * conjg(PhiInt(k1)+PhiExt(k1)+PhiIn(k1,m)+PhiEx(k1,m))*rcent(k1)**2. )) 
                ! alternativni zpusob. analogie vypoctu obdobneho clenu pro topografie (Egrav), spocivajici v pouziti plosne hustoty na stupni 2
          ENDDO
          ymin=0.; ymax=0.;
          ! the minus sign is because Phi has opposite sign than in usual notations (here force = grad Phi)
          EdrhoPhi = EdrhoPhi - faktor*integrate(rcent(2:N),rmin,rmax,ymin,ymax,velicinas,2)
          
      ENDDO         ! end of m=0,1,2 cycle
      
      ! Rotational energy of the core is included here through Ih(3)      
      ErotHarm = ErotHarm + 3.*sqrt(pi)*Ih(3)*real(phi(0,0,omega))                 
      Edisip = advance(Edisip,DD,dt,t,rule)
      CMBTopoEn = advance(CMBTopoEn,CMBTopoPo,dt,t,rule)
      CMBRoPrEn = advance(CMBRoPrEn,CMBRoPrPo,dt,t,rule)
      CMBSgPrEn = advance(CMBSgPrEn,CMBSgPrPo,dt,t,rule)*selfg
      SuTrEn = advance(SuTrEn,SuTrPo,dt,t,rule)
      BoFoEn = advance(BoFoEn,BoFoPo,dt,t,rule)*selfg
      BougrEn = advance(BougrEn,BougrPo,dt,t,rule)
      BoRoEn = advance(BoRoEn,BoRoPo,dt,t,rule)
      Eellost = advance(Eellost,DDdot,dt,t,rule)         
      
      IF ((t==0.).or.(pr==1.and.icetest.and.iceslow.and.t>tgrowth.and.nullit)) THEN
          IF(write_totalE) THEN
             ! E0 values are subtracted from the total values on output
             ErotHarm0=0.; Egrav0=0.; Eel0=0.; EnG0=0.; EdrhoPhi0=0.; EdrhoV0=0.;
             SEng0=0.;  SEng_CMB0=0.;  Ssg0=0.;  Ssg_CMB0=0.; Erot0=0.;
          ELSE
             ErotHarm0=ErotHarm; Egrav0=Egrav; Eel0=Eel; EnG0=EnG(1); Erot0=Erot; EdrhoV0=EdrhoV;
             SEng0=SEng; SEng_CMB0=SEng_CMB; Ssg0=Ssg; Ssg_CMB0=Ssg_CMB; EdrhoPhi0=EdrhoPhi;
             ! In order to get matching terms
             IF (pr==1.and.icetest.and.iceslow.and.t>tgrowth.and.nullit) THEN
                 CMBTopoEn=0.; CMBRoPrEn=0.; CMBSgPrEn=0.; SuTrEn=0.; BoFoEn=0.; 
                 BougrEn=0.; BoRoEn=0.; Eellost=0.; Edisip=0.;
                 nullit=.false.
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
        CMBRoPrPo(q) = CMBRoPrPo(q+1);
      ENDDO      
      ylastprev(1,:,:)=ylastprev(2,:,:)
      ylastprev(2,:,:)=ylast      
      Ellastfuck=Eel
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
            dJdt = (f2 - f1) / dt
        ELSE            
            dJdt = f0*dx2/(dx1*dx) - f1*dx/(dx1*dx2) + f2*(dx2+dx)/(dx*dx2)
        ENDIF
        IF(step/=1 .and. dt<=0.) stop 'dt should be zero only for t==0.'
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
            DO k1=1,5
              imprnt=1;
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
              IF (k1/=5) THEN; imprnt=1; call derivs(t,yn,fn0); fnini(:,k1)=fn0; imprnt=0; fn=fn0;  ENDIF;
            ENDDO
        CASE(1)
        ! rozjezd, kdy AB formule nizsich radu pouziji jen na prvni ctyri kratsi kroky
        ! problem je, ze ten jeden krok ynp1=yn+h*fn jiz zpusobi onen zkoumany skok v Erot
            imprnt=1;
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
            imprnt=1; call derivs(t,yn,fn); fnini(:,1)=fn; imprnt=0;
            DO k1=1,20
                ynp4=yn; fn0=fn1; fn1=fn2; fn2=fn3; fn3=fn4; fn4=fn; 
                ynp5=ynp4+h*(1901./720.*fn4 - 1387./360.*fn3 + 109./30.*fn2 - 637./360.*fn1 + 251./720.*fn0)
                yn=ynp5; t=t+h;
                imprnt=1; 
                IF(k1/=20) call derivs(t,yn,fn);
                IF(k1==5) fnini(:,2)=fn
                IF(k1==10) fnini(:,3)=fn
                IF(k1==15) fnini(:,4)=fn
            ENDDO
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
      END SUBROUTINE
      
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
        SELECT CASE(icemode)
        CASE(0) !cap
            legendre = -rhoi*h/(4.-4.*cos(alfa)) * ((cos((n+1)*alfa)-cos((n+2)*alfa))/(n+3./2.) - (cos((n-1)*alfa)-cos(n*alfa))/(n-1./2.))
        CASE(1) !disc
            legendre = rhoi*h/2.*(-(5./2.*x**3-3./2.*x) + (x))
        CASE(2) !point
            legendre = Mpoint/(4.*pi*rmax**2.)*(2.*n+1)
        END SELECT
      END FUNCTION legendre

      FUNCTION absv(vektor)
      REAL  :: absv
      REAL , INTENT(IN) :: vektor(:)
      absv=sqrt(dot_product(vektor,vektor))
      END FUNCTION
      
      FUNCTION fluid_axis(J, Itpold, omega, dt, t)
      IMPLICIT NONE
      REAL, INTENT(IN) :: J(3,3), omega(3), Itpold(3,3,2), dt, t
      REAL :: Itp(3,3), Jrate(3,3), ItpPW(3,3), JratePW(3,3), Qmat(3,3), nomg(3), fluid_axisPW(3)
      REAL :: fluid_axis(3), Imax, dotprod, denom, eigenumber(3), eigevec(3,3), col, lon, A2B, CmA, CmB
      INTEGER :: idImax, k1, tmpimprnt
        Itp = J/sclJ
        nomg = omega/absv(omega)
        call eige_problem(Itp,eigenumber,eigevec,1)
        Imax=maxval(abs(eigenumber))
        DO k1=1,3
            IF(abs(eigenumber(k1))==Imax) idImax=k1
        ENDDO       

        ! Flipping the eigenvector: idImaxdotprod must be within the definition domain of acos <-1,1.> 
        dotprod = min(dot_product(eigevec(:,idImax),nomg), 1.)
        dotprod = max(dot_product(eigevec(:,idImax),nomg), -1.)
        IF(acos(dotprod)/pi*180. > 90.) THEN
            eigevec(:,idImax) = - eigevec(:,idImax)
        ENDIF

        SELECT CASE(quasi_fluid)        
        CASE(1)     ! Change of amplitude via conservation of angular momentum (normalized by sclJ, sclOmega)
            fluid_axis(1:3) = eigevec(1:3,idImax)
            fluid_axis(1:3) = fluid_axis(1:3) * absH0 / absv(matmul(Itp,fluid_axis))
        CASE(2)     ! Change of amplitude accounted via m3
            fluid_axis(1:2) = eigevec(1:2,idImax)
            fluid_axis(3)   = 1. - (Itp(3,3)-J0(3,3))/J0(3,3)            
        CASE(3)     ! Linearized Liouville equation (LLE)
            fluid_axis(1:2) = Itp(1:2,3) / (J0(3,3)-J0(1,1))            
            fluid_axis(3)   = 1. - (Itp(3,3)-J0(3,3)) / J0(3,3)
        CASE(4)     ! Extended LLE from Hu et al., 2017a,b
            tmpimprnt = imprnt;  imprnt = 0;
            Jrate = dJdt(Itp,Itpold,dt,t)
            imprnt = tmpimprnt
            ! Initial condition on dJ/dt, generally it can be non-zero
            IF(dt==0.) Jrate=0.
            IF (PW_Hu) THEN
                Qmat(1,1) = nomg(3) + nomg(2)**2/(1.+nomg(3))
                Qmat(2,2) = 1. - nomg(2)**2/(1.+nomg(3))
                Qmat(3,3) = nomg(3)
                Qmat(1,2) = - nomg(1)*nomg(2)/(1.+nomg(3))
                Qmat(1,3) = nomg(1)
                Qmat(2,3) = nomg(2)
                Qmat(2,1) = Qmat(1,2) 
                Qmat(3,1) = -nomg(1)
                Qmat(3,2) = -nomg(2)
                ItpPW = matmul(TRANSPOSE(Qmat), matmul(Itp,Qmat))
                JratePW = matmul(TRANSPOSE(Qmat), matmul(Jrate,Qmat))
            ELSE
                ItpPW = Itp; JratePW = Jrate;
            ENDIF           
            denom = sclOmega * (J0(3,3)-J0(1,1)) * (J0(3,3)-J0(2,2))
            CmA = J0(3,3)-J0(1,1)
            CmB = J0(3,3)-J0(2,2)
            A2B = 2.*J0(1,1)*J0(2,2)
            IF (Hu_Venus) THEN
                ! Equation 23 from Hu et al.,2017b
                fluid_axisPW(1)   = (CmB*ItpPW(1,3)*(sclOmega*dt)**2 + 2.*J0(2,2)*sclOmega*dt*ItpPW(2,3)) / A2B
                fluid_axisPW(2)   = (CmA*ItpPW(2,3)*(sclOmega*dt)**2 - 2.*J0(1,1)*sclOmega*dt*ItpPW(1,3)) / A2B
                fluid_axisPW(3)   = 1. - (ItpPW(3,3)-J0(3,3)) / J0(3,3)           
            ELSE 
                ! Equation 11 from Hu et al., 2017a
                fluid_axisPW(1)   = ItpPW(1,3) / CmA + J0(3,3)*JratePW(2,3)/denom
                fluid_axisPW(2)   = ItpPW(2,3) / CmB + J0(3,3)*JratePW(1,3)/denom
                fluid_axisPW(3)   = 1. - (ItpPW(3,3)-J0(3,3)) / J0(3,3)
            ENDIF
            IF (PW_Hu) THEN
                fluid_axis = matmul(Qmat, fluid_axisPW)
            ELSE
                fluid_axis = fluid_axisPW 
            ENDIF
        CASE DEFAULT
            stop 'unrecognized quasi_fluid option'
        END SELECT        
        
        col = acos(fluid_axis(3)/absv(fluid_axis))/pi*180
        lon = (atan(fluid_axis(2)/fluid_axis(1))/pi*180)
        !print *,'Colatitude, Longitude  ', colat, longit        
        
      END FUNCTION fluid_axis

      SUBROUTINE print_eige(A,label,wo)
      REAL, INTENT(IN) :: A(3,3)
      CHARACTER(*), INTENT(IN) :: label
      LOGICAL, INTENT(IN) :: wo
      REAL :: eigen(3),eigev(3,3),zunit(3)=(/0.,0.,1./)
          call eige_problem(A,eigen,eigev,1)
          print *,label
          print *,'eigenumbers: '
          print *,eigen(:)          
          print *,'eigenvectors: '
          print *,'1: ', eigev(:,1), 'colatitude: ', acos(dot_product(eigev(:,1),zunit))/pi*180.
          print *,'2: ', eigev(:,2), 'colatitude: ', acos(dot_product(eigev(:,2),zunit))/pi*180.
          print *,'3: ', eigev(:,3), 'colatitude: ', acos(dot_product(eigev(:,3),zunit))/pi*180.
          print *,'Input tensor: '
          print *,A(:,1)
          print *,A(:,2)
          print *,A(:,3)
          print *

          IF(wo) THEN
            open(94,file='run/principal.dat',Access='append')
            write(94,'(7e25.15,a,a)') 0.,0.,0.,eigev(1,1),eigev(2,1),eigev(3,1),eigen(1), '  ', label
            write(94,'(7e25.15,a,a)') 0.,0.,0.,eigev(1,2),eigev(2,2),eigev(3,2),eigen(2), '  ', label
            write(94,'(7e25.15,a,a)') 0.,0.,0.,eigev(1,3),eigev(2,3),eigev(3,3),eigen(3), '  ', label
            close(94)
          END IF
      END SUBROUTINE print_eige
      
      INCLUDE 'check_solver.f90'
      INCLUDE '3Dvisualization.f90'
      INCLUDE 'equations.f90'
      
      END MODULE
