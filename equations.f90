      SUBROUTINE BC_quasifreesurf(j,surf)     !skalovane konstantou sclBC
      INTEGER :: j      
      CHARACTER(5) :: surf
      IF(surf(1:5)=='outer') THEN
        !okrajova podminka BC 1
		AA(1,1)=+0.5*j*rho0(1)*g0(1)/(2.*j+1.)*sclBC
        AA(1,7)=AA(1,1)
        AA(1,2)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        AA(1,8)=AA(1,2)
        AA(1,3)=-sqrt(j/(3.*(2.*j+1.)))*sclBC
        AA(1,4)=sqrt((j-1.)/(2.*j-1.))*sclBC
        AA(1,5)=-sqrt((j+1.)*(2.*j+3)/(6.*(2.*j+1.)*(2.*j-1.)))*sclBC		
        !okrajova podminka BC 2 
        AA(2,1)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        AA(2,7)=AA(2,1)	
        AA(2,2)=+0.5*(j+1.)*rho0(1)*g0(1)/(2.*j+1.)*sclBC
        AA(2,8)=AA(2,2)
		AA(2,3)=sqrt((j+1.)/(3.*(2.*j+1.)))*sclBC
        AA(2,5)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3)))*sclBC
		AA(2,6)=-sqrt((j+2.)/(2.*j+3))*sclBC
      ELSEIF(surf(1:5)=='inner') THEN      
        !okrajova podminka BC 1  
        AA(6*N+1,6*N-5)=-0.5*j*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
        AA(6*N+1,6*N+1)=AA(6*N+1,6*N-5)
        AA(6*N+1,6*N-4)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        AA(6*N+1,6*N+2)=AA(6*N+1,6*N-4)
        AA(6*N+1,6*N-3)=-sqrt(j/(3.*(2.*j+1.)))*sclBC
        AA(6*N+1,6*N-2)=sqrt((j-1.)/(2.*j-1.))*sclBC
        AA(6*N+1,6*N-1)=-sqrt((j+1.)*(2.*j+3)/(6.*(2.*j+1.)*(2.*j-1.)))*sclBC
        !okrajova podminka BC 2 
        AA(6*N+2,6*N-5)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        AA(6*N+2,6*N+1)=AA(6*N+2,6*N-5)	
        AA(6*N+2,6*N-4)=-0.5*(j+1.)*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
        AA(6*N+2,6*N+2)=AA(6*N+2,6*N-4)
		AA(6*N+2,6*N-3)=sqrt((j+1.)/(3.*(2.*j+1.)))*sclBC
        AA(6*N+2,6*N-1)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3)))*sclBC
		AA(6*N+2,6*N)=-sqrt((j+2.)/(2.*j+3))*sclBC
      ELSE; stop 'unrecognized BC';
      ENDIF
      END SUBROUTINE    
    
      SUBROUTINE BC_noslip(surf)
      CHARACTER(5) :: surf
      IF(surf(1:5)=='outer') THEN
		AA(1,1)=1./2							!okrajova podminka BC 1 
		AA(1,7)=1./2
        AA(2,2)=1./2		                    !okrajova podminka BC 2							
		AA(2,8)=1./2
      ELSEIF(surf(1:5)=='inner') THEN      
   		AA(6*N+1,6*N+1)=1./2
		AA(6*N+1,6*N-5)=1./2
		AA(6*N+2,6*N+2)=1./2
		AA(6*N+2,6*N-4)=1./2  
      ELSE; stop 'unrecognized BC';
      ENDIF
      END SUBROUTINE

	  SUBROUTINE napln_A(j,A)					!naplni matici A(6x6) na kazdem rozhrani, s vyjimkou diagonalnich prvku
	  integer,intent(in) :: j									!poradi rovnic je: P.R. 1 - 2, R.K, Reo 1 - 3
	  real,intent(out) :: A(6,6)		
	  
			A(1,3)=(-1.)*(j+1.)*sqrt(j/(6.*j+3.))
			A(1,4)=(-1.)*(j-2.)*sqrt((j-1.)/(2.*j-1.))
			A(1,5)=(-1.)*(j+1.)*sqrt((j+1.)*(2.*j+3.)/((12*j-6.)*(2.*j+1.)))
			A(2,3)=(-1.)*(j)*sqrt((j+1.)/(6.*j+3.))
			A(2,5)=(-1.)*(j)*sqrt((j*(2.*j-1.))/((12*j+18.)*(2.*j+1.)))
			A(2,6)=(-1.)*(j+3.)*sqrt((j+2.)/(2.*j+3.))
			A(3,1)=(-1.)*(j-1.)*sqrt(j/(2.*j+1.))*sclRK
			A(3,2)=(-1.)*(j+2.)*sqrt((j+1.)/(2.*j+1.))*sclRK
			A(4,1)=2.*j*Wgnr_smbl(j,j-2,j-1)*sqrt(5*j-5.)*sclREO				!viskozitou jsme reorovnice podelili
			A(5,1)=2.*(j-1)*Wgnr_smbl(j,j,j-1)*sqrt(5.*j)*sclREO
			A(5,2)=2.*(j+2)*Wgnr_smbl(j,j,j+1)*sqrt(5*j+5.)*sclREO
			A(6,2)=2.*(j+1)*Wgnr_smbl(j,j+2,j+1)*sqrt(5*j+10.)*sclREO
	
	  END SUBROUTINE

	  SUBROUTINE napln_B(j,A)					!naplni matici B(6x6) na kazdem rozhrani, s vyjimkou diagonalnich prvku
	  integer,intent(in) :: j									!poradi rovnic je: P.R. 1 - 2, R.K, Reo 1 - 3
	  real,intent(out) :: A(6,6)
		
			A(1,3)=(-1.)*sqrt(j/(6.*j+3.))
			A(1,4)=sqrt((j-1.)/(2.*j-1.))
			A(1,5)=(-1.)*sqrt((j+1.)*(2.*j+3.)/((12*j-6.)*(2.*j+1.)))
			A(2,3)=sqrt((j+1.)/(6.*j+3.))
			A(2,5)=sqrt((j*(2.*j-1.))/((12*j+18.)*(2.*j+1.)))
			A(2,6)=(-1.)*sqrt((j+2.)/(2.*j+3.))
			A(3,1)=sqrt(j/(2.*j+1.))*sclRK
			A(3,2)=(-1.)*sqrt((j+1.)/(2.*j+1.))*sclRK
			A(4,1)=2.*Wgnr_smbl(j,j-2,j-1)*sqrt(5*j-5.)*sclREO			!viskozitou jsme reorovnice podelili
			A(5,1)=-2.*Wgnr_smbl(j,j,j-1)*sqrt(5.*j)*sclREO
			A(5,2)=2.*Wgnr_smbl(j,j,j+1)*sqrt(5*j+5.)*sclREO
			A(6,2)=-2.*Wgnr_smbl(j,j+2,j+1)*sqrt(5*j+10.)*sclREO
	
	  END SUBROUTINE

	  SUBROUTINE rovkon(i)							!rovnice kontinuity - v AA je polozena na prvni radek za BC
	  integer i
			AA(6*i-3,6*i-5)=A(3,1)/(2.*redge(i))+B(3,1)/dr
			AA(6*i-3,6*i+1)=A(3,1)/(2.*redge(i))-B(3,1)/dr
			AA(6*i-3,6*i-4)=A(3,2)/(2.*redge(i))+B(3,2)/dr
			AA(6*i-3,6*i+2)=A(3,2)/(2.*redge(i))-B(3,2)/dr
	  END SUBROUTINE

	  SUBROUTINE reovztah(i)						!Reologicke vztahy - v AA jsou polozeny na trech radcich za R.K.
	  integer i			
			AA(6*i-2,6*i-5)=A(4,1)/(2.*redge(i))+B(4,1)/dr		!Reo vztah 1
			AA(6*i-2,6*i+1)=A(4,1)/(2.*redge(i))-B(4,1)/dr
			AA(6*i-2,6*i-2)=A(4,4)

			AA(6*i-1,6*i-5)=A(5,1)/(2.*redge(i))+B(5,1)/dr		!Reo vztah 2
			AA(6*i-1,6*i+1)=A(5,1)/(2.*redge(i))-B(5,1)/dr
			AA(6*i-1,6*i-4)=A(5,2)/(2.*redge(i))+B(5,2)/dr		
			AA(6*i-1,6*i+2)=A(5,2)/(2.*redge(i))-B(5,2)/dr
			AA(6*i-1,6*i-1)=A(5,5)

			AA(6*i,6*i-4)=A(6,2)/(2.*redge(i))+B(6,2)/dr        !Reo vztah 3
			AA(6*i,6*i+2)=A(6,2)/(2.*redge(i))-B(6,2)/dr
			AA(6*i,6*i)=A(6,6)
	  END SUBROUTINE

	  SUBROUTINE pohrov(i)
	  integer i
			AA(6*i+1,6*i-3)=A(1,3)/(2.*rcent(i+1))+B(1,3)/dr		! P.R. 1
			AA(6*i+1,6*i+3)=A(1,3)/(2.*rcent(i+1))-B(1,3)/dr
			AA(6*i+1,6*i-2)=A(1,4)/(2.*rcent(i+1))+B(1,4)/dr
			AA(6*i+1,6*i+4)=A(1,4)/(2.*rcent(i+1))-B(1,4)/dr
			AA(6*i+1,6*i-1)=A(1,5)/(2.*rcent(i+1))+B(1,5)/dr
			AA(6*i+1,6*i+5)=A(1,5)/(2.*rcent(i+1))-B(1,5)/dr
            AA(6*i+1,6*i+1)=A(1,1); AA(6*i+1,6*i+2)=A(1,2);
            
			AA(6*i+2,6*i-3)=A(2,3)/(2.*rcent(i+1))+B(2,3)/dr		! P.R. 2
			AA(6*i+2,6*i+3)=A(2,3)/(2.*rcent(i+1))-B(2,3)/dr
			AA(6*i+2,6*i-1)=A(2,5)/(2.*rcent(i+1))+B(2,5)/dr
			AA(6*i+2,6*i+5)=A(2,5)/(2.*rcent(i+1))-B(2,5)/dr
			AA(6*i+2,6*i)=A(2,6)/(2.*rcent(i+1))+B(2,6)/dr
			AA(6*i+2,6*i+6)=A(2,6)/(2.*rcent(i+1))-B(2,6)/dr
            AA(6*i+2,6*i+1)=A(2,1); AA(6*i+2,6*i+2)=A(2,2);
      END SUBROUTINE      
      
      SUBROUTINE rovkonab(i)							!rovnice kontinuity - v AA je polozena na prvni radek za BC
	  integer i,f,i1
            f=kl+ku+1;  i1=6*i-3;
			abini(f+i1-(6*i-5),6*i-5)=A(3,1)/(2.*redge(i))+B(3,1)/dr
			abini(f+i1-(6*i+1),6*i+1)=A(3,1)/(2.*redge(i))-B(3,1)/dr
			abini(f+i1-(6*i-4),6*i-4)=A(3,2)/(2.*redge(i))+B(3,2)/dr
			abini(f+i1-(6*i+2),6*i+2)=A(3,2)/(2.*redge(i))-B(3,2)/dr
	  END SUBROUTINE

	  SUBROUTINE reovztahab(i)						!Reologicke vztahy - v AA jsou polozeny na trech radcich za R.K.
	  integer i,f,i1
			f=kl+ku+1;  i1=6*i-2;
			abini(f+i1-(6*i-5),6*i-5)=A(4,1)/(2.*redge(i))+B(4,1)/dr		!Reo vztah 1
			abini(f+i1-(6*i+1),6*i+1)=A(4,1)/(2.*redge(i))-B(4,1)/dr
			abini(f+i1-(6*i-2),6*i-2)=A(4,4)
            i1=6*i-1;
			abini(f+i1-(6*i-5),6*i-5)=A(5,1)/(2.*redge(i))+B(5,1)/dr		!Reo vztah 2
			abini(f+i1-(6*i+1),6*i+1)=A(5,1)/(2.*redge(i))-B(5,1)/dr
			abini(f+i1-(6*i-4),6*i-4)=A(5,2)/(2.*redge(i))+B(5,2)/dr		
			abini(f+i1-(6*i+2),6*i+2)=A(5,2)/(2.*redge(i))-B(5,2)/dr
			abini(f+i1-(6*i-1),6*i-1)=A(5,5)
            i1=6*i;
			abini(f+i1-(6*i-4),6*i-4)=A(6,2)/(2.*redge(i))+B(6,2)/dr        !Reo vztah 3
			abini(f+i1-(6*i+2),6*i+2)=A(6,2)/(2.*redge(i))-B(6,2)/dr
			abini(f+i1-(6*i),6*i)=A(6,6)
	  END SUBROUTINE

	  SUBROUTINE pohrovab(i)
	  integer i,f,i1
            f=kl+ku+1;  i1=6*i+1;
			abini(f+i1-(6*i-3),6*i-3)=A(1,3)/(2.*rcent(i+1))+B(1,3)/dr		! P.R. 1
			abini(f+i1-(6*i+3),6*i+3)=A(1,3)/(2.*rcent(i+1))-B(1,3)/dr
			abini(f+i1-(6*i-2),6*i-2)=A(1,4)/(2.*rcent(i+1))+B(1,4)/dr
			abini(f+i1-(6*i+4),6*i+4)=A(1,4)/(2.*rcent(i+1))-B(1,4)/dr
			abini(f+i1-(6*i-1),6*i-1)=A(1,5)/(2.*rcent(i+1))+B(1,5)/dr
			abini(f+i1-(6*i+5),6*i+5)=A(1,5)/(2.*rcent(i+1))-B(1,5)/dr
            abini(f+i1-(6*i+1),6*i+1)=A(1,1); 
            abini(f+i1-(6*i+2),6*i+2)=A(1,2);
            i1=6*i+2;
			abini(f+i1-(6*i-3),6*i-3)=A(2,3)/(2.*rcent(i+1))+B(2,3)/dr		! P.R. 2
			abini(f+i1-(6*i+3),6*i+3)=A(2,3)/(2.*rcent(i+1))-B(2,3)/dr
			abini(f+i1-(6*i-1),6*i-1)=A(2,5)/(2.*rcent(i+1))+B(2,5)/dr
			abini(f+i1-(6*i+5),6*i+5)=A(2,5)/(2.*rcent(i+1))-B(2,5)/dr
			abini(f+i1-(6*i),6*i)=A(2,6)/(2.*rcent(i+1))+B(2,6)/dr
			abini(f+i1-(6*i+6),6*i+6)=A(2,6)/(2.*rcent(i+1))-B(2,6)/dr
            abini(f+i1-(6*i+1),6*i+1)=A(2,1); 
            abini(f+i1-(6*i+2),6*i+2)=A(2,2);
    END SUBROUTINE
    
    SUBROUTINE BCab_quasifreesurf(j,surf)     !skalovane konstantou sclBC
      INTEGER :: j,f,i1            
      CHARACTER(5) :: surf
      f=kl+ku+1;
      IF(surf(1:5)=='outer') THEN
        !okrajova podminka BC 1
        i1=1      
		abini(f+i1-1,1)=+0.5*j*rho0(1)*g0(1)/(2.*j+1.)*sclBC
        abini(f+i1-7,7)=+0.5*j*rho0(1)*g0(1)/(2.*j+1.)*sclBC
        abini(f+i1-2,2)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        abini(f+i1-8,8)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        abini(f+i1-3,3)=-sqrt(j/(3.*(2.*j+1.)))*sclBC
        abini(f+i1-4,4)=sqrt((j-1.)/(2.*j-1.))*sclBC
        abini(f+i1-5,5)=-sqrt((j+1.)*(2.*j+3)/(6.*(2.*j+1.)*(2.*j-1.)))*sclBC		
        !okrajova podminka BC 2 
        i1=2
        abini(f+i1-1,1)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        abini(f+i1-7,7)=-0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(1)*g0(1)*sclBC
        abini(f+i1-2,2)=+0.5*(j+1.)*rho0(1)*g0(1)/(2.*j+1.)*sclBC
        abini(f+i1-8,8)=+0.5*(j+1.)*rho0(1)*g0(1)/(2.*j+1.)*sclBC
		abini(f+i1-3,3)=sqrt((j+1.)/(3.*(2.*j+1.)))*sclBC
        abini(f+i1-5,5)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3)))*sclBC
		abini(f+i1-6,6)=-sqrt((j+2.)/(2.*j+3))*sclBC
      ELSEIF(surf(1:5)=='inner') THEN      
        !okrajova podminka BC 1  
        i1=6*N+1
        abini(f+i1-(6*N-5),6*N-5)=-0.5*j*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
        abini(f+i1-(6*N+1),6*N+1)=-0.5*j*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
        abini(f+i1-(6*N-4),6*N-4)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        abini(f+i1-(6*N+2),6*N+2)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        abini(f+i1-(6*N-3),6*N-3)=-sqrt(j/(3.*(2.*j+1.)))*sclBC
        abini(f+i1-(6*N-2),6*N-2)=sqrt((j-1.)/(2.*j-1.))*sclBC
        abini(f+i1-(6*N-1),6*N-1)=-sqrt((j+1.)*(2.*j+3)/(6.*(2.*j+1.)*(2.*j-1.)))*sclBC
        !okrajova podminka BC 2 
        i1=6*N+2
        abini(f+i1-(6*N-5),6*N-5)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        abini(f+i1-(6*N+1),6*N+1)=+0.5*sqrt(j*(j+1.))/(2.*j+1.)*rho0(N+1)*g0(N+1)*sclBC
        abini(f+i1-(6*N-4),6*N-4)=-0.5*(j+1.)*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
        abini(f+i1-(6*N+2),6*N+2)=-0.5*(j+1.)*rho0(N+1)*g0(N+1)/(2.*j+1.)*sclBC
		abini(f+i1-(6*N-3),6*N-3)=sqrt((j+1.)/(3.*(2.*j+1.)))*sclBC
        abini(f+i1-(6*N-1),6*N-1)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3)))*sclBC
		abini(f+i1-(6*N),6*N)=-sqrt((j+2.)/(2.*j+3))*sclBC
      ELSE; stop 'unrecognized BC';
      ENDIF
    END SUBROUTINE    
    
    SUBROUTINE BCab_noslip(surf)
      INTEGER :: f
      CHARACTER(5) :: surf
      f=kl+ku+1;
      IF(surf(1:5)=='outer') THEN
		abini(f+1-1,1)=1./2							!okrajova podminka BC 1 
		abini(f+1-7,7)=1./2
        abini(f+2-2,2)=1./2		                    !okrajova podminka BC 2							
		abini(f+2-8,8)=1./2
      ELSEIF(surf(1:5)=='inner') THEN              
   		abini(f+6*N+1-(6*N+1),6*N+1)=1./2
		abini(f+6*N+1-(6*N-5),6*N-5)=1./2
		abini(f+6*N+2-(6*N+2),6*N+2)=1./2
		abini(f+6*N+2-(6*N-4),6*N-4)=1./2  
      ELSE; stop 'unrecognized BC';
      ENDIF
    END SUBROUTINE

	  FUNCTION Wgnr_smbl(a,b,c)							!vypocet Wignerovych symbolu, obsazenych v maticich A,B
	  integer a,b,c
	  real Wgnr_smbl
		if (a==0) then
			Wgnr_smbl=(merge(1.,0.,c==1))*(merge(1.,0.,b==2))/(3.)
		else
			if ((b==a-2).and.(c==a-1)) then
				Wgnr_smbl=1./(sqrt(5.*(2.*a-1)))
			elseif ((b==a-1).and.(c==a)) then
				Wgnr_smbl=-sqrt((a-1.)/(10.*(a)*(2.*a+1)))
			elseif ((b==a).and.(c==a-1)) then
				Wgnr_smbl=sqrt(((a+1.)*(2.*a+3))/(30.*(a)*(2.*a-1)*(2.*a+1)))
			elseif ((b==a).and.(c==a+1)) then
				Wgnr_smbl=sqrt(((a)*(2.*a-1))/(30.*(2.*a+3)*(a+1)*(2.*a+1)))
			elseif ((b==a+1).and.(c==a)) then
				Wgnr_smbl=-sqrt((a+2.)/(10.*(a+1)*(2.*a+1)))
			elseif ((b==a+2).and.(c==a+1)) then
				Wgnr_smbl=1./sqrt(5.*(2.*a+3))
			else
				print *,'nedefinovany Wigneruv symbol'
				stop
			endif
		endif
	  END FUNCTION
