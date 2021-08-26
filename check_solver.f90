	  SUBROUTINE test_presnosti(y)
	  real :: tiny,abserr,relerr
      COMPLEX :: y(NN)
	  integer :: k1
		tiny=1.e-5				!pokud projde, pak soustava je zrejme uspesne resena s presnosti na tiny			
		do k1=1,6*N+2			
			if (abs(ycheck(k1)-dot_product(AAcheck(k1,:),y(:)))>tiny) then
                abserr=abs(ycheck(k1)-dot_product(AAcheck(k1,:),y(:)))
                relerr=abs(ycheck(k1)-dot_product(AAcheck(k1,:),y(:)))/abs(dot_product(y,y)/NN)
                print *,'Nepresne reseno - hloubka, stupen a rad:',k1,j
                print *,'abs. err: ',abserr,'rel. err:',relerr	
                stop
            endif
		enddo
	  END SUBROUTINE

	  SUBROUTINE zpetna_kontrola(y)
	  COMPLEX :: y(NN)
	  integer k1
	  open (20,file='rozhrani.dat')
	  do k1=1,N
		write(20,'(3f25.10)') redge(k1),2*shear(k1)*sqrt(5.*(j-1.))*Wgnr_smbl(j,j-2,j-1)*((y(6*k1-5)-y(6*k1+1))/dr &
		& + j*(y(6*k1-5)+y(6*k1+1))/redge(k1)/2)-y(6*k1-2)
		write(20,'(3f25.10)') redge(k1),2*shear(k1)*((-sqrt(5.*j))*Wgnr_smbl(j,j,j-1)*((y(6*k1-5)-y(6*k1+1))/dr-&
		&(j-1)*(y(6*k1-5)+y(6*k1+1))/redge(k1)/2)+(sqrt(5.*(j+1)))*Wgnr_smbl(j,j,j+1)*((y(6*k1-4)-y(6*k1+2))/dr+(j+2)*(y(6*k1-4)+&
		&y(6*k1+2))/redge(k1)/2))-y(6*k1-1)
		write(20,'(3f25.10)') redge(k1),(-2.)*shear(k1)*sqrt(5.*(j+2.))*Wgnr_smbl(j,j+2,j+1)*((y(6*k1-4)-y(6*k1+2))/dr &
		&-(j+1)*(y(6*k1-4)+y(6*k1+2))/redge(k1)/2)-y(6*k1)
		write(20,'(3f25.10)') redge(k1),sqrt(j*1.)*((y(6*k1-5)-y(6*k1+1))/dr-(y(6*k1-5)+y(6*k1+1))*(j-1)/redge(k1)/2)-&
		&sqrt(j+1.)*((y(6*k1-4)-y(6*k1+2))/dr+(y(6*k1-4)+y(6*k1+2))*(j+2)/redge(k1)/2)		
	  enddo
	  close (20)
	  open (21,file='stredy.dat')
	  do k1=2,N
		write (21,'(3f18.10)') rcent(k1),sqrt((j+1)/(6.*j+3))*((y(6*k1-9)-y(6*k1-3))/dr-(y(6*k1-9)+y(6*k1-3))*j/rcent(k1)/2)+&
		&sqrt((j*(2*j-1))/(6.*(2*j+3)*(2*j+1)))*((y(6*k1-7)-y(6*k1-1))/dr-(y(6*k1-7)+y(6*k1-1))*j/rcent(k1)/2)-&
		&sqrt((j+2)/(2.*j+3))*((y(6*k1-6)-y(6*k1))/dr+(y(6*k1-6)+y(6*k1))*(j+3)/rcent(k1)/2)
		write (21,'(3f18.10)') rcent(k1),-sqrt(j/(6.*j+3))*((y(6*k1-9)-y(6*k1-3))/dr+(y(6*k1-9)+y(6*k1-3))*(j+1)/rcent(k1)/2)+&
		&sqrt((j-1)/(2.*j-1))*((y(6*k1-8)-y(6*k1-2))/dr-(y(6*k1-8)+y(6*k1-2))*(j-2)/rcent(k1)/2)-&
		&sqrt((j+1.)*(2.*j+3)/(6.*(2.*j-1)*(2*j+1)))*((y(6*k1-7)-y(6*k1-1))/dr+(y(6*k1-7)+y(6*k1-1))*(j+1)/rcent(k1)/2)	
	  enddo
	  close (21)
      END SUBROUTINE
      
      SUBROUTINE printAA(y)
      COMPLEX :: y(NN)
      INTEGER :: k1
        open (25,file='AA.dat'); open (26,file='y.dat');
        do k1=1,NN ;
            write (25,'(1000e15.4)') AAcheck(k1,:); 
        enddo; 
        write (26,'(6e15.4)') real(y) 
        close(25); close(26);
      END SUBROUTINE
