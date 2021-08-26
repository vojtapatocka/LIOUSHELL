	  SUBROUTINE vec_sphharm_to_space(theta,phi,vjmin1,vjplus1,j,m,vx,vy,vz)
	  integer,intent(IN) :: j,m
	  complex :: vr1,vth1,vphi1,vr2,vth2,vphi2                   !koeficient pro v_jm_j-1, koeficient pro v_jm_j+1	 
      real :: e=10.**(-8)		
	  real, intent(OUT) :: vx,vy,vz		
	  real,intent(IN) :: theta,phi
      complex,intent(IN) :: vjmin1,vjplus1
	  	  
	  vr1=(j*vjmin1/sqrt((j)*(2.*j+1.)))*sphHarm(theta,phi,j,m)		!aplikace definicniho vztahu A.14
	  vth1=(sphHarm(theta+e,phi,j,m)-sphHarm(theta-e,phi,j,m))/(2*e)*(vjmin1/sqrt((j)*(2.*j+1.)))
	  vphi1=(sphHarm(theta,phi+e,j,m)-sphHarm(theta,phi-e,j,m))/(2*e*sin(theta))*(vjmin1/sqrt((j)*(2.*j+1.)))
	  vr2=(-(j+1.)*vjplus1/sqrt((j+1.)*(2.*j+1.)))*sphHarm(theta,phi,j,m)
	  vth2=(sphHarm(theta+e,phi,j,m)-sphHarm(theta-e,phi,j,m))/(2*e)*(vjplus1/sqrt((j+1.)*(2.*j+1.)))
	  vphi2=(sphHarm(theta,phi+e,j,m)-sphHarm(theta,phi-e,j,m))/(2*e*sin(theta))*(vjplus1/sqrt((j+1.)*(2.*j+1.)))
      vphi1=merge(vphi1,(0.,0.),abs(vphi1)>=0.); vphi2=merge(vphi2,(0.,0.),abs(vphi2)>=0.);
      vth1=merge(vth1,(0.,0.),abs(vth1)>=0.); vth2=merge(vth2,(0.,0.),abs(vth2)>=0.);
      vr1=merge(vr1,(0.,0.),abs(vr1)>=0.); vr2=merge(vr2,(0.,0.),abs(vr2)>=0.);
	  vx=cos(theta)*cos(phi)*(real(vth1+vth2))+sin(theta)*cos(phi)*real(vr1+vr2)-sin(phi)*real(vphi1+vphi2)
	  vy=cos(phi)*real(vphi1+vphi2)+sin(theta)*sin(phi)*real(vr1+vr2)+cos(theta)*sin(phi)*real(vth1+vth2)
	  vz=-sin(theta)*real(vth1+vth2)+cos(theta)*real(vr1+vr2)
	  
      END SUBROUTINE
      
      SUBROUTINE zapis3D(NR,Ntheta,Nphi,nid)
      real :: theta,phi,x1,x2,x3
      integer :: k1,k2,k3,nid
      integer, intent(IN) :: NR,Ntheta,Nphi
      
      select case(modd)
		case(1)
			open (2,file='para_skalar'//trim(str(nid))//'.vtk')
			call para_hlavicka(NR,Ntheta,Nphi,2)	
		case(2)
			open (2,file='para_vek.vtk')
			call para_hlavicka(NR,Ntheta,Nphi,2)	
		case(3)
			open (42,file='vrad'//trim(str(nid))//'.vtk')
			open (46,file='rych_pole'//trim(str(nid))//'.vtk')
			call para_hlavicka(NR,Ntheta,Nphi,42)
			call para_hlavicka(NR,Ntheta,Nphi,46)
		end select					 

		do k1=2,NR
			do k2=0,Ntheta					
			theta=pi*k2/Ntheta
				do k3=0,Nphi 
					phi=2*pi*k3/Nphi
					x1=rcent(k1)*sin(theta)*cos(phi)
					x2=rcent(k1)*sin(theta)*sin(phi)
					x3=rcent(k1)*cos(theta)
					select case(modd)
					case(1,2)
						write(2,*) x1,x2,x3					
					case(3)
						write(42,*) x1,x2,x3; 
						write(46,*) x1,x2,x3; 
					end select
				enddo
			enddo
		enddo	

		select case(modd)
		case(1)
			write (2,'(A,1I8)') 'POINT_DATA ',(Nphi+1)*(Ntheta+1)*(NR-1)
			write (2,'(A)') 'SCALARS skalar float'
            write (2,'(A)') 'LOOKUP_TABLE default'
        case(2)
            write (2,'(A,1I8)') 'POINT_DATA ',(Nphi+1)*(Ntheta+1)*(NR-1)
			write (2,'(A)') 'VECTORS rychlost float'
		case(3)
			write (42,'(A,1I8)') 'POINT_DATA ',(Nphi+1)*(Ntheta+1)*(NR-1)
			write (42,'(A)') 'SCALARS vrad float'
            write (42,'(A)') 'LOOKUP_TABLE default'
			write (46,'(A,1I8)') 'POINT_DATA ',(Nphi+1)*(Ntheta+1)*(NR-1)
			write (46,'(A)') 'VECTORS rychlost float'
        end select
        
        do k1=2,NR
			do k2=0,Ntheta					!bereme jen realnou cast harmonik
			theta=pi*k2/Ntheta
				do k3=0,Nphi 
					phi=2*pi*k3/Nphi					
					select case (modd)
					case (1)				
						write(2,*) D3vrad(k1,k2,k3)	
					case (2)							
						write(2,*) D3vx(k1,k2,k3),D3vy(k1,k2,k3),D3vz(k1,k2,k3)
                    case (3)
                        write(42,*) D3vrad(k1,k2,k3)	
						write(46,*) D3vx(k1,k2,k3),D3vy(k1,k2,k3),D3vz(k1,k2,k3)
					end select
				enddo
			enddo
        enddo	
        
        D3vx=0. ; D3vy=0. ; D3vz=0. ; D3vrad=0. ;
        select case (modd)
		case (1,2)
			close(2)		
		case(3)
			close(42); close(45);   close(46)
		end select
      
      
      END SUBROUTINE
      
	  SUBROUTINE rez_radialni_slozky_v(y,Ntheta,Nphi,NR,j,m,modd)
	  complex :: y(6*N+2),v_jmin1(N+1),v_jplus1(N+1)
	  real :: theta,phi,x1,x2,x3,vx,vy,vz
	  integer, intent(IN) :: NR,Ntheta,Nphi,j,m,modd
	  integer :: k1,k2,k3	  		
		v_jmin1=y(1::6)
		v_jplus1=y(2::6)		
		do k1=2,NR
			do k2=0,Ntheta					
			theta=(pi*k2)/Ntheta
				do k3=0,Nphi 
					phi=(2*pi*k3)/Nphi					
					select case (modd)
                    case(1)
                        !D3vrad(k1,k2,k3)=D3vrad(k1,k2,k3)+real((sqrt(j/(2.*j+1))*v_jmin1(k1)-sqrt((j+1)/(2.*j+1))*v_jplus1(k1))*sphHarm(theta,phi,j,m))
                        !D3vrad(k1,k2,k3)=D3vrad(k1,k2,k3) + rho0(k1)/3
                        D3vrad(k1,k2,k3)=D3vrad(k1,k2,k3)+real(extload(m)*sphHarm(theta,phi,j,m))
					case(2)				
						call vec_sphharm_to_space(theta,phi,y(6*k1-5),y(6*k1-4),j,m,vx,vy,vz)
						D3vx(k1,k2,k3)=D3vx(k1,k2,k3)+vx ; D3vy(k1,k2,k3)=D3vy(k1,k2,k3)+vy ; D3vz(k1,k2,k3)=D3vz(k1,k2,k3)+vz ;
					case(3)
						D3vrad(k1,k2,k3)=D3vrad(k1,k2,k3) + real((sqrt(j/(2.*j+1))*v_jmin1(k1)-sqrt((j+1)/(2.*j+1))*v_jplus1(k1))*sphHarm(theta,phi,j,m))
						call vec_sphharm_to_space(theta,phi,y(6*k1-5),y(6*k1-4),j,m,vx,vy,vz)
						D3vx(k1,k2,k3)=D3vx(k1,k2,k3)+vx ; D3vy(k1,k2,k3)=D3vy(k1,k2,k3)+vy ; D3vz(k1,k2,k3)=D3vz(k1,k2,k3)+vz ;
					end select
				enddo
			enddo
		enddo			
      END SUBROUTINE
	  
      SUBROUTINE para_hlavicka(NR,Ntheta,Nphi,idfile)
      INTEGER, INTENT(IN) :: NR,Ntheta,Nphi,idfile
        write (idfile,'(A)') '# vtk DataFile Version 3.0'
        write (idfile,'(A)')
        write (idfile,'(A)') 'ASCII'
        write (idfile,'(A)') 'DATASET STRUCTURED_GRID'
        write (idfile,'(A,3I5)') 'DIMENSIONS ',Nphi+1,Ntheta+1,NR-1
        write (idfile,'(A,1I8,A)') 'POINTS ',(Nphi+1)*(Ntheta+1)*(NR-1),' float'        
      END SUBROUTINE
      
	  SUBROUTINE amira_hlavicka(NR,Ntheta,Nphi,form,file)
	  integer, intent(in) :: NR,Ntheta,Nphi,form,file
		write (file,'(A)') '# AmiraMesh ASCII 1.0'
		write (file,'(A)')
		write (file,'(A,3I5)') 'define Lattice',Nphi+1,Ntheta+1,NR-1
		write (file,'(A)')
		write (file,'(A)') 'Parameters {'
		write (file,'(A)') '    CoordType "curvilinear"'
		write (file,'(A)') '}'
		write (file,'(A)')
		write (file,'(A)') 'Lattice { float[3] Coordinates } = @1'
		if (form==1) write (file,'(A)') 'Lattice { float ScalarField } = @2'
		if (form==2) write (file,'(A)') 'Lattice { float[3] VectorField } = @2'
		write (file,'(A)')
		write (file,'(A)') '@1'		
      END SUBROUTINE

	  FUNCTION SphHarm(theta,phi,j,m)           !pozor, je trik s m<0 udelany spravne?????????????
	  REAL :: theta,phi,plm,faktor
      COMPLEX :: SphHarm
	  INTEGER :: j,m        
        if (m/=0) then ; faktor=2. ; else ; faktor =1. ; endif        !umele vyscitani i pres harmoniky s m<0
		plm=plgndr(j,m,cos(theta))
		SphHarm=faktor*plm*sqrt((2*j+1)*(fac(j-m))/(4*pi*fac(j+m)))*cos(m*phi)*(1.,0.)
        SphHarm=SphHarm+faktor*plm*sqrt((2*j+1)*(fac(j-m))/(4*pi*fac(j+m)))*sin(m*phi)*(0.,1.)        
      END FUNCTION
        
	  RECURSIVE FUNCTION fac(n) RESULT (result) ! rekurzivni funkce s povinnou navratovou hodnotou
	  REAL result
	  INTEGER n
	  result=1.
	  if (n>1) result=n*fac(n-1)        ! rekurzivni volani
      END FUNCTION  
      
    CHARACTER(len=20) FUNCTION str(k)
!   "Convert an integer to string."
    INTEGER, INTENT(IN) :: k
    write (str, '(I3.3)') k
    str = adjustl(str)
    END FUNCTION str
