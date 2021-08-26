! Numerical Recipe pro pridruzene Legendrovy funkce (double precision)
! nrbook kap. 6.8

	  MODULE nr
	  USE mConst
	  IMPLICIT NONE
	  CONTAINS
	  
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END SUBROUTINE

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=1000,TINY=1.0d-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
		print *,'singular matrix in ludcmp'
		crash=1
		return
	 endif
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE


      FUNCTION plgndr(l,m,x)
      INTEGER l,m
      REAL plgndr,x
      INTEGER i,ll
      REAL fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0)stop
     *'bad arguments in plgndr'
      pmm=1.d0
      if(m.gt.0) then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END FUNCTION
      
      FUNCTION dpythag(a,b)
      REAL a,b,dpythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END FUNCTION

            
      
      SUBROUTINE dsvdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=1000)
                      !CU    USES dpythag
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) !,dpythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=dpythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) stop 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=dpythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=dpythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=dpythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END SUBROUTINE
            
      SUBROUTINE dsvbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=1000)
      INTEGER i,j,jj
      REAL s,tmp(NMAX)
      do 12 j=1,n
        s=0.0d0
        if(w(j).ne.0.0d0)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.0d0
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END SUBROUTINE

      SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL d,a(np,mp),al(np,mpl),TINY
      PARAMETER (TINY=1.d-20)
      INTEGER i,j,k,l,mm
      REAL dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) stop 'bad args in bandec'
      l=m1
      do 13 i=1,m1
        do 11 j=m1+2-i,mm
          a(i,j-l)=a(i,j)
11      continue
        l=l-1
        do 12 j=mm-l,mm
          a(i,j)=0.d0
12      continue
13    continue
      d=1.d0
      l=m1
      do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
          if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
          endif
14      continue
        indx(k)=i
        if(dum.eq.0.d0) a(k,1)=TINY
        if(i.ne.k)then
          d=-d
          do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
15        continue
        endif
        do 17 i=k+1,l
          dum=a(i,1)/a(k,1)
          al(k,i-k)=dum
          do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
16        continue
          a(i,mm)=0.d0
17      continue
18    continue
      return
      END SUBROUTINE
            
      SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL a(np,mp),al(np,mpl),b(n)
      INTEGER i,k,l,mm
      REAL dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) stop 'bad args in banbks'
      l=m1
      do 12 k=1,n
        i=indx(k)
        if(i.ne.k)then
          dum=b(k)
          b(k)=b(i)
          b(i)=dum
        endif
        if(l.lt.n)l=l+1
        do 11 i=k+1,l
          b(i)=b(i)-al(k,i-k)*b(k)
11      continue
12    continue
      l=1
      do 14 i=n,1,-1
        dum=b(i)
        do 13 k=2,l
          dum=dum-a(i,k)*b(k+i-1)
13      continue
        b(i)=dum/a(i,1)
        if(l.lt.mm) l=l+1
14    continue
      return
      END SUBROUTINE
	  
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=3)
CU    USES derivs
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
     *,
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31
     *=3.d0/40.d0,
     *B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52
     *=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0
     */512.d0,
     *B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1
     *=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0
     */27648.d0,
     *DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0
     */14336.d0,
     *DC6=C6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END SUBROUTINE
	  
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=3)
CU    USES derivs,rkck
      INTEGER i
      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY
     *,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.d0)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x) then
		print *, 'stepsize underflow in rkqs'
		crash=1
		return
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END SUBROUTINE
      
      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=5)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      print *,'too many iterations in jacobi'
	  crash=1;
      return
      END SUBROUTINE
      
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv)
     *,SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=3,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     *REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
     *,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x) then
		print *,'step size underflow in bsstep'
		crash=1
		return
	 endif
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
		if (crash==1) return
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END SUBROUTINE
      
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=3)
      INTEGER i,n
      REAL h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END SUBROUTINE
      
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END SUBROUTINE
	  
	  
      SUBROUTINE stiff(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs,jac)
      INTEGER n,NMAX,MAXTRY
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
     *,SAFETY,GROW,
     *PGROW,SHRNK,PSHRNK,ERRCON,GAM,A21,A31,A32,A2X,A3X,C21,C31,C32,C41,
     *C42,C43,B1,B2,B3,B4,E1,E2,E3,E4,C1X,C2X,C3X,C4X
      EXTERNAL derivs,jac
      PARAMETER (NMAX=50,SAFETY=0.9d0,GROW=1.5d0,PGROW=-.25d0,SHRNK
     *=0.5d0,
     *PSHRNK=-1.d0/3.d0,ERRCON=.1296d0,MAXTRY=40)
      PARAMETER (GAM=1.d0/2.d0,A21=2.d0,A31=48.d0/25.d0,A32=6.d0/25.d0
     *,C21=-8.d0,
     *C31=372.d0/25.d0,C32=12.d0/5.d0,C41=-112.d0/125.d0,C42=-54.d0
     */125.d0,C43=-2.d0/5.d0,
     *B1=19.d0/9.d0,B2=1.d0/2.d0,B3=25.d0/108.d0,B4=125.d0/108.d0,E1
     *=17.d0/54.d0,E2=7.d0/36.d0,
     *E3=0.d0,E4=125.d0/108.d0,C1X=1.d0/2.d0,C2X=-3.d0/2.d0,C3X=121.d0
     */50.d0,C4X=29.d0/250.d0,
     *A2X=1.d0,A3X=3.d0/5.d0)
CU    USES derivs,jac,lubksb,ludcmp
      INTEGER i,j,jtry,indx(NMAX)
      REAL d,errmax,h,xsav,a(NMAX,NMAX),dfdx(NMAX),dfdy(NMAX
     *,NMAX),
     *dysav(NMAX),err(NMAX),g1(NMAX),g2(NMAX),g3(NMAX),g4(NMAX),
     *ysav(NMAX)
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      call jac(xsav,ysav,dfdx,dfdy,n)      
      h=htry
      do 23 jtry=1,MAXTRY
        do 13 i=1,n
          do 12 j=1,n
            a(i,j)=-dfdy(i,j)
12        continue
          a(i,i)=1.d0/(GAM*h)+a(i,i)
13      continue
        call ludcmp(a,n,NMAX,indx,d)
        do 14 i=1,n
          g1(i)=dysav(i)+h*C1X*dfdx(i)
14      continue
        call lubksb(a,n,NMAX,indx,g1)
        do 15 i=1,n
          y(i)=ysav(i)+A21*g1(i)
15      continue
        x=xsav+A2X*h
        call derivs(x,y,dydx)
        do 16 i=1,n
          g2(i)=dydx(i)+h*C2X*dfdx(i)+C21*g1(i)/h
16      continue
        call lubksb(a,n,NMAX,indx,g2)
        do 17 i=1,n
          y(i)=ysav(i)+A31*g1(i)+A32*g2(i)
17      continue
        x=xsav+A3X*h
        call derivs(x,y,dydx)
        do 18 i=1,n
          g3(i)=dydx(i)+h*C3X*dfdx(i)+(C31*g1(i)+C32*g2(i))/h
18      continue
        call lubksb(a,n,NMAX,indx,g3)
        do 19 i=1,n
          g4(i)=dydx(i)+h*C4X*dfdx(i)+(C41*g1(i)+C42*g2(i)+C43*g3(i))/h
19      continue
        call lubksb(a,n,NMAX,indx,g4)
        do 21 i=1,n
          y(i)=ysav(i)+B1*g1(i)+B2*g2(i)+B3*g3(i)+B4*g4(i)
          err(i)=E1*g1(i)+E2*g2(i)+E3*g3(i)+E4*g4(i)
21      continue
        x=xsav+h
        if(x.eq.xsav) then
		print *,'stepsize not significant in stiff'
		crash=1
		return
	 endif
        errmax=0.d0
        do 22 i=1,n
          errmax=max(errmax,abs(err(i)/yscal(i)))
22      continue
        errmax=errmax/eps
        if(errmax.le.1.d0)then
          hdid=h
          if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*errmax**PGROW
          else
            hnext=GROW*h
          endif
          return
        else
          hnext=SAFETY*h*errmax**PSHRNK
          h=sign(max(abs(hnext),SHRNK*abs(h)),h)
        endif
23    continue
      stop 'exceeded MAXTRY in stiff'
      END SUBROUTINE
      
      
      SUBROUTINE stifbs(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs
     *,jac)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv)
     *,SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      EXTERNAL derivs,jac
      PARAMETER (NMAX=50,KMAXX=7,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     *REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0)
CU    USES derivs,jac,simpr,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nvold,nseq(IMAX)
      REAL eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
     *,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),dfdx(NMAX),dfdy(NMAX,NMAX),err(KMAXX),
     *yerr(NMAX),ysav(NMAX),yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,nvold,xnew
      DATA first/.true./,epsold/-1.d0/,nvold/-1/
      DATA nseq /2,6,10,14,22,34,50,70/
      if(eps.ne.epsold.or.nv.ne.nvold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        nvold=nv
        a(1)=nv+a(1)
        do 14 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
14      continue
        do 15 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
15      continue
1       kmax=kopt
      endif
      h=htry
      do 16 i=1,nv
        ysav(i)=y(i)
16    continue
      call jac(x,y,dfdx,dfdy,nv,nmax)
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 18 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)stop 'stepsize underflow in stifbs'
        call simpr(ysav,dydx,dfdx,dfdy,nmax,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 17 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
17        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
18    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 19 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
19    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END SUBROUTINE
      
      SUBROUTINE simpr(y,dydx,dfdx,dfdy,nmax,n,xs,htot,nstep,yout,
     *derivs)
      INTEGER n,nmax,nstep,NMAXX
      REAL htot,xs,dfdx(n),dfdy(nmax,nmax),dydx(n),y(n),yout
     *(n)
      EXTERNAL derivs
      PARAMETER (NMAXX=50)
CU    USES derivs,lubksb,ludcmp
      INTEGER i,j,nn,indx(NMAXX)
      REAL d,h,x,a(NMAXX,NMAXX),del(NMAXX),ytemp(NMAXX)
      h=htot/nstep
      do 12 i=1,n
        do 11 j=1,n
          a(i,j)=-h*dfdy(i,j)
11      continue
        a(i,i)=a(i,i)+1.d0
12    continue
      call ludcmp(a,n,NMAXX,indx,d)
      do 13 i=1,n
        yout(i)=h*(dydx(i)+h*dfdx(i))
13    continue
      call lubksb(a,n,NMAXX,indx,yout)
      do 14 i=1,n
        del(i)=yout(i)
        ytemp(i)=y(i)+del(i)
14    continue
      x=xs+h
      call derivs(x,ytemp,yout)
      do 17 nn=2,nstep
        do 15 i=1,n
          yout(i)=h*yout(i)-del(i)
15      continue
        call lubksb(a,n,NMAXX,indx,yout)
        do 16 i=1,n
          del(i)=del(i)+2.d0*yout(i)
          ytemp(i)=ytemp(i)+del(i)
16      continue
        x=x+h
        call derivs(x,ytemp,yout)
17    continue
      do 18 i=1,n
        yout(i)=h*yout(i)-del(i)
18    continue
      call lubksb(a,n,NMAXX,indx,yout)
      do 19 i=1,n
        yout(i)=ytemp(i)+yout(i)
19    continue
      return
      END SUBROUTINE
	  
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END SUBROUTINE
      
      subroutine manexp(x, e) 
      integer e 
      real xx, x, m 
      if (x .lt. 0.) then 
         xx = -x 
         e = int(log10(xx)) 
         if (e .gt. 0) e = e + 1 
         m = - xx * 1.e1 ** (-e) 
      else if (x .gt. 0.) then 
         e = int(log10(x)) 
         if (e .gt. 0) e = e + 1       
         m = x * 1.e1 ** (-e) 
      else 
         e = 0. 
         m = 0. 
      endif 
      end subroutine
      
	  END MODULE
