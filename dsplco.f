c-----------------------------------------------------------
c  Double Precision version of splco.f
c  Routine from Numerical Recipes
c  x(n): array of n abszissae
c  y(n): array of n function values
c  n:    number of samples
c  yp1 : value of first derivative at x(1)
c  yp2 : value of first derivative at x(n)
c  If yp1 or ypn are set to 1.e30 the routine is
c  signaled to compute a natural spline with zero
c  second derivative at that boundary.
c  Returns the second derivatives y2(n).
c  Evaluation is done by the formula
c  y=a*y(j)+b*y(j+1)+c*y2(j)+d*y2(j+1)
c  with a=(x(j+1)-x)/h, b=1-a
c  c=h**2*(a**3-a)/6 and d=h**2*(b**3-b)/6
c  with h=x(j+1)-x(j)
c------------------------------------------------------------
      SUBROUTINE dsplco(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))
     1        /(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.
