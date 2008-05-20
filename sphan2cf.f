c ----------------------------------------------------------------------------
c   ($Id:$)
c ============================================================================
c
c  Computes ratio of (complex) spherical Hankel-functions
c  of the second kind using a continued fraction
c 
c    z(n):= x*h(n+1)/h(n)
c
c    z(n)=n+1+zi*x+zi*CF
c    CF = a1/(b1 +) a2/(b2 + )......
c    a(j)=(j-1/2)**2-(n+1/2)**2, b(j)=2*(x+j*zi)
c
c  The continued fraction
c  is evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
c  al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
c  p. 169 for further information.
c--------------------------------------------------------------------
	subroutine sphan2cf(zx,el,zf,iprint)
	include 'pi.h'
	double precision tiny,el,eps,elph2
	double complex zi
	integer maxit,j,iprint
	parameter(maxit=10000,tiny=1.d-30,eps=1.d-6,zi=(0.d0,1.d0))
	double complex zx,zf,zd,zc,zdelta,za,zb
c
c  Evaluate CF
c
	zf=tiny
	zc=zf
	zd=0.d0
	elph2=(el+.5d0)**2
	do j=1,maxit
		zb=2.d0*(zx+j*zi)
		za=(dble(j)-.5d0)**2-elph2
		zd=zb+za*zd
		if(zabs(zd).lt.tiny) zd=tiny
		zc=zb+za/zc
		if(zabs(zc).lt.tiny) zc=tiny
		zd=1.d0/zd
		zdelta=zc*zd
		zf=zf*zdelta
		if (zabs(zdelta-1.d0).lt.eps) then
			if(iprint.gt.0) print *,'sphan2cf: Iterationen = ', j,real(zx)/el
			goto 11
		endif
	enddo	
c
c  no convergence
c
	print *, '<sphancf>: zx,el,eps,zf,maxit,tiny=',zx,el,eps,zf,maxit,tiny
	stop '*** Iteration failed to converge in sphancf ! ***'
c
c  converged
c
 11	continue
	zf=el+1.d0+zi*zx+zi*zf
      end
