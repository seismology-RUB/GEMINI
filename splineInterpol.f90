!--------------------------------------------------------------
!  spline interpolation
!--------------------------------------------------------------
 module splineInterpol
	interface splineInterpol2P
		module procedure splineInterpol2P_real
		module procedure splineInterpol2P_double_precision
	end interface
 contains
!--------------------------------------------------------------------
	subroutine splineInterpol2P_real(x,xl,xr,yl,yr,ypl,ypr,y,yp)
	real, parameter :: zero = 1.e-6
	real :: x,xl,xr,yl,yr,ypl,ypr,y,yp
	real :: h,h2,p,q,a2,a3,t
!
	h=xr-xl
	if(abs(h).lt.zero) then
		y=yl
		yp=ypl
		return
	endif
	h2=h*h
!
!  quadratic and cubic coefficients for
!  spline interpolation of eigenfunctions
!
	p=yr-yl-h*ypl
	q=ypr-ypl
	a2=(3.d0*p-h*q)/h2
	a3=(h*q-2.d0*p)/(h2*h)
!
	t=x-xl
	y=yl+t*(ypl+t*(a2+t*a3))
	yp=ypl+t*(2.*a2+t*3.*a3)
!
	end subroutine splineInterpol2P_real
!---------------------------------------------------------------------
	subroutine splineInterpol2P_double_precision(x,xl,xr,yl,yr,ypl,ypr,y,yp)
	double precision, parameter :: zero = 1.d-10 
	double precision :: x,xl,xr
	double precision :: h,h2,t
	real :: yl,yr,ypl,ypr,y,yp
	real p,q,a2,a3
!
	h=xr-xl
	if(dabs(h).lt.zero) then
		y=yl
		yp=ypl
		return
	endif
	h2=h*h
!
!  quadratic and cubic coefficients for
!  spline interpolation of eigenfunctions
!
	p=yr-yl-h*ypl
	q=ypr-ypl
	a2=(3.d0*p-h*q)/h2
	a3=(h*q-2.d0*p)/(h2*h)
!
	t=x-xl
	y=yl+t*(ypl+t*(a2+t*a3))
	yp=ypl+t*(2.*a2+t*3.*a3)
!
	end subroutine splineInterpol2P_double_precision
!
 end module splineInterpol
