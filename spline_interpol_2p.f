c------------------------------------------------------------------------
c  cubic interpolation of function between two points
c  using derivative information 
c  x:		location where function is evaluated
c  xl:		left point where function is known
c  xr:		right point where function is known
c  yl:		left function value
c  yr:		right function value		
c  ypl,ypr:	same for derivative of function
c  y,yp:	interpolated value and derivative at x
c-----------------------------------------------------------
	subroutine spline_interpol_2p(x,xl,xr,yl,yr,ypl,ypr,y,yp)
	real zero
	parameter(zero=1.e-6)
	real x,xl,xr,yl,yr,ypl,ypr,y,yp
	real h,h2,p,q,a2,a3,t
c
	h=xr-xl
	if(abs(h).lt.zero) then
		y=yl
		yp=ypl
		return
	endif
	h2=h*h
c
c  quadratic and cubic coefficients for
c  spline interpolation of eigenfunctions
c-
	p=yr-yl-h*ypl
	q=ypr-ypl
	a2=(3.d0*p-h*q)/h2
	a3=(h*q-2.d0*p)/(h2*h)
c
	t=x-xl
	y=yl+t*(ypl+t*(a2+t*a3))
	yp=ypl+t*(2.*a2+t*3.*a3)
c
	return
	end
