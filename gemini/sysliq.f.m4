c----------------------------------------------------------------------
c  Evalate system matrix for ODE in liquid medium
c----------------------------------------------------------------------
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	define(m4_function_type, `double complex')
	define(m4_elcon_type, `double complex')
c----------------------------------------------------------------------
	subroutine sysliq(r,a,b)
	include 'nvmax.h'
	m4_function_type a(nvmax,nvmax),b(nvmax),zom
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zomro,zmu,zkap
	double precision ro,om,dll1,r
	common/degree/dll1
	common/omega/zom,om
c
	call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
	zomro=zom*zom*ro
	a(1,1)=-2.d0/r
	a(1,2)=1.d0/zpc - dll1/(zomro*r*r)
	a(2,1)=-zomro
	a(2,2)=0.d0
c
	b(1)=0.d0
	b(2)=0.d0
c
	return
	end
