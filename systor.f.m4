c----------------------------------------------------------------------
c  $Id: systor.f,v 1.1 2003/05/30 23:26:31 wolle Exp $
c
c  Evalate system matrix for ODE in liquid medium
c----------------------------------------------------------------------
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	define(m4_function_type, `double complex')
	define(m4_elcon_type, `double complex')
c----------------------------------------------------------------------
	subroutine systor(r,a,b)
	include 'nvmax.h'
	m4_function_type a(nvmax,nvmax),b(nvmax),zom,zomro
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zmu,zkap
	double precision ro,om,dll1,r,rr,rr2
	common/degree/dll1
	common/omega/zom,om
c
	call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
	zomro=zom*zom*ro
	rr=1.d0/r
	rr2=rr*rr
	a(1,1)=rr
	a(1,2)=1.d0/zpl
	a(2,1)=-zomro-zpn*rr2*(2.d0-dll1)
	a(2,2)=-3.*rr
c
	b(1)=0.d0
	b(2)=0.d0
c
	return
	end
