c------------------------------------------------------------------
c  $Id: syssolidt.f,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c   Negative transpose of system matrix A
c------------------------------------------------------------------
	subroutine syssolidt(r,a,b)
	include 'nvmax.h'
	double complex a(nvmax,nvmax),b(nvmax)
	double precision dll1,ro,r,rr,rr2,om
	double complex zom,zpa,zpc,zpf,zpl,zpn,zkap,zmu,zrpc,z1,z2,zdlz1,zomro
	common/degree/dll1
	common/omega/zom,om
c
	rr=1.d0/r
	rr2=rr*rr
	call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
	zrpc=1.d0/zpc
	z1=zpa-zpf*zpf*zrpc-zpn
	z2=zpf*rr*zrpc
	zdlz1 = dll1*z1
	zomro = zom*zom*ro
c
	a(1,1)=2.d0*z2
	a(2,1)=-zrpc
	a(3,1)=-dll1*z2
	a(4,1)=0.d0
	a(1,2)=zomro-4.d0*z1*rr2
	a(2,2)=-a(1,1)+2.d0*rr
	a(3,2)=2.d0*rr2*zdlz1
	a(4,2)=-dll1*rr
	a(1,3)=rr
	a(2,3)=0.d0
	a(3,3)=-a(1,3)
	a(4,3)=-1.d0/zpl
	a(1,4)=2.d0*z1*rr2
	a(2,4)=z2
	a(3,4)=zomro - ( zdlz1 + (dll1-2.d0)*zpn )*rr2
	a(4,4)=3.d0*rr
c
	b(1)=0.d0
	b(2)=0.d0
	b(3)=0.d0
	b(4)=0.d0
c
	return
	end
