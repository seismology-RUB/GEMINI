c--------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
c------------------------------------------------------------------
c   Spheroidal system matrix A, solid medium
c------------------------------------------------------------------
	subroutine sysspher(r,a,b)
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
	a(1,1)=-2.d0*z2
	a(1,2)=zrpc
	a(1,3)=dll1*z2
	a(1,4)=0.d0
	a(2,1)=-zomro+4.d0*z1*rr2
	a(2,2)=-a(1,1)-2.d0*rr
	a(2,3)=-2.d0*rr2*zdlz1
	a(2,4)=dll1*rr
	a(3,1)=-rr
	a(3,2)=0.d0
	a(3,3)=-a(3,1)
	a(3,4)=1.d0/zpl
	a(4,1)=-2.d0*z1*rr2
	a(4,2)=-z2
	a(4,3)=-zomro + ( zdlz1 + (dll1-2.d0)*zpn )*rr2
	a(4,4)=-3.d0*rr
c
	b(1)=0.d0
	b(2)=0.d0
	b(3)=0.d0
	b(4)=0.d0
c
	return
	end
