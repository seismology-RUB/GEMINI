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
c----------------------------------------------------------------------
c  Evaluate system matrix for ODE of minors
c.......................................................................
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	define(m4_function_type, `double complex')
	define(m4_elcon_type, `double complex')
c----------------------------------------------------------------------
	subroutine sysminors(r,a,b)
	include 'nvmax.h'
	m4_function_type a(nvmax,nvmax),b(nvmax),zom,zomro
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zkap,zmu,zrpc,z1,z2,z3
	double precision om,dll1,r,rr,rr2,ro
	common/omega/zom,om
	common/degree/dll1
c
c  need to compute matrix elements
c
	rr = 1.d0/r
	rr2 = rr*rr
	call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c	write(6,'(7f10.3)') r,ro,zpa,zpc,zpf,zpl,zpn
	zomro = zom*zom*ro
	zrpc = 1.d0/zpc
	z1 = zpa - zpf*zpf*zrpc - zpn
	z2 = 2.d0*zpf*zrpc
	z3 = 4.d0*z1*rr2
c
c Evaluate system matrix of SODE
c
	a(1,1)=-2.d0*rr
	a(1,2)=-2.d0*dll1*z1*rr2
	a(1,3)=dll1*rr
	a(1,4)=-dll1*zpf*zrpc*rr
	a(1,5)=0.d0
	a(2,1)=0.d0
	a(2,2)=(1.d0-z2)*rr
	a(2,3)=1.d0/zpl
	a(2,4)=zrpc
	a(2,5)=0.d0
	a(3,1)=-z2*rr
	a(3,2)=-zomro+( dll1*(z1+zpn) - 2.d0*zpn )*rr2
	a(3,3)=-(3.d0+z2)*rr
	a(3,4)=0.d0
	a(3,5)=a(2,4)
	a(4,1)=-a(1,1)
	a(4,2)=-zomro+z3
	a(4,3)=0.d0
	a(4,4)=-a(2,2)
	a(4,5)=a(2,3)
	a(5,1)=z3
	a(5,2)=0.d0
	a(5,3)=a(4,2)
	a(5,4)=a(3,2)
	a(5,5)=(z2-5.d0)*rr
c
	b(1)=0.d0
	b(2)=0.d0
	b(3)=0.d0
	b(4)=0.d0
	b(5)=0.d0
c
	return
	end
