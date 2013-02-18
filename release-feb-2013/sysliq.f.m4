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
