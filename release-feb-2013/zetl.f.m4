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
c ----------------------------------------------------------------------------
c  Computes ratio of (complex) spherical Bessel-functions 
c    z(n):= x*j(n+1)/j(n)
c
c  by calculation of the recurrence formula
c    z(n)=x**2/(2*n+3 -z(n+1)).
c  (Takeuchi and Saito (1972), Methods in computional Physics, p.243)
c
c  The recurrence formula is treated as a continued fraction
c  and evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
c  al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
c  p. 169 for further information.
c
c  Author: J.R. Dalkolmo
c--------------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_function_type,`double complex')
	define(m4_abs,ifelse(m4_function_type,`double complex',`zabs($*)',`dabs($*)'))
	define(m4_sqrt,ifelse(m4_function_type,`double complex',`zsqrt($*)',`dsqrt($*)'))
c
      subroutine zetl(zx2,el,z,iprint)
	include 'pi.h'
      double precision tiny,el,eps
      integer maxit,i,iprint
      parameter(maxit=10000,tiny=1.d-30,eps=1.d-6)
      m4_function_type zx2,z,zd,zc,zdelta 
c
c First iteration
c
      z=tiny
      zd=2.d0*el+3.d0
      zc=zd+zx2/z
      zd=1.d0/zd
      zdelta=zc*zd
      z=z*zdelta
c
c Remaining iterations
c
      do i=2,maxit
        zd=2.d0*(el+dble(i))+1.d0-zx2*zd
        if (zd.eq.0.d0) zd=tiny
        zc=2.d0*(el+dble(i))+1.d0-zx2/zc
        if (zc.eq.0.d0) zc=tiny
        zd=1.d0/zd
        zdelta=zc*zd
        z=z*zdelta
        if (dabs(m4_abs(zdelta)-1.d0).lt.eps) then
		if(iprint.gt.2) print *,'zetl: Iterationen = ', i,real(m4_sqrt(zx2))/el
		return
        endif
      enddo
      print *, 'zx2,el,eps,z,maxit,tiny=',zx2,el,eps,z,maxit,tiny
      stop '*** Iteration failed to converge in ZETL ! ***'
      end
