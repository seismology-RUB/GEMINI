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
c---------------------------------------------------------------------
c  Calculate coefficients for a (natural) cubic spline interpolation
c  Follows algorithm described in Numerical Recipes 3.3.
c  Uses Lapack routine for solution of tridiagonal system
c  x:     array with abscissa values at which 
c  y:     function values are given
c  n:     number of function values
c  y2:    values of second derivative (output) with y2(1) = y2(n) = 0
c  work:  (n x 3)- matrix for work space
c---------------------------------------------------------------------
c  double precision version
c---------------------------------------------------------------------
	subroutine dcubsplco(x,y,n,y2,work)
	implicit none
	integer n,info,j
	double precision x(n),y(n),y2(n),work(n,3)
c
c  set up tridiagonal matrix
c
	do j = 2,n-1
		work(j-1,1) = (x(j)-x(j-1))/6.
		work(j,2) = (x(j+1)-x(j-1))/3.
		work(j,3) = (x(j+1)-x(j))/6.
		y2(j) = (y(j+1)-y(j))/(x(j+1)-x(j)) - (y(j)-y(j-1))/(x(j)-x(j-1))
	enddo
	work(1,2) = x(2)-x(1)
	work(n,2) = x(n)-x(n-1)
	work(1,3) = 0.0
	work(n-1,1) = 0.0
	y2(1) = 0.0
	y2(n) = 0.0
c
c  solve tridiagonal system for y2
c
	call dgtsv(n,1,work(1,1),work(1,2),work(1,3),y2,n,info)
	if (info .ne. 0) then
		print *,'CUBSPLCO: tridiagonal system could not be solved'
		stop
	endif
	return
	end
c-----------------------------------------------------------------------