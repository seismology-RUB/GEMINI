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
c-------------------------------------------------------------------
c  For given value of xx and an ordered sequence x(i), i=1...n
c  find the index j such that xx lies between x(j) and x(j+1).
c  If xx sits exactly on one of the x(i), j is chosen such that
c  xx sits on the greater one of x(j) and x(j+1). This means:
c  x(i) ascending:  x(j) < xx <= x(j+1)
c  x(i) descending  x(j) >= xx > x(j+1)
c  double precision version
c-------------------------------------------------------------------
	subroutine dlocate(xx,n,x,j)
	integer n,j
	double precision xx,x(n)
	integer lidx,ridx,midx
c
c  start with lower and upper index bounds 0 and n+1
c
	lidx = 0
	ridx = n+1
	do while (ridx-lidx .gt. 1)
		midx = (lidx+ridx)/2
c
c  x(i) is an ascending sequence
c
		if (x(n) .gt. x(1)) then
			if (xx .gt. x(midx)) then
				lidx = midx
			else
				ridx = midx
			endif
c
c  x(i) is a descending sequence
c
		else
			if (xx .gt. x(midx)) then
				ridx = midx
			else
				lidx = midx
			endif
		endif
	enddo
	j = lidx
	return
	end
