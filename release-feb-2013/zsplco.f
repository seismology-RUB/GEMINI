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
c----------------------------------------------------------------
c  $Id: zsplco.f,v 1.2 2003/04/03 12:09:38 wolle Exp $
c
c  Spline complex function values by splining real and imaginary
c  part separately. Calls dcubsplco.f: see comments there.
c  x:    abscissa values
c  zy:   complex function values
c  n:    number of points
c  zy2:  complex valued second derivatives (output)
c  y:    workspace of dimension n
c  y2:   workspace of dimension n
c  work: workspace of dimension (n x 3)
c-----------------------------------------------------------------
	subroutine zsplco(x,zy,n,zy2,y,y2,work)
	integer n,i
	double precision x(n),y(n),y2(n),work(n,3)
	double complex zy(n),zy2(n),zi
c 
	zi=dcmplx(0.d0,1.d0)
	do i=1,n
		y(i)=real(zy(i))
c		print *, y(i)
	enddo
	call dcubsplco(x,y,n,y2,work)
	do i=1,n
		zy2(i)=dcmplx(y2(i),0.d0)
c		print *,y2(i)
	enddo
	do i=1,n
		y(i)=dimag(zy(i))
c		print *,y(i)
	enddo
	call dcubsplco(x,y,n,y2,work)
	do i=1,n
		zy2(i)=zy2(i)+zi*y2(i)
	enddo
	return
	end

