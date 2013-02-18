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
c--------------------------------------------------------
c  calculate sqrt(x**2+y**2) in a robust way
c
	real function pythag(x,y)
	real x,y,ax,ay,r
c
c  absolute values
c
	ax = abs(x)
	ay = abs(y)
c
c  case where abs(x) >= abs(y)
c
	if (ax .ge. ay) then
		if (ax .eq. 0.) then
			pythag = 0.0
		else
			r = ay/ax
			pythag = ax*sqrt(r*r+1.)
		endif
c
c  case where abs(x) < abs(y)
c
	else
		r = ax/ay
		pythag = ay*sqrt(r*r+1.)
	endif
	return
	end