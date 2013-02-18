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
c---------------------------------------------------------
c  trim a string from left, return index k of first
c  non--blank character, returns k=len(string)+1 if 
c  string is all blanks
c---------------------------------------------------------
	subroutine trimleft(string,k)
	character string*(*)
	integer k,nc
c
	k=1
	nc=len(string)
 1	if(string(k:k).ne.' ') return
	if(k.lt.nc) then
		k=k+1
		goto 1
	else
		k=k+1
	endif
	return
	end
c----------------------------------------------------------
c  trim a string from right and return index k of last
c  non--blank chracter, returns k=0 if string is all blanks
c----------------------------------------------------------
	subroutine trimright(string,k)
	character string*(*)
	integer k
c
	k=len(string)
 1	if(string(k:k).ne.' ') return
	if(k.gt.1) then
		k=k-1
		goto 1
	else
		k=k-1
	endif
	return
	end
