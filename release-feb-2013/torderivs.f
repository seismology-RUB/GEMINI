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
c                           TORDERIVS
c
c  compute derivatives w.r.t. r from DSV and overwrite DSV with W,WP	
c---------------------------------------------------------------------------
	subroutine torderivs(gf)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'layer.h'
	include 'zero.h'
	include 'nvmax.h'
	double complex gf(nnd,3),yh
	double complex a(nvmax,nvmax),b(nvmax)
	double precision r
	integer j,top
	integer nr,nl
c
c  loop over nodes
c
	do nr=jstu,nnod
		r=rnod(nr)
		call layer_getindex(r,nl,top)
		if(top.eq.1.and.dabs(rnod(nr-1)-rnod(nr)).lt.zero) nl=nl+1
		call layer_setnlactive(nl)
		if(layer_iflso(nl).eq.1) then
			continue
		else
			call systor(r,a,b)
			yh=0.d0
			do j=1,2
				yh=yh+a(1,j)*gf(nr,j)
			enddo
			gf(nr,3) = yh
		endif
	enddo
c
	return
	end
