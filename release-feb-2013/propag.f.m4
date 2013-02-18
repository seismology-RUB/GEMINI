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
c  Propagate ODE through layers from rbeg to rend
c  and store solution at nodes
c
c  rbeg:       starting point
c  rend:       end point
c  ystart:     initial value of solution vector
c  sysmat:     routine to evaluate system matrix
c  nvar:       number of components of solution vector
c  csy:        type of ODE
c  nstore:     number of stored values of solution vector (output)
c  xstore:     x-values where solution is stored (output)
c  ystore:     stored values of solution vector (output)
c---------------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_function_type,`double complex')
c---------------------------------------------------------------------
	subroutine propag(rbeg,rend,ystart,nvar,csy,sysmat)
	include 'zero.h'
	integer nvar,nl,top,idir,nla,nle,ifl
	double precision rbeg,rend,x1,x2,om,h1
	external sysmat,sysliq
	character*1 csy
	m4_function_type ystart(nvar),zom

ifdef(`MODES',`
	integer rt_docount,rt_eif,rt_count
 	common/roots/rt_docount,rt_eif,rt_count
')
	common/omega/zom,om
c
c  do nothing if rend=rbeg
c
	if(dabs(rend-rbeg).lt.zero) return
	idir=int(sign(1.d0,rend-rbeg))
c
c  find layer index of rbeg
c  if rbeg coincides with a discontinuity, assume layer above if idir=+1
c  and layer below if idir=-1
c
	call layer_getindex(rbeg,nla,top)
	if(top.eq.1.and.idir.gt.0) nla=nla+1
c
c  find layer of rend. If rend coincides with a discontinuity
c  stop below if idir>0 and above if idir<0
c
	call layer_getindex(rend,nle,top)
	if(top.eq.1.and.idir.lt.0) nle=nle+1
c
c  Step through layers 
c  boundary conditions are not checked in the first layer.
c  After each layer, the solution is evaluated at the nodes
c
	x1=rbeg
c	print *,'propag: nla,nle,idir: ',nla,nle,idir
	do nl=nla,nle,idir
		call layer_stepsize(nl,x1,om,h1)
		call layer_setnlactive(nl)
		call layer_isfluid(nl,ifl)
		if(idir.gt.0) then
			call layer_getrb(nl,x2)
			x2=min(rend,x2)
		else
			call layer_getrb(nl-1,x2)
			x2=max(rend,x2)
		endif
c		print *,'propag: nl,x1,x2,h1 = ',nl,x1,x2,h1
		if(nl.ne.nla) call boundary(nl,ystart,nvar,csy,idir)
		if(ifl.eq.1) then
			call bsint_driver(ystart,2,x1,x2,h1,sysliq)
ifdef(`GREEN',`
			call nodes_getsolution(nl,2,x1,x2,idir)
')
ifdef(`MODES',`
			if(rt_docount.eq.1) call rootcount(2)
			if(rt_eif.gt.0) call nodes_getsolution(nl,2,x1,x2,idir)
')
		else
			call bsint_driver(ystart,nvar,x1,x2,h1,sysmat)
ifdef(`GREEN',`
			call nodes_getsolution(nl,nvar,x1,x2,idir)
')
ifdef(`MODES',`
			if(rt_docount.eq.1) call rootcount(nvar)
			if(rt_eif.gt.0) call nodes_getsolution(nl,nvar,x1,x2,idir)
')
		endif
		x1=x2
	enddo
	return
	end
c------------------------------------------------------------------------     
c   check for liquid <-> solid transitions and apply boundary conditions
c   transition is from layer nl to nl-idir
c
	subroutine boundary(nl,zy,nvar,csy,idir)
	integer nl,idir,ifl1,ifl2,nvar
	m4_function_type zy(nvar)
	character csy*1
c
	call layer_isfluid(nl,ifl1)
	call layer_isfluid(nl-idir,ifl2)
	if(ifl1+ifl2.eq.1) then
c
c Boundary solid -> liquid
c
		if(ifl1.gt.ifl2) then
			if(csy.eq.'M') then
				zy(1)=zy(3)
				zy(2)=zy(5)
			endif
c
c Boundary liquid -> solid
c
		else
			if(csy.eq.'M') then
				zy(4)=zy(2)
				zy(2)=zy(1)
				zy(1)=0.d0
				zy(3)=0.d0
				zy(5)=0.d0
			else if(csy.eq.'C' .or. csy.eq.'P') then
				zy(3)=1.d0
				zy(4)=0.d0
			endif
		endif
	endif
	return
	end
