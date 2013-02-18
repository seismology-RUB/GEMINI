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
c----------------------------------------------------------
c  Procedures to implement layer description of earth model
c----------------------------------------------------------
c        subroutine earthmodel_layer
c        subroutine layer_getiktop(nl,j)
c        subroutine layer_getikbot(nl,j)
c        subroutine layer_printinfo
c        subroutine layer_getindex(r,nl,top)
c        subroutine layer_stepsize(nl,r,om,h)
c        subroutine layer_set_stepsize_wlfrac(wlfrac)
c        subroutine layer_setnlactive(nl)
c        subroutine layer_getnlactive(nl)
c        subroutine layer_isfluid(nl,j)
c        subroutine layer_getnlay(n)
c        subroutine layer_getrb(nl,r)
c        subroutine layer_gettopnode(nl,j)
c        subroutine layer_getbotnode(nl,j)
c---------------------------------------------------------------------
c   Create layer information from earthmodel
c   Layers are defined by discontinuities
c   Deepest rb is earth's center
c---------------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_elcon_type,`double complex')
	define(m4_abs,ifelse(m4_function_type,`double complex',`zabs($*)',`dabs($*)'))
c---------------------------------------------------------------
	subroutine earthmodel_layer
	include 'earthmodel.h'
	include 'layer.h'
	include 'zero.h'
	integer i
c
	layer_rb(0)=0.d0
	layer_nlay=0
	layer_iktop(0)=0
	layer_ikbot(1)=0
	do i=1,earthmodel_nk-1
		if(dabs(earthmodel_rk(i+1)-earthmodel_rk(i)).lt.zero) then
			layer_nlay=layer_nlay+1
			if(m4_abs(earthmodel_zl(i)).lt.zero) then
				layer_iflso(layer_nlay)=1
			else
				layer_iflso(layer_nlay)=0
			endif
			layer_iktop(layer_nlay)=i
			layer_ikbot(layer_nlay+1)=i+1
			layer_rb(layer_nlay)=earthmodel_rk(i)
		endif
	enddo
c
c  add top layer
c
	layer_nlay=layer_nlay+1
	layer_rb(layer_nlay)=earthmodel_rk(earthmodel_nk)
	layer_iktop(layer_nlay)=earthmodel_nk
	if(m4_abs(earthmodel_zl(earthmodel_nk)).lt.zero) then
		layer_iflso(layer_nlay)=1
	else
		layer_iflso(layer_nlay)=0
	endif
c
c  if iflso(2)=0 we have a fluid core and an inner core
c
	if(layer_nlay.gt.2.and.layer_iflso(2).eq.1.and.layer_iflso(1).eq.0) then
		layer_noc=2
		layer_nic=1
	else
		layer_nic=0
		layer_noc=0
	endif
	return
	end
c--------------------------------------------------------------
c   Get model node index of top layer boundary
c
	subroutine layer_getiktop(nl,j)
	include 'layer.h'
	integer nl,j
	j=layer_iktop(nl)
	return
	end
c--------------------------------------------------------------
c   Get model node index of bottom layer boundary
c
	subroutine layer_getikbot(nl,j)
	include 'layer.h'
	integer nl,j
	j=layer_ikbot(nl)
	return
	end	
c-------------------------------------------------------------
c   print layer info
c
	subroutine layer_printinfo
	include 'layer.h'
	integer nl
c
	print *,'RB: ',layer_nlay,(layer_rb(layer_nlay)-layer_rb(nl),nl=0,layer_nlay)
	print *,'iflso: ',(layer_iflso(nl),nl=1,layer_nlay)
	print *,'nic,noc: ',layer_nic,layer_noc
	return
	end
c--------------------------------------------------------------
c   determine layer from radius
c   top=1 if r=rb(nl)
c
	subroutine layer_getindex(r,nl,top)
	include 'layer.h'
	include 'zero.h'
	integer nl,top
	double precision r
c
	top=0
	call dlocate(r,layer_nlay,layer_rb(0),nl)
	if(dabs(r-layer_rb(nl)).lt.zero) top=1
c
	return
	end
c------------------------------------------------------------
c   determine initial stepsize for integration from velocity
c
	subroutine layer_stepsize(nl,r,om,h)
	include 'layer.h'
	include 'pi.h'
	integer nl
	double precision r,h,alf,om,vel
c
	call earthmodel_isovel(nl,r,alf,vel)
	if(layer_iflso(nl).eq.1) then
		vel=alf*layer_bali
	endif
	h=layer_stepsize_wlfrac*2.d0*pi*vel/om
c
	return
	end
c------------------------------------------------------------
c  set wavelength fraction for calculation of initial
c  stepsize
c
	subroutine layer_set_stepsize_wlfrac(wlfrac)
	include 'layer.h'
	double precision wlfrac
c
	layer_stepsize_wlfrac=wlfrac
	return	
	end
c------------------------------------------------------------
c   set active layer index
c
	subroutine layer_setnlactive(nl)
	include 'layer.h'
	integer nl
	layer_nlactive=nl
	return
	end
c------------------------------------------------------------
c   get active layer index
c
	subroutine layer_getnlactive(nl)
	include 'layer.h'
	integer nl
	nl=layer_nlactive
	return
	end
c--------------------------------------------------------------
c  is layer fluid ?
c
	subroutine layer_isfluid(nl,j)
	include 'layer.h'
	integer nl,j
c
	j=layer_iflso(nl)
	return
	end
c----------------------------------------------------------------
c  number of layers
c
	subroutine layer_getnlay(n)
	include 'layer.h'
	integer n
	n=layer_nlay
	return
	end
c----------------------------------------------------------------
c   query top layer boundary
c
	subroutine layer_getrb(nl,r)
	include 'layer.h'
	integer nl
	double precision r
	r=layer_rb(nl)
	return
	end
c-------------------------------------------------------------------
c   query node index of top node of layer
c
	subroutine layer_gettopnode(nl,j)
	include 'layer.h'
	integer nl,j
c
	j=layer_nrne(nl)
	return
	end
c-------------------------------------------------------------------
c   query node index of bottom node of layer
c
	subroutine layer_getbotnode(nl,j)
	include 'layer.h'
	integer nl,j
c
	j=layer_nrna(nl)
	return
	end
