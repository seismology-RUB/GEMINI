c---------------------------------------------------------------
c  Compute minor vectors from starting radius to surface or
c  top of solid part of medium (if water layer is on top)
c  or downward from surface to starting radius
c  or downward from surface to bottom of water layer (if present)
c---------------------------------------------------------------
	subroutine tormot(el,nlstrt,rstart,ystart,knt,vertadd,iprint,ieif,docount)
	include 'nvmax.h'
	double precision rstart,el,om,elp1,re,vertadd
	double precision ystart(nvmax)
	integer istone,knt,iprint
	integer docount,nlstrt,ieif,nlay,ifliq,nvar
	integer rt_docount,rt_eif,rt_count
	common/roots/rt_docount,rt_eif,rt_count
	common/omega/om
	common/degree/elp1
	external systor
c
c Get bottom radius of any integration
c
	nvar=2
	rt_docount=docount
	rt_eif=ieif
	elp1=el*(el+1.d0)
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	if(ifliq.eq.0) then
		call layer_getrb(nlay,re)
	else
		call layer_getrb(nlay-1,re)
	endif
	call startr(vertadd,rstart,nlstrt,istone,iprint)
c-
	if(iprint.gt.2) print *,'Tormot: ieif,nlstrt,el,rstart: ',ieif,nlstrt,el,rstart
c
c  jump to 5 for upward integration of eigenfunction
c
	if(ieif.eq.1) goto 5
c---------------------------------------------------
c  this part for root bracketing
c  no root if starting radius is in water layer
c--------------------------------------------------- 
	rt_count=0
	if(ifliq.eq.1.and.nlstrt.eq.nlay) return
	call stavani(ystart,nvar,rstart,nlstrt,'T',iprint)
	call propag(rstart,re,ystart,nvar,'T',systor)
	knt=rt_count
	return
c----------------------------------------------
c  eigenfunction calculation upwards -> m^(1)
c  beginning with start radius
c  m1 in liquid layer is not used anymore
c----------------------------------------------
 5	if(ifliq.eq.1.and.nlstrt.eq.nlay) return
	call stavani(ystart,nvar,rstart,nlstrt,'T',iprint)
	call propag(rstart,re,ystart,nvar,'T',systor)
	return
c
	end
