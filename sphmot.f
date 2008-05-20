c---------------------------------------------------------------
c  Compute minor vectors from starting radius to surface or
c  top of solid part of medium (if water layer is on top)
c  or downward from surface to starting radius
c  or downward from surface to bottom of water layer (if present)
c---------------------------------------------------------------
	subroutine sphmot(el,nlstrt,rstart,ystart,knt,vertadd,iprint,ieif,docount)
	include 'nvmax.h'
	double precision rstart,el,om,elp1,re,rsolid,vertadd
	double precision ystart(nvmax),yh(nvmax)
	integer istone,knt,kount,iprint
	integer docount,nlstrt,ieif,nlay,ifliq,nvar
	integer rt_docount,rt_eif,rt_count
	integer j
	common/roots/rt_docount,rt_eif,rt_count
	common/omega/om
	common/degree/elp1
	external sysminors
c
c Get bottom radius of any integration
c also determine if Stoneley mode is hanging around
c
	nvar=5
	rt_docount=docount
	rt_eif=ieif
	elp1=el*(el+1.d0)
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	call layer_getrb(nlay,re)
	call startr(vertadd,rstart,nlstrt,istone,iprint)
c
	if(iprint.gt.2) print *,'Sphmot: ieif,nlstrt,el,x1: ',ieif,nlstrt,el,rstart
c
c  jump to 15 for downward integration of minors
c  for eigenfunction construction
c
	if(ieif.eq.2) goto 15
c
c  jump to 5 for upward integration of eigenfunction
c
	if(ieif.eq.1) goto 5
c---------------------------------------------------
c  this part for root bracketing
c  in case of Scholte wave, use root count
c  from integration in solid part of medium
c  also get decaying solution in water by integrating
c  from the surface to the water bottom
c--------------------------------------------------- 
	call layer_getrb(nlay-1,rsolid)
	rt_count=0
	if(ifliq.eq.1.and.istone.eq.1) then
		if(iprint.gt.1) then
			print *,'<sphmot>: Scholte wave detected!',kount
			print *,'<sphmot>: no turning point in model! l= ',el
		endif
		call stavani(ystart,nvar,rstart,nlstrt,'M',iprint)
		call propag(rstart,rsolid,ystart,nvar,'M',sysminors)
		rt_docount=0
		yh(1)=1.d0
		yh(2)=0.d0
		call propag(re,rsolid,yh,nvar,'M',sysminors)
		ystart(2)=-ystart(3)*yh(2)+ystart(5)*yh(1)
		ystart(1)=0.d0
		rt_docount=docount
	else
		call stavani(ystart,nvar,rstart,nlstrt,'M',iprint)
		call propag(rstart,re,ystart,nvar,'M',sysminors)
		if(nlstrt.eq.nlay.and.ifliq.eq.1) then
			ystart(2)=-ystart(2)
		endif
	endif
	knt=rt_count
	return
c----------------------------------------------
c  eigenfunction calculation upwards -> m^(1)
c  beginning with start radius
c  m1 in liquid layer is not used anymore
c----------------------------------------------
 5	call stavani(ystart,nvar,rstart,nlstrt,'M',iprint)
	call propag(rstart,re,ystart,nvar,'M',sysminors)
	return
c-----------------------------------------------
c  eigenfunction calculation downwards -> m^(2)
c  starting at the surface
c-----------------------------------------------
 15	continue
	if(ifliq.eq.0) then
 		ystart(1)=0.d0
		ystart(2)=1.d0
		ystart(3)=0.d0
		ystart(4)=0.d0
		ystart(5)=0.d0
	else
		ystart(1)=1.d0
		ystart(2)=0.d0
	endif
	call propag(re,rstart,ystart,nvar,'M',sysminors)
c
	return
	end
