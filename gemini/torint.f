c----------------------------------------------------------------------------------------
c  $Id: torint.f,v 1.2 2003/06/01 07:44:50 wolle Exp $
c
c                         TORINT
c
c  Compute two solutions by integrating
c  from a starting radius to max(re,rs) and by integrating
c  from the surface to min(rs,re).
c  The solutions are later used to construct Green functions
c  at the receiver nodes.
c
c  rs:			radius of source
c  re:			radius of receiver
c  gh:			solution receiver nodes gh(2,below-above-source,node) -> output
c  det:			det(rs-or-re), determinant at rs and re formed from gh -> output
c  rst:			starting radius
c  nnsol:			uppermost solid node
c  iprint:			print level
c-----------------------------------------------------------------------------------------
	subroutine torint(nls,rs,nle,re,nlstrt,rst,gh,det,nnsol,iprint)
	include 'nvmax.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	double precision rs,re,rbeg,rend,rst,elp1
	double complex gh(2,2,nnd),det(2),ystart(nvmax)
	integer nvar,nlstrt,nr,i,j,iprint,nnsol
	integer nls,nle,nlay,ifla,ifls,ifle
	character csyst*1
	common/degree/elp1
	external systor
c
c  Get starting values for integration UPWARDS to source/receiver
c  minors and mtil get second and third index of 1,respectively
c
	nvar=2
	csyst = 'T'
	call layer_getnlay(nlay)
	call layer_isfluid(nls,ifls)
	call layer_isfluid(nle,ifle)
	call layer_isfluid(nlay,ifla)
c
c  get starting values
c
	call stavani(ystart,nvar,rst,nlstrt,csyst,iprint)
c
c  zero gh (first index: component, second index: below (1) or above (2) rs,re)
c
	do nr=1,nnod
		do j=1,2
			do i=1,2
				gh(i,j,nr)=0.d0
			enddo
		enddo
	enddo
c
c  solution BELOW the source
c  integrate upwards and stop at nodes below rend
c
	rend=max(rs,re)
	call propag(rst,rend,ystart,nvar,csyst,systor)
c
c  fill gh
c  use only nodes above starting radius
c
	call dlocate(rst,nnod,rnod,jstu)
	if(jstu.eq.0) jstu=1
 1	if(rnod(jstu).le.rst) then
		jstu=jstu+1
		goto 1
	endif
	do nr=jstu,max(jsd,jed)
		do i=1,nvar
			gh(i,1,nr)  = ynod(i,nr)
		enddo
		if(iprint.gt.0) write(6,'(i5,5e13.3)')  nr,(zabs(ynod(i,nr)),i=1,nvar)
	enddo
c
c  Solution ABOVE the source
c  Integration downwards to the source/receiver
c  stavsurf returns rbeg, changed if water layer on top
c
	if(ifla.eq.1) then
		call stavsurf(ystart,nvar,rbeg,csyst,iprint)
		call layer_getrb(nlay-1,rbeg)
		nnsol=jwd
	else
		call stavsurf(ystart,nvar,rbeg,csyst,iprint)
		nnsol=nnod
	endif
c
c  stop at receivers above rend
c
	rend=min(rs,re)
	call propag(rbeg,rend,ystart,nvar,csyst,systor)
	do nr=nnsol,min(jeu,jsu),-1
		do i=1,2
			gh(i,2,nr)  =  ynod(i,nr)
		enddo
		if(iprint.gt.0) write(6,'(i5,5e13.3)') nr,(zabs(ynod(i,nr)),i=1,nvar)
	enddo
c
c  Determinant at rs and re: det(rs-or-re)
c
	det(1)=gh(1,1,jsd)*gh(2,2,jsu)-gh(1,2,jsu)*gh(2,1,jsd)
	det(2)=gh(1,1,jed)*gh(2,2,jeu)-gh(1,2,jeu)*gh(2,1,jed)
	if(iprint.gt.0) print *,'<torint>: Dets: ',det(1),det(2)
c
	return
	end
