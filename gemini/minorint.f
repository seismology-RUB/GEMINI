c----------------------------------------------------------------------------------------
c  $Id: minorint.f,v 1.2 2003/05/30 23:26:31 wolle Exp $
c
c                         MINORINT
c
c  Compute two minor-vectors by integrating
c  from a starting radius to max(re,rs) and by integrating
c  from the surface to min(rs,re) but at most the starting radius.
c  The minor vectors are later used to construct Green functions
c  at the receiver nodes.
c
c  rs:			radius of source
c  re:			radius of receiver
c  minors:			minors(6,below-above-source,rs-or-re)  -> output
c  mtil:			minor matrix at receiver nodes mtil(1-4,1-4,below-above-source,node) -> output
c  det:			det(rs-or-re), determinant at rs and re formed from minors -> output
c  rst:			starting radius
c  iprint:			print level
c-----------------------------------------------------------------------------------------
	subroutine minorint(nls,rs,nle,re,nlstrt,rst,minors,mtil,det,iprint)
	include 'nvmax.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	double precision rs,re,rbeg,rend,rst,elp1
	double complex minors(6,2,2),mtil(4,4,2,nnd),det(2),ystart(nvmax)
	integer nvar,nlstrt,nr,i,j,iprint
	integer index(4,4),nls,nle,nlay,ifla,ifls,ifle
	character csyst*1
	common/degree/elp1
	external sysminors
c
c
c Evaluate index-array for mtil/dmtil-matrices. Only 6 values
c according to the number of minors are necessary.

	index(1,2) = 1
	index(1,3) = 2
	index(1,4) = 3
	index(2,3) = 4
	index(2,4) = 5
	index(3,4) = 6
c
c  Get starting values for integration UPWARDS to source/receiver
c  minors and mtil get second and third index of 1,respectively

	nvar=5
	csyst = 'M'
	call layer_getnlay(nlay)
	call layer_isfluid(nls,ifls)
	call layer_isfluid(nle,ifle)
	call layer_isfluid(nlay,ifla)
c
c  get starting values

	call stavani(ystart,nvar,rst,nlstrt,csyst,iprint)
c
c  zero minors and mtil

	do nr=1,nnod
		do j=1,4
			do i=1,4
				mtil(i,j,1,nr)=0.d0
				mtil(i,j,2,nr)=0.d0
			enddo
		enddo
	enddo
	do nr=1,2
		do i=1,6
			minors(i,1,nr)=0.d0
			minors(i,2,nr)=0.d0
		enddo
	enddo
c
c  Minors BELOW the source
c  integrate upwards and stop at nodes below rend
c  but above starting radius
c  Evaluate (antisymmetric) Matrix M~
c  and store minors
c  minors(component,up-or-down,rs-or-re)

	rend=max(rs,re)
	call propag(rst,rend,ystart,nvar,csyst,sysminors)
	do i=1,nvar
		minors(i,1,1)=ynod(i,jsd)
		minors(i,1,2)=ynod(i,jed)
	enddo
	minors(6,1,1)=-ynod(1,jsd)/elp1
	minors(6,1,2)=-ynod(1,jed)/elp1
c
	if(iprint.gt.0) then
		print *,'Minors below source'
		write(6,'(a5,5a13)') 'nr','y1','y2','y3','y4','y5'
	endif
c
c  fill mtil
c  use only nodes above starting radius
c
	call dlocate(rst,nnod,rnod,jstu)
	if(jstu.eq.0) jstu=1
 1	if(rnod(jstu).le.rst) then
		jstu=jstu+1
		goto 1
	endif
	do nr=jstu,max(jsd,jed)
		ynod(6,nr)=-ynod(1,nr)/elp1
		do i=1,4
			do j=i,4
				if (i.eq.j) then
					mtil(i,j,1,nr)  = 0.d0
				else
					mtil(i,j,1,nr)  = ynod(index(i,j),nr)
					mtil(j,i,1,nr)  = -mtil(i,j,1,nr)
				endif
			enddo
		enddo
		if(iprint.gt.0) write(6,'(i5,5e13.3)')  nr,(zabs(ynod(i,nr)),i=1,nvar)
	enddo
c
c  Minors ABOVE the source
c  Integration downwards to the source/receiver or to the starting radius
c  minors and mtil get second and third index of 2,respectively
c  stavsurf returns rbeg!

	rend=min(rs,re)
	if(ifla.eq.1) then
		ystart(1)=1.d0
		ystart(2)=0.d0
		call layer_getrb(nlay,rbeg)
	else
		call stavsurf(ystart,nvar,rbeg,csyst,iprint)
	endif
c
c  stop at receivers above rend
c  Evaluate (antisymmetric) Matrix M~

	call propag(rbeg,rend,ystart,nvar,csyst,sysminors)
	do i=1,nvar
		minors(i,2,1)=ynod(i,jsu)
		minors(i,2,2)=ynod(i,jeu)
	enddo
	minors(6,2,1)=-ynod(1,jsu)/elp1
	minors(6,2,2)=-ynod(1,jeu)/elp1
	if(iprint.gt.0) then
		print *,'Minors above source'
		write(6,'(a5,5a13)') 'nr','y1','y2','y3','y4','y5'
	endif
	do nr=nnod,min(jeu,jsu),-1
		ynod(6,nr)=-ynod(1,nr)/elp1
		do i=1,4
			do j=i,4
				if (i.eq.j) then
					mtil(i,j,2,nr)  =  0.d0
				else
					mtil(i,j,2,nr)  =  ynod(index(i,j),nr)
					mtil(j,i,2,nr)  = -mtil(i,j,2,nr)
				endif
			enddo
		enddo
		if(iprint.gt.0) write(6,'(i5,5e13.3)') nr,(zabs(ynod(i,nr)),i=1,nvar)
	enddo
c
c  Determinant at rs and re: det(rs-or-re)
c  minors(component,up-or-down,rs-or-re)
c
	do i=1,2
		det(i) =  minors(1,1,i)*minors(6,2,i) + minors(6,1,i)*minors(1,2,i) -
     &			minors(2,1,i)*minors(5,2,i) + minors(3,1,i)*minors(4,2,i) +
     &			minors(4,1,i)*minors(3,2,i) - minors(5,1,i)*minors(2,2,i)
	enddo
c
c  exceptions when either rs or re are in liquid layer
c
	if(ifls.eq.1) then
		det(1)=minors(1,1,1)*minors(2,2,1)-minors(2,1,1)*minors(1,2,1)
	endif
	if(ifle.eq.1) then
		det(2)=minors(1,1,2)*minors(2,2,2)-minors(2,1,2)*minors(1,2,2)
	endif
	if(iprint.gt.0) print *,'Dets: ',det(1),det(2)
	if(iprint.gt.0) print *,'Minors(1,1): ',minors(1,1,1),minors(1,1,2)
c
	return
	end
