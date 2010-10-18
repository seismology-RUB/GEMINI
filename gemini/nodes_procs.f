c----------------------------------------------------------
c  $Id: nodes_procs.f,v 1.4 2003/05/04 07:30:52 wolle Exp $
c
c  Procedures for nodes-class
c----------------------------------------------------------
c  Establish radial nodes where solutions are stored.
c  Discontinuities get a twin storing node because
c  system matrix needs not be continuous at a discontinuity.
c  Source and receiver depths get two nodes each,
c  because solution G(n,n,r,r_s) is discontinuous at r_s.
c  Exceptions happen when rs,re=rb(nlay) or rs,re < rstore.
c  
c  nlstore is the layer index of the smallest radius
c  for which the basis solutions are to be stored
c  Actual storage will start from rstore.
c
c  !! The code assumes that nlstore == 1 !!!
c  This is automatically set in flgevas_init.f
c
c  f:          Frequency
c  zstore:     Bottom depth down to which nodes are created
c  zs:         Source depth
c  ze:         Receiver depth
c  nnmin:      min number of nodes per layer
c------------------------------------------------------------
	subroutine nodes_create(f,wlratio,zstore,zs,ze,nnmin,iprint)
	include 'pi.h'
	include 'layer.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	include 'zero.h'
	double precision f,om,vel,dr,dro,drmax,rs,drub,drnd,wlratio
	double precision rstore,re,rmin,rmax,alf1,vel1,alf2,vel2
	double precision zstore,zs,ze
	integer nlstore,nls,nle,top,nnmin,iprint,nadd
	integer nival,nl,j,jd,ju,na,n
c
	om=2.*pi*f
	drmax=30.
c
	re=layer_rb(layer_nlay)-ze
	rs=layer_rb(layer_nlay)-zs
	rstore=layer_rb(layer_nlay)-zstore
c
c  get layer of rstore and rs and re
c
	call layer_getindex(rstore,nlstore,top)
	if(top.eq.1) nlstore=nlstore+1
	call layer_getindex(rs,nls,top)
	call layer_getindex(re,nle,top)
c
	if(iprint.gt.0) print *,'<nodes_create>: nls,nle,nlstore,rstore: ',nls,nle,nlstore,rstore
c
	nnod=0
c
c  rs or re might be less than rstore but we need nodes there too
c  so we pack some nodes below rstore
c
	if(rs.lt.rstore.and.re.lt.rstore) then
		if(dabs(rs-re).lt.zero) then
			jsd=1
			jsu=2
			jed=jsd
			jeu=jsu
			nnod=nnod+2
			rnod(jsd)=rs
			rnod(jsu)=rs
		else if(rs.lt.re) then
			jsd=1
			jsu=2
			jed=3
			jeu=4
			nnod=nnod+4
			rnod(jsd)=rs
			rnod(jsu)=rs
			rnod(jed)=re
			rnod(jeu)=re
		else if(rs.gt.re) then
			jsd=3
			jsu=4
			jed=1
			jeu=2
			rnod(jsd)=rs
			rnod(jsu)=rs
			rnod(jed)=re
			rnod(jeu)=re
			nnod=nnod+4
		endif
	else if(rs.lt.rstore.and.re.ge.rstore) then
		jsd=1
		jsu=2
		nnod=nnod+2
		rnod(jsd)=rs
		rnod(jsu)=rs
	else if(rs.ge.rstore.and.re.lt.rstore) then
		jed=1
		jeu=2
		nnod=nnod+2
		rnod(jed)=re
		rnod(jeu)=re
	endif
c
c  now step upwards from rstore through all layers and assign nodes
c  without regard of rs and re
c
	do nl=nlstore,layer_nlay
c
c  set dr in layer either from number of minimum nodes per layer
c  or from average wavelength in layer
c
c  layer_iktop(nl) is the model node at the top of layer nl
c  layer_iktop(nl-1)+1 is the model node at the bottom of layer nl
c  check if rk(na) = rmin. This happens if rstore is on a node!
c
		if(nl.eq.nlstore) then
			layer_nrna(nl)=1
			rmin=rstore
			call earthmodel_getuppernodeindex(nl,rmin,na)
 12			if(dabs(earthmodel_rk(na)-rmin)/rmin.lt.1.e-6) then
				na=na+1
				goto 12
			endif
		else
			layer_nrna(nl)=nnod+1
			na=layer_iktop(nl-1)+1+1
			rmin=earthmodel_rk(na-1)
		endif
c
c  In a new layer, the bottom model node is always a solution node
c
		nnod=nnod+1					
		rnod(nnod)=rmin
c
c  Loop here over model nodes within the layer
c  The last solution node will be the model node at the top of the layer
c
		do n=na,layer_iktop(nl)
			rmax=earthmodel_rk(n)
			drnd=rmax-rmin
			drub=drnd/nnmin
			call earthmodel_isovel(nl,rmin,alf1,vel1)
			call earthmodel_isovel(nl,rmax,alf2,vel2)
			vel=0.5*(vel1+vel2)
			if(layer_iflso(nl).eq.1) then
				vel=0.5*(alf1+alf2)*layer_bali
			endif
			dr=2.*pi*vel/om/wlratio
			dro=min(drub,dr)
			nival=nint((rmax-rmin)/dro)
			dr=(rmax-rmin)/nival		
			if(nnod+nival.gt.nnd) then
				print *,'<nodes_create>: more nodes than dimensioned, max = ',nnd
				print *,'<nodes_create>: layer: ',nl,' with rtop = ',rmax
				print *,'<nodes_create>: dr = ',dr
				stop
			endif
			do j=nnod+1,nnod+nival
				rnod(j)=rmin+(j-nnod)*dr
			enddo
			nnod=nnod+nival
			rmin=rmax
		enddo			
		layer_nrne(nl)=nnod	
	enddo
	if(iprint.gt.2) then
		print *,'Rnodes after first step:'
		print *,nnod,((layer_rb(layer_nlay)-rnod(j))*1.d3,j=1,nnod)
	endif
c
c  jump over fitin if re and rs are below rstore
c
	if(rs.lt.rstore.and.re.lt.rstore) goto 11
c
c  fit in source and receiver if rs != re
c  first fit in the lower one to avoid a change
c  of indices when the second one is fit in
c  check if the lower one is in layer 1, then it
c  is already checked in
c
	if(dabs(rs-re).gt.zero) then
		if(re.gt.rs) then
			if(rs.gt.rstore) then
				call nodes_fitin(rs,jd,ju,nadd)
				jsd=jd
				jsu=ju
				call nodes_adjustnrae(nls,nadd)
			endif
			call nodes_fitin(re,jd,ju,nadd)
			jed=jd
			jeu=ju
			call nodes_adjustnrae(nle,nadd)
		else
			if(re.gt.rstore) then
				call nodes_fitin(re,jd,ju,nadd)
				jed=jd
				jeu=ju
				call nodes_adjustnrae(nle,nadd)
			endif
			call nodes_fitin(rs,jd,ju,nadd)
			jsd=jd
			jsu=ju
			call nodes_adjustnrae(nls,nadd)
		endif
c
c  rs=re
c
	else
		call nodes_fitin(rs,jd,ju,nadd)
		jsd=jd
		jed=jd
		jsu=ju
		jeu=ju
		call nodes_adjustnrae(nls,nadd)
	endif
c
 11	if(iprint.gt.2) then 
		print *,'After nodes_fitin of rs and re: '
		print *,nnod,((layer_rb(layer_nlay)-rnod(j))*1.d3,j=1,nnod)
	endif
c
c  indices of bottom of water layer rw- and rw+ 
c  set always rw=rb(nlay-1)
c  dummies if top layer is solid
c
	if(layer_iflso(layer_nlay).eq.1) then
		call dlocate(layer_rb(layer_nlay-1),nnod,rnod,j)
		jwd=j+1
		jwu=j+2
		rnod(jwd)=layer_rb(layer_nlay-1)
		rnod(jwu)=layer_rb(layer_nlay-1)
	else
		jwd=0
		jwu=0
	endif
c
	if(iprint.gt.2) then
		print *,'Final Rnodes: ',(layer_rb(layer_nlay)-rnod(j),j=1,nnod)
		print *,'jsd,jsu,jed,jeu,jwd,jwu: ',jsd,jsu,jed,jeu,jwd,jwu
	endif
	if(iprint.gt.0) print *,'nodes_create: ',nnod,' Nodes created!' 
c
	return
	end
c-------------------------------------------------------------------
c   fit in source or receiver into node order
c
	subroutine nodes_fitin(rs,jd,ju,nadd)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'zero.h'
	double precision rs
	integer jd,ju,nadd,i,j
c
	nadd=0
c
	call dlocate(rs,nnod,rnod,j)
c
c  source coincides with upper node, add one node.
c  index of upper node and following increases by one
c  ju=j+2 and jd=j+1
c
	if(dabs(rs-rnod(j+1)).lt.zero) then
		ju=j+2
		jd=j+1
		nadd=1
c
c  source coincides with lower node, add one node.
c  Lower node keeps index but node above shift by one
c
	else if(dabs(rs-rnod(j)).lt.zero) then
		jd=j
		ju=j+1
		nadd=1
c
c  source lies between the nodes, add two
c
	else
		jd=j+1
		ju=j+2
		nadd=2
	endif
c
c  shift nodes
c
	if(nnod+nadd.gt.nnd) then
		print *,'<nodes_create>: more nodes than dimensioned'
	endif
	do i=nnod+nadd,j+nadd+1,-1
		rnod(i)=rnod(i-nadd)
	enddo
	nnod=nnod+nadd
	rnod(jd)=rs
	rnod(ju)=rs
c
	return
	end	
c------------------------------------------------------------------
c  adjust indices of first and last node in layer
c
	subroutine nodes_adjustnrae(nl1,nadd)
	include 'layer.h'
	integer nl1,nadd,nl
c
	do nl=nl1,layer_nlay
		layer_nrne(nl)=layer_nrne(nl)+nadd
	enddo
	do nl=nl1+1,layer_nlay
		layer_nrna(nl)=layer_nrna(nl)+nadd
	enddo
	return
	end
c---------------------------------------------------------------------
c  print nodes to screen or file
c
	subroutine nodes_print(iunit,filename)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'layer.h'
	integer i,iunit
	character*(*) filename
	double precision z
c
	if(iunit.ne.6) then
		open(iunit,file=filename)
	endif
	write(iunit,'(a20)') ' RADIAL NODES  '
	write(iunit,'(a5,a15,a12)') 'nr','Radius (km)','Depth (m)'
	write(iunit,'(8i8)') nnod,jsd,jsu,jed,jeu,jwd,jwu,jstu
	do i=nnod,1,-1
		z=(layer_rb(layer_nlay)-rnod(i))*1.d3
		if(jsd.eq.jed) then
			if(i.eq.jsd) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' S- R-'
			else if(i.eq.jsu) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' S+ R+'
			else
				write(iunit,'(i5,f15.6,f12.3)') i,rnod(i),z
			endif
		else
			if(i.eq.jsd) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' S-'
			else if(i.eq.jsu) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' S+'
			else if(i.eq.jed) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' R-'
			else if(i.eq.jeu) then
				write(iunit,'(i5,f15.6,f12.3,a6)') i,rnod(i),z,' R+'
			else 
				write(iunit,'(i5,f15.6,f12.3)') i,rnod(i),z
			endif
		endif
	enddo
	write(iunit,'(a6,a15,2a6)') 'Layer','Depth in m','nrna','nrne'
	do i=1,layer_nlay
		write(iunit,'(i6,f15.2,2i6)') i,(layer_rb(layer_nlay)-layer_rb(i))*1.d3,
     1	                              layer_nrna(i),layer_nrne(i)
	enddo
	if(iunit.ne.6) close(iunit)
	return
	end
c---------------------------------------------------------------------
c  read nodes from file
c
	subroutine nodes_read(iunit,filename)
	include 'nodesdim.h'
	include 'nodes.h'
	integer i,iunit,k
	character*(*) filename
	character cdum*6
	real z
c
	open(iunit,file=filename)
	read(iunit,'(1x)') 
	read(iunit,'(1x)') 
	read(iunit,'(8i8)') nnod,jsd,jsu,jed,jeu,jwd,jwu,jstu
	if(nnod.gt.nnd) then
		print *,'<nodes_read>: More nodes than dimensioned',nnod,nnd	
		stop
	endif
	do i=nnod,1,-1
		if(jsd.eq.jed) then
			if(i.eq.jsd) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else if(i.eq.jsu) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else
				read(iunit,'(i5,f15.6,f12.3)') k,rnod(i),z
			endif
		else
			if(i.eq.jsd) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else if(i.eq.jsu) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else if(i.eq.jed) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else if(i.eq.jeu) then
				read(iunit,'(i5,f15.6,f12.3,a6)') k,rnod(i),z,cdum
			else 
				read(iunit,'(i5,f15.6,f12.3)') k,rnod(i),z
			endif
		endif
	enddo
	close(iunit)
	return
	end
c----------------------------------------------------------------------------
c  interpolate stored function vector at nodes in a given layer
c
	subroutine nodes_getsolution(nl,nvar,x1,x2,idir)
	include 'layer.h'
	include 'nvmax.h'
	include 'bsstore.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	include 'zero.h'
	integer i,nl,k,nr,nvar,idir,firsttime,klo,inc,khi
	double precision x1,x2,a,b,h
c
c  no node in layer
c  this may happen for small starting radii
c
	if(layer_nrna(nl).eq.0) return
c
c  Loop over nodes in layer
c
	firsttime=1
	do nr=layer_nrna(nl),layer_nrne(nl)
c
c  no function value available, do nothing,
c  because either x1 is above node (idir>0) or x1 is below node (idir<0)
c
		if(idir*(x1-rnod(nr)).gt.zero) then
			goto 1
c
c  no function value available, do nothing,
c  because either x2 is below node (idir>0) or x2 is above node (idir<0)
c
		else if(idir*(rnod(nr)-x2).gt.zero) then
			goto 1
c
c  normal case
c  for ascending  xstore we might get k=0 if rnod(nr)=xstore(1)
c  for descending xstore we might get k=jstore if rnod(nr)=xstore(jstore)
c  Interpolation formula is correct for both directions because
c  signs of h and xstore(k+1)-rnod(nr) change.
c
		else
			if(firsttime.eq.1) then
				call dlocate(rnod(nr),jstore,xstore,k)
				if(idir.gt.0) then
					if(k.eq.0) k=1
				else
					if(k.eq.jstore) k=k-1
				endif
				firsttime=0
c
c  perform a hunt, starting with previously found k
c
			else
				inc=idir
				if(idir.gt.0) then
					klo=k
 11					khi=min(klo+inc,jstore)
					if(khi.lt.jstore.and.xstore(khi).lt.rnod(nr)) then
						klo=khi
						inc=2*inc 
						goto 11
					else
						call dlocate(rnod(nr),khi-klo+1,xstore(klo),k)
						if(k.eq.0) k=1
						k=k+klo-1
					endif
				else if (idir.lt.0) then
					khi=k+1
 12					klo=max(khi+inc,1)
					if(klo.gt.1.and.xstore(klo).lt.rnod(nr)) then
						khi=klo
						inc=2*inc 
						goto 12
					else
						call dlocate(rnod(nr),khi-klo+1,xstore(klo),k)
						if(k.eq.khi-klo+1) k=k-1
						k=k+klo-1
					endif
				endif
			endif
			h=xstore(k+1)-xstore(k)
			a=(xstore(k+1)-rnod(nr))/h
			b=1.d0-a
c			print *,'nodes_getsolution:',nr,k,xstore(k),rnod(nr),xstore(k+1),h,a,b,ystore(2,k),ystore(2,k+1)
			do i=1,nvar
				ynod(i,nr)=a*ystore(i,k)+b*ystore(i,k+1)
			enddo
c			print *,(ynod(i,nr),i=1,4)
		endif
 1		continue
	enddo
c
	return
	end
c-------------------------------------------------------
c  find index of DSV node for given radius
c  take node above given radius
c  if selected node is a double node, take the upper one
c
	subroutine nodes_getnodeindex(r,n)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'zero.h'
	integer n,k
	double precision r
c
	call dlocate(r,nnod,rnod,k)
	if(k.eq.0) then
		n=1
	else if(k.eq.nnod) then
		n=nnod
	else
		n=k+1
		if(n.lt.nnod .and. dabs(rnod(n)-rnod(n+1)).lt.zero) n=n+1
	endif
	return
	end
