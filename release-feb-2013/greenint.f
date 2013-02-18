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
c-----------------------------------------------------------------------------------
c                             GRREENINT
c
c  compute a Green function solution for given source terms and available
c  minors using Woodhouse's method
c
c  nls:			layer index of source
c  zsrc:			jump vector at source
c  minors:			minors(1-6,below-above-rs)
c  mtil:			mtilde-matrix(1-4,1-4,below-above-rs,receiver-node)
c  det:			determinant formed from minors at source
c  yout:			output Green 4-vector(1-4,receiver-node)
c  jd,ju:			rnode index of rs- and rs+, or re- and re+
c  jstu:			first receiver above starting radius
c  sflag:			if set, get a solution for source at rs, else for source at re
c  iprint:              print level	
c------------------------------------------------------------------------------------
	subroutine greenint(nls,zsrc,minors,mtil,det,yout,sflag,iprint)
	include 'nvmax.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	double complex zsrc(4),minors(6,2),mtil(4,4,2,nnd),det,deth,yout(4,nnd)
	double complex ystart(nvmax),fac
	double precision x1,x2
	integer nls,sflag,ifls,ifla,nlay
	integer jend,nr,i,j,nra,nvar,iprint,jd,ju
	character csyst*1
	external syssolidt
c
	if(sflag.eq.1) then
		jd=jsd
		ju=jsu
	else
		jd=jed
		ju=jeu
	endif
c
	nvar=4
	csyst='P'
c
c  zero yout
c
	do nr=1,nnod
		do i=1,4
			yout(i,nr)=(0.d0,0.d0)
		enddo
	enddo
c
	call layer_isfluid(nls,ifls)
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifla)
c
c   Upward integration from rs to nodes. If rs is in the water layer, we can directly
c   assign vlalues there. If not, we integrate ODEPCT to either
c   the surface or the bottom of the water layer and the make the transition to the
c   water layer if present.
c 
	if(iprint.gt.0) then
		print *,'              GREENINT above rs/re'
		write(6,'(a5,4a13)') 'nr','y1','y2','y3','y4'
	endif
c
c  source is in fluid layer
c	
	if(ifls.eq.1) then
		if(iprint.gt.0) print *,'in water'
		do nr=ju,nnod
			yout(1,nr)=mtil(1,2,2,nr)*(-zsrc(1)*minors(2,1)+zsrc(2)*minors(1,1))/det
			yout(2,nr)=mtil(1,3,2,nr)*(-zsrc(1)*minors(2,1)+zsrc(2)*minors(1,1))/det
			if(iprint.gt.0) write(6,'(i5,5e13.3)') nr,(zabs(yout(i,nr)),i=1,nvar)
		enddo
c
c  source is in solid layer
c
	else
		ystart(1)=-minors(6,1)*zsrc(2)+minors(5,1)*zsrc(3)-minors(4,1)*zsrc(4)
		ystart(2)=+minors(6,1)*zsrc(1)-minors(3,1)*zsrc(3)+minors(2,1)*zsrc(4)
		ystart(3)=-minors(5,1)*zsrc(1)+minors(3,1)*zsrc(2)-minors(1,1)*zsrc(4)
		ystart(4)=+minors(4,1)*zsrc(1)-minors(2,1)*zsrc(2)+minors(1,1)*zsrc(3)
		x1=rnod(ju)
		x2=rnod(nnod)
		jend=nnod
		if(ifla.eq.1) then
			x2=rnod(jwd)
			jend=jwd
		endif
		call propag(x1,x2,ystart,nvar,csyst,syssolidt)
		if(iprint.gt.0) print *,'in solid'
		do nr=ju,jend
			do i=1,4
				yout(i,nr)=0.d0
				do j=1,4
					yout(i,nr)=yout(i,nr)+mtil(i,j,2,nr)*ynod(j,nr)
				enddo   
				yout(i,nr)=+yout(i,nr)/det
			enddo
			if(iprint.gt.0) write(6,'(i5,5e13.5)') nr,(zabs(yout(i,nr)),i=1,nvar)
		enddo
c
c  Transition from solid to water layer
c
		if(ifla.eq.1) then			
			if(zabs(mtil(1,1,2,jwu)).gt.zabs(mtil(2,1,2,jwu))) then
				fac=yout(1,jwd)/mtil(1,2,2,jwu)
			else
				fac=yout(2,jwd)/mtil(1,3,2,jwu)
			endif
			if(iprint.gt.0) print *,'in water'
			do nr=jwu,nnod
				yout(1,nr)=fac*mtil(1,2,2,nr)
				yout(2,nr)=fac*mtil(1,3,2,nr)
			enddo
		endif
	endif
	if(iprint.gt.1) print *,'yout(i,jend): ',(yout(i,nnod),i=1,4)
c
c  Downward integration from rs to nodes below the source. If rs is in the water layer
c  we first step through this layer and the continue with ODEPCT in the solid layers
c  by using the appropriate starting values.
c
	if(iprint.gt.0) then
		print *,'           GREENINT below rs/re'
		write(6,'(a5,4a13)') 'nr','y1','y2','y3','y4'
	endif
	if(ifls.eq.1) then
		if(iprint.gt.0) print *,'in water'
		do nr=jd,jwu,-1
			yout(1,nr)=mtil(1,2,1,nr)*(-zsrc(1)*minors(2,2)+zsrc(2)*minors(1,2))/det
			yout(2,nr)=mtil(1,3,1,nr)*(-zsrc(1)*minors(2,2)+zsrc(2)*minors(1,2))/det
			if(iprint.gt.0) write(6,'(i5,5e13.5)') nr,(zabs(yout(i,nr)),i=1,nvar)
		enddo
c
c  initial values for transition from water to solid
c
		ystart(1)=-yout(2,jwu)
		ystart(2)= yout(1,jwu)
		ystart(3)=0.d0
		ystart(4)=0.d0
		nra=jwd
		deth=-mtil(1,2,1,jwd)
	else
		ystart(1)=-minors(6,2)*zsrc(2)+minors(5,2)*zsrc(3)-minors(4,2)*zsrc(4)
		ystart(2)=+minors(6,2)*zsrc(1)-minors(3,2)*zsrc(3)+minors(2,2)*zsrc(4)
		ystart(3)=-minors(5,2)*zsrc(1)+minors(3,2)*zsrc(2)-minors(1,2)*zsrc(4)
		ystart(4)=+minors(4,2)*zsrc(1)-minors(2,2)*zsrc(2)+minors(1,2)*zsrc(3)
		nra=jd
		deth=det
	endif
	if(iprint.gt.0) print *,'in solid'
	x1=rnod(nra)
	x2=rnod(jstu)
	call propag(x1,x2,ystart,nvar,csyst,syssolidt)
	do nr=nra,jstu,-1
c		print *,nr,(zabs(ynod(i,nr)),i=1,nvar)
c		print *,nr,(zabs(mtil(2,i,1,nr)),i=1,4)
		do i=1,4
			yout(i,nr)=0.d0
			do j=1,4
				yout(i,nr)=yout(i,nr)+mtil(i,j,1,nr)*ynod(j,nr)
			enddo   
			yout(i,nr)=-yout(i,nr)/deth
		enddo
		if(iprint.gt.0) write(6,'(i5,5e13.5)') nr,(zabs(yout(i,nr)),i=1,nvar)
	enddo
c
	return
	end
