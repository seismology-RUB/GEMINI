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
c  Find the starting radius.
c  This routine builds on the convention that
c  (1) the earthmodel is parameterized with nodes
c      and linear dependence in between is a 
c      good approximation and
c  (2) that the earth is homogeneous below
c      the deepest node with the properties
c      specified for the deepest node.
c  In that case we may look for the turning point
c  by starting at the top of the halfspace and
c  searching for a node where F=elp1-(om*r/v)**2
c  is negative.
c
c  Let r+,v+ be radius and velocity of the negative node
c  and r-,v- of the positive node below. Then F goes through
c  zero at r0=b/(om/elp1 - m) where
c  m = (v+ - v-)/(r+ - r-)   and
c  b=(r+v- - r-v+)/(r+ - r-).
c
c  The starting radius is found by integrating F(r)
c  downward from r0 until the integral gets greater
c  than vertno or the starting radius drops below rstore.
c
c  Special treatment if turning point is in halfspace  
c------------------------------------------------------------
	subroutine startr(rn1,rstart,nla,iprint)
	include 'earthmodel.h'
	include 'layer.h'
	include 'pi.h'
	include 'zero.h'
	double precision vertno,third,delj
	parameter(vertno=14.d0,third=0.333333333333d0,delj=1.0d0)
	double precision om,elp1,vel,rstart,q,p,pmin
	double precision rn1,r,dr,rtop,rbot,alf
	double precision qinttop,qintbot
	double precision rold,rturn,velold,slope,vcept,h
	integer j,nla,jmin,iktop1,jturn,top,nl,iprint
	double complex zom
	common/degree/elp1
	common/omega/zom,om
c
	call layer_getiktop(1,iktop1)
	j=iktop1
	pmin=1.d200
	jmin=j
 1	r=earthmodel_rk(j)
	vel=dsqrt( real(earthmodel_zmue(j))/earthmodel_ro(j) )
!	print *,'startr: vel = ',vel
	if(vel.lt.zero) vel=dsqrt( real(earthmodel_zkap(j))/earthmodel_ro(j) )
	p=elp1-(om*r/vel)**2
	if(p.lt.pmin) then
		pmin=p
		jmin=j
	endif
	if(iprint.gt.1) print *,'startr: ',j,p,r,rold
c
c  p is negative
c
	if(p.lt.0.d0) then		
c
c  turning point is in halfspace, set rstart,nla and return
c  If rn1 < rturn linearize r around rturn and
c  simplify integrand to sqrt(elp1)/rturn*sqrt(2u/rturn) with u=rturn-r.
c  This can be integrated analytically. The distance d below rturn where
c  the integral exceeds vertno is
c  d = ((3/(2*sqrt(2))*(1/sqrt(elp1))*vertno)**3/2*rturn
c
		if(j.eq.iktop1) then
c			print *,'WARNING form startr: turning point is in half space'
c			print *,'Increase the minimum slowness to slightly > 1/vs of the halfspace'
			rturn=dsqrt(elp1)/om*vel
			if(rn1.ge.r) then
				rstart=r-1.d-6
				nla=1
				goto 99
			else if(rn1.lt.r.and.rn1.ge.rturn) then
				rstart=rn1-1.d-6
				nla=1
				goto 99
			else if(rn1.lt.rturn) then
				rstart=rturn*(1.d0-(1.5d0/dsqrt(2.d0*elp1)*vertno)**1.5)
				if(rstart.le.rn1) rstart=rn1-1.d-6
				nla=1
				goto 99
			endif
c
c  turning point above halfspace, set rturn and goto integration of p
c  jturn gives upper node of interval containing rturn
c
		else
			h=r-rold
c
c  check if sign change occurs over discontinuity
c
			if(h.lt.1.d-12) then
				rturn=r
				jturn=j-1
				dr=delj*r/dsqrt(elp1-(om*r/velold)**2)
				qinttop=1.d0/dr
			else
				slope=(vel-velold)/h
				vcept=(r*velold-rold*vel)/h
				rturn=vcept/(om/dsqrt(elp1)-slope)
				jturn=j
				dr=(2.d0*delj**2*vel/elp1/vcept)**third*rturn
				qinttop=0.d0
			endif
			goto 2
		endif
c
c  p is still positive, next node
c
	else
c
c  reached the surface and no turning point found
c  take radius of minimum p as turning point
c
		if(j.eq.earthmodel_nk) then
			rturn=earthmodel_rk(jmin)
			jturn=jmin
			dr=0.5*(r-rold)
			qinttop=1.d0/dr
			goto 2
		else
			rold=r
			velold=vel
			j=j+1
			goto 1
		endif
	endif
c
c  integrate 1/r*dsqrt(p(r)) downwards with given dr until integral exceeds vertno
c
 2	continue
	q=0
	j=jturn
	rtop=rturn
	call layer_getindex(rtop,nl,top)
 3	rbot=rtop-dr
	if(rbot.lt.layer_rb(nl-1)) rbot=layer_rb(nl-1)
	call earthmodel_isovel(nl,rbot,alf,vel)
	if(vel.lt.1.d-12) vel=alf
	qintbot=dsqrt(elp1-(om*rbot/vel)**2)/rbot
	q=q+0.5*dr*(qinttop+qintbot)
	if(iprint.gt.1) print *,'startr: ',j,nl,rtop,rbot,qinttop,qintbot,q,dr
	if(nl.eq.1.and.rn1.ge.layer_rb(1)) then
			rstart=layer_rb(1)-1.d-6
			nla=1
			goto 99
	endif
	if(nl.eq.1.and.rn1.lt.layer_rb(1).and.rbot.le.rn1) then
		rstart=rn1-1.d-6
		nla=1
		goto 99
	endif
	if(q.ge.vertno) then
		rstart=rbot
		nla=nl
		goto 99
	else
		rtop=rbot
		qinttop=qintbot
		if(dabs(rtop-layer_rb(nl-1)).lt.1.d-12) nl=nl-1
		goto 3
	endif
c
 99	continue
c
	return
	end
