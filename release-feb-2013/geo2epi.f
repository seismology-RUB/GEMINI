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
c--------------------------------------------------------------
c  computes epicentral coordinates of a point on the unit
c  sphere with respect to pole thpol,phipol that has 
c  geographical coordinates thgeo, phigeo
c  returns theta between 0 and pi
c  and phi between 0 and 2*pi
c--------------------------------------------------------------
	subroutine geo2epi(thgeo,phigeo,thpol,phipol,theta,phi)
	real theta,phi,thpol,phipol,thgeo,phigeo
	double precision cthpol,sthpol,ctheta,stheta,cthgeo,sthgeo
	double precision pi,sa,ca,alfa,one
c
	pi=4.d0*datan(1.d0)
	one=1.d0
	cthpol=dcos(dble(thpol))
	sthpol=dsin(dble(thpol))
	cthgeo=dcos(dble(thgeo))
	sthgeo=dsin(dble(thgeo))
	ctheta=cthpol*cthgeo+sthpol*sthgeo*dcos(dble(phigeo)-dble(phipol))
	if(abs(ctheta).gt.1.d0) ctheta=sign(1.d0,ctheta)
	theta=dacos(ctheta)
	if(theta.lt.1.d-4*pi.or.theta.gt.0.999d0*pi) then
		phi=phipol
		return
	endif
	stheta=dsin(dble(theta))
	sa=sthgeo*dsin(dble(phigeo)-dble(phipol))/stheta
	ca=(cthgeo-ctheta*cthpol)/(stheta*sthpol)
	if(abs(sa).gt.1.d0) then
		sa=sign(one,sa)
	endif
	alfa=dasin(sa)
	if(ca.ge.0.d0) phi=alfa
	if(ca.lt.0.d0) then
		if(sa.ge.0.d0) phi=pi-alfa
		if(sa.lt.0.d0) phi=-pi-alfa
	endif
	phi=pi-phi
c
	return
	end

