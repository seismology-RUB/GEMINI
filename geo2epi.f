c--------------------------------------------------------------
c  computes epicentral coordinates of a point on the unit
c  sphere with respect to pole thpol,phipol that has 
c  geographical coordinates thgeo, phigeo
c  returns theta between 0 and pi
c  and phi between 0 and 2*pi
c--------------------------------------------------------------
	subroutine geo2epi(thgeo,phigeo,thpol,phipol,theta,phi)
	real theta,phi,thpol,phipol,thgeo,phigeo
	real cthpol,sthpol,ctheta,stheta,cthgeo,sthgeo
	real pi,sa,ca,alfa,one
c
	pi=4.*atan(1.)
	one=1.
	cthpol=cos(thpol)
	sthpol=sin(thpol)
	cthgeo=cos(thgeo)
	sthgeo=sin(thgeo)
	ctheta=cthpol*cthgeo+sthpol*sthgeo*cos(phigeo-phipol)
	if(abs(ctheta).gt.1.) ctheta=sign(1.,ctheta)
	theta=acos(ctheta)
	if(theta.lt.1.e-4*pi.or.theta.gt.0.999*pi) then
		phi=phipol
		return
	endif
	stheta=sin(theta)
	sa=sthgeo*sin(phigeo-phipol)/stheta
	ca=(cthgeo-ctheta*cthpol)/(stheta*sthpol)
	if(abs(sa).gt.1.0) then
		sa=sign(one,sa)
	endif
	alfa=asin(sa)
	if(ca.ge.0.) phi=alfa
	if(ca.lt.0.) then
		if(sa.ge.0.) phi=pi-alfa
		if(sa.lt.0.) phi=-pi-alfa
	endif
	phi=pi-phi
c
	return
	end

