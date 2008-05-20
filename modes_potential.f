c--------------------------------------------------------------------------
c  Compute displacement potential for spherically symmetric earth 
c  for one receiver.
c
c  iswitch:  1 = potential itself
c            2 = d/dtheta of potential
c            3 = 1/sin(theta) d/dphi of potential
c-------------------------------------------------------------------------
	subroutine potential(dis,phi,f,p,cg,qatt,psi,re,pot,iswitch)
	include 'pis.h'
c
	integer k,iswitch
	real pi,wu2,vl,vl2,cph,cg,eps,deltainv,re,arg,arg2,expo
	real dis,phi,f,p,cg,qatt
	real r1psph,r2psph,cosph,sinph,cos2ph,sin2ph,dfr1psph,dfr2psph
	complex psi(0:2),ci,han(0:2),feldfac0,feldfac,dhan(0:2),pot
	double precision delta,sdelta
	real bessj, bessy, bessj0, bessj1, bessy0, bessy1
c
c
	ci=cmplx(0.,1.)
	wu2=sqrt(2.)
c
	cph = 1./(re*p)
	eps=(pi*f)/(qatt*cg/re)
	feldfac0=-ci/(4.*cph*cg/re)
	vl=2.*pi*f/cph
	vl2=vl*vl
c
c  Compute reference wavefield of exciting mode at receivers
c  
	delta=dble(dis/re)
	deltainv=1./delta
	sdelta=dsin(delta)
	feldfac=feldfac0*dsqrt(delta/sdelta)
	arg=dble(vl)*delta
	expo=dble(eps)*delta
c
c  Hankel functions
c
	if(arg.lt.6.*pi) then
		han(0)=bessj0(arg)-ci*bessy0(arg)
		han(1)=bessj1(arg)-ci*bessy1(arg)
		han(2)=bessj(2,arg)-ci*bessy(2,arg)
	else 
		do k=0,2
			arg2=arg-.5*pi*(k+.5)
			han(k)=sqrt(2./(pi*arg))*cmplx(cos(arg2),-sin(arg2))*exp(-expo)
		enddo
	endif
c
	dhan(0)=-vl*han(1)
	dhan(1)=-deltainv*han(1)+vl*han(0)
	dhan(2)=-2.*deltainv*han(2)+vl*han(1)
c
	cosph=cos(phi)
	sinph=sin(phi)
	cos2ph=cos(2.*phi)
	sin2ph=sin(2.*phi)
	r1psph=-wu2*vl*( real(psi(1))*cosph-aimag(psi(1))*sinph )
	dfr1psph=-wu2*vl*( -real(psi(1))*sinph-aimag(psi(1))*cosph )
	r2psph=2.*vl2*( real(psi(2))*cos2ph-aimag(psi(2))*sin2ph )
	dfr2psph=2.*vl2*( -2.*real(psi(2))*sin2ph-2.*aimag(psi(2))*cos2ph )
	if(iswitch.eq.1) then
 		pot=feldfac*(han(0)*real(psi(0))+han(1)*r1psph+han(2)*r2psph)
	else if(iswitch.eq.2) then
		pot=feldfac*(dhan(0)*real(psi(0))+dhan(1)*r1psph+dhan(2)*r2psph)
	else if(iswitch.eq.3) then
		pot=feldfac/sdelta*(han(1)*dfr1psph+han(2)*dfr2psph)
	endif
c
	return
	end
