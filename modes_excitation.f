c---------------------------------------------------------
c  excitation coefficients for spheroidal modes
c
	subroutine excoso_sph(rm,rs,vl,usf,vsf,upsf,vpsf,psi)
c
	complex ci,psi(0:2)
	real rm(6),usf,vsf,upsf,vpsf,fs,xs,wu2,vl,vl2
	double precision rs
c
	ci=cmplx(0.,1.)
	wu2=sqrt(2.)
c
	vl2=vl*vl
	fs=(2.*usf-vl2*vsf)/rs
	xs=vpsf+(usf-vsf)/rs
	psi(0)=rm(1)*upsf+(rm(2)+rm(3))*.5*fs
	psi(1)=(-rm(4)+ci*rm(5))*xs/wu2
	psi(2)=(.5*(rm(2)-rm(3))-ci*rm(6))*.5*vsf/rs
c
	return
	end
c--------------------------------------------------------
c  excitation coefficients for toroidal modes
c
	subroutine excoso_tor(rm,rs,vl,wsf,wpsf,psi)
c
	complex ci,psi(0:2)
	real rm(6),wsf,wpsf,zs,wu2,vl,vl2
	double precision rs
c
	ci=cmplx(0.,1.)
	wu2=sqrt(2.)
c
	vl2=vl*vl
	zs=wpsf-wsf/rs
	psi(0)=0.
	psi(1)=(-rm(4)+ci*rm(5))*ci*zs/wu2
	psi(2)=(.5*(rm(2)-rm(3))-ci*rm(6))*.5*ci*wsf/rs
c
	return
	end
