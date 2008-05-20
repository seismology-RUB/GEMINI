c-------------------------------------------------------------------------------------
c  $Id: greenkfsph.f,v 1.2 2003/06/01 07:44:50 wolle Exp $
c
c  SPHEROIDAL MOTION 
c
c  Compute depth dependent Green functions
c  for given values of frequency and wavenumber
c
c  wn:		wavenumber in 1/km
c  f:			frequency in Hz
c  sigma:		negative imaginary part of omega
c  zs:		source depth in km
c  ze:		receiver depth in km
c  ncomp:		number of desired components
c  comps:		components array
c  npar:		number of model parameters for which Frechet kernels are computed
c  pars:		indices of model parameters for which Frechet kernels are computed
c			1 = density, 2 = P-velocity, 3 = S-velocity, 4 = qkappa, 5 = qmu
c  gfs:		gfs(4,nr,jump): 4-component Greens function at radial nodes nr
c		                      for a source at rs using specific jumps at zs
c  gfsr:		gfsr(4,jump): 4-component Greens function at receiver ze
c  gfe:		gfe(4,comp,nr):  4 basis Greens functions at receiver re for
c		                       given components comp and sources at the radial nodes
c  sk:		sk(nr,comp,par,jump): Hat kernel for each component, earth model node, parameter, jump
c  skd:		skd(disco,comp,jump): Kernel for perturbation of layer boundaries
c			Here nkk is used as dimension instead of nlayer to avoid including layer.h
c  selflag:		a selection flag that tells greenkfsph what to do.
c			The following options are available:
c			(1) compute both gfs and gfe (needed for partial derivatives)
c			(2) compute only gfs (if a Green function for a given source
c			    at various receivers depths is needed)
c  nsrc:       	number of sources (max 4)
c  zsrc:       	jump at the source
c  jump:          index of jump basis vector associated with zsrc (see gemini-displ.tex)
c  iprint:		print flag
c--------------------------------------------------------------------------------
c  Note: only request the components you really need at re because there will be
c        one integration per component
c
c  Theory outlined in gemini-displ.tex, rezifrechet.tex and flgevask.tex
c-------------------------------------------------------------------------------------- 
	subroutine greenkfsph(wn,f,sigma,zs,ze,ncomp,comps,npar,pars,
     1	                   gfsr,sk,selflag,nsrc,zsrc,jump,iprint)
	include 'pi.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	include 'parameter.h'
	integer ncomp,comps(ncomp),npar,pars(nparr),jump(nsrc)
	integer top,nls,nle,jstu,ico,i,nb,selflag,iprint
	integer nr,ipar,nlay,nsrc,nlstrt,m
	integer nb4comp(4)
	double precision wn,f,zs,ze,sigma
	double precision om,elp1,rst,rs,re,rearth
	double complex det(2),gh(4,nnd),zsrc(4,4),gfs(4,nnd,4)
	double complex zsh(4),sk(nkk,4,nparr,4),skd(nkk,4,4)
	double complex gfe(4,4,nnd),minors(6,2,2),mtil(4,4,2,nnd),pdev(nnd,4),gfsr(4,4)
	double complex zom
	common/degree/elp1
	common/omega/zom,om
	data nb4comp/2,1,4,3/
c
	if(iprint.gt.0) write(6,'(a,f10.3)') ' wn/(2*pi) = ',wn/(2.*pi)
	om=2.*pi*f
	zom=dcmplx(om,-sigma)
	rearth=earthmodel_rearth
	elp1=(rearth*wn)**2
	re=rearth-ze
	rs=rearth-zs
c
c  get layer of rs and re
c
	call layer_getindex(rs,nls,top)
	call layer_getindex(re,nle,top)
c
c  zero gfs and gfe and gfsr and sk and skd
c
	do m=1,nsrc
		do i=1,4
			gfsr(i,m)=0.d0
			do nr=1,earthmodel_nk
				skd(nr,comps(i),m)=0.d0
			enddo
		enddo
		do ipar=1,npar
			do i=1,4
				do nr=1,earthmodel_nk
					sk(nr,i,ipar,m)=0.d0
				enddo
			enddo
		enddo
	enddo
	do nr=1,nnod
		do i=1,4
			do m=1,nsrc
				gfs(i,nr,m)=0.d0
			enddo
			do nb=1,4
				gfe(nb,i,nr)=0.d0
			enddo
		enddo
	enddo
c
c  return if wn=0, gfs, gfe and pdev = 0
c
	if(wn.lt.1.d-8) return
c
c  compute starting radius
c  for a global application, more help to find the starting radius
c  is needed (see lmspher)
c
	call startr(rnod(1),rst,nlstrt,iprint)
	if(iprint.gt.0) print *,'<greenkfsph>: rst,nlstrt :',rst,nlstrt
c
c  if starting radius is above the source or receiver
c  stop working , because motion is either not excited or
c  not observed for this k,om.
c  The starting radius may not imply a zero solution there if
c  it is in the halfspace. But startr.f ensures that the starting radius
c  is below the deepest node in that case which may also be source or receiver.
c  Thus, if rst happens to be above rs,re it does imply a zero solution there
c  and we can safely skip the integration.

	if (rst.ge.re.or.rst.ge.rs) return
c
c  From now on we know that the starting radius is below source AND receiver !
c
c  integrate minors from starting radius up to max(rs,re) and down from
c  surface to min(rs,re)

	call minorint(nls,rs,nle,re,nlstrt,rst,minors,mtil,det,iprint)
c
c  now get a 4-component Green function for a source at rs using source terms
c  by integrating up and downward from rs with Woodhouse's method
c
	if(selflag.eq.1.or.selflag.eq.2) then
		do m=1,nsrc
			call greenint(nls,zsrc(1,jump(m)),minors(1,1,1),mtil,det(1),gfs(1,1,jump(m)),1,iprint)
		enddo
		if(iprint.eq.1) print *,'gfs: ',(gfs(i,nnod-6,1),i=1,4)
	endif
c
c  Integrate here selected 4-component basis Green vectors 
c  for a source at re.
c  The reciprocity principle is used to determine from
c  these 4 basis Green vectors the desired Green components at re 
c  for sources at the r-nodes.
c
	if(selflag.eq.1) then
		do i=1,ncomp
			ico=comps(i)
			do nb=1,4
				zsh(nb)=0.d0
			enddo
			nb=nb4comp(ico)
			zsh(nb)=1.d0
			call greenint(nle,zsh,minors(1,1,2),mtil,det(2),gh,0,iprint)
			do nr=jstu,nnod
				call reciprocity(nb,rnod(nr),re,elp1,gh(1,nr),gfe(1,ico,nr))
			enddo
		enddo
	endif
c
c  frechet returns the integrand pdev = G(re,r) dA(r) y(r,rs)
c  and jump of G(re,r)*dA(r)*y(r,rs) at the interfaces
c  pdev=0 for nr < jstu
c
c  compute change of DSV for each node of the earth model. Assume a hat function
c  perturbation for each node with relative perturbation of 1. Let rk(k) be
c  the target earthmodel node. Then we integrate pdev(r)*hat(r) from rk(k-1)
c  to rk(k+1) and hat(rk(k-1))=hat(rk(k+1))=0. If rk(k) is a on discontinuity
c  a half hat is used on the appropriate side.
c  We return sk(model-node,comp,ipar,jump) containing the hat kernels.
c
c  TODO:
c  better split off computation of disco kernels because loop over ipar is 
c  not required for them
c
	if(selflag.eq.1) then
		call layer_getnlay(nlay)
		do m=1,nsrc
			do nb=1,npar
				ipar=pars(nb)
				call frechet(ipar,ncomp,comps,gfe,gfs(1,1,jump(m)),pdev,skd(1,1,jump(m)),iprint)
				call hatkernel(ncomp,comps,pdev,sk(1,1,ipar,jump(m)),iprint)
			enddo
		enddo
	endif
c
c  compute gfsr
c  gfs at other nodes was only needed for the hat kernels
c
	if(selflag.eq.1.or.selflag.eq.2) then
		do m=1,nsrc
			do i=1,4
				gfsr(i,jump(m))=gfs(i,jed,jump(m))
			enddo
		enddo
	endif
c
	return
	end
