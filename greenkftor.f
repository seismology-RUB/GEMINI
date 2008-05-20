c-------------------------------------------------------------------------------------
c  $Id:$
c
c  TOROIDAL MOTION
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
c  gfs:		gfs(2,nr,jump): 2-component Greens function at radial nodes nr
c		                      for a source at rs using specific jumps at rs
c  gfsr:		gfsr(2,jump): 2-component Greens function at receiver ze
c  gfe:		gfe(2,comp,nr): 2 basis Greens functions at receiver re for
c		                      given components comp and sources at the radial nodes
c  sk:		sk(nk,comp,par,jump): Hat kernel for each component, earth model node, parameter and jump
c  skd:		skd(disco,comp,jump): Kernel for perturbation of layer boundaries
c			Here nkk is used as dimension instead of nlayer to avoid including layer.h
c  selflag:		a selection flag that tells greenpf what to do.
c			The following options are available:
c			(1) compute both gfs and gfe (needed for partial derivatives)
c			(2) compute only gfs (if a Green function for a given source
c			    at various receivers depths is needed)
c  nsrc:       	number of sources (max 2)
c  zsrc:       	jump vector at the source
c  jump:          index of jump vector associated with zsrc (see gemini-displ.tex)
c  iprint:		print flag
c--------------------------------------------------------------------------------
c  Note: only request the components you really need at re because there will be
c        one integration per component
c
c  Theory: gemini-displ.tex, flgevask.tex, rezifrechet.tex 
c-------------------------------------------------------------------------------------
	subroutine greenkftor(wn,f,sigma,zs,ze,ncomp,comps,npar,pars,
     1	                gfsr,sk,selflag,nsrc,zsrc,jump,iprint)
	include 'pi.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	include 'parameter.h'
	double precision wn,f,zs,ze,sigma
	integer ncomp,comps(ncomp),npar,pars(nparr),jump(nsrc)
	integer top,nls,nle,jstu,i,nb,selflag,iprint,nr,m,mm
	integer ipar,nlay,nlstrt,nsrc,nnsol
	double precision om,elp1,rst,rs,re,rearth,ratio
	double complex det(2),gh(2,2,nnd),zsrc(2,2),gfs(2,nnd,2)
	double complex sk(nkk,2,nparr,2),skd(nkk,2,2)
	double complex gfe(2,2,nnd),pdev(nnd,2),gfsr(2,2),grse11,grse12,grse21,grse22
	double complex zom
	common/degree/elp1
	common/omega/zom,om
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
		do i=1,ncomp
			gfsr(comps(i),m)=0.d0
			do nr=1,earthmodel_nk
				skd(nr,comps(i),m)=0.d0
			enddo
		enddo
		do ipar=1,npar
			do nr=1,earthmodel_nk
				do i=1,ncomp
					sk(nr,comps(i),pars(ipar),m)=0.d0
				enddo
			enddo
		enddo
	enddo
	do nr=1,nnod
		do i=1,2
			do m=1,nsrc
				gfs(i,nr,m)=0.d0
			enddo
			do nb=1,2
				gfe(nb,i,nr)=0.d0
			enddo
		enddo
	enddo
c
c  return if wn=0, gfs, gfe and pdev = 0
c
	if(wn.lt.1.d-8) return
c
c  return if source or receiver are in a top fluid layer
c
	call layer_isfluid(nls,i)
	if(i.eq.1) return
	call layer_isfluid(nle,i)
	if(i.eq.1) return
c
c  compute starting radius
c  for a global application, more help to find the starting radius
c  is needed (see lmspher)
c
	call startr(rnod(1),rst,nlstrt,iprint)
	if(iprint.gt.0) print *,'<greenkftor>: rst,nlstrt :',rst,nlstrt
c
c  if starting radius is above the source or receiver
c  stop working , because motion is either not excited or
c  not observed for this k,om.
c  The starting radius may not imply a zero solution there if
c  it is in the halfspace. But startr.f ensures that the starting radius
c  is below the deepest node in that case which may also be source or receiver.
c  Thus, if rst happens to be above rs,re it does imply a zero solution there
c  and we can safely skip the integration.
c
	if (rst.ge.re.or.rst.ge.rs) return
c
c  return if starting radius is in a top fluid layer
c
	call layer_isfluid(nlstrt,i)
	if(i.eq.1) return
c
c  From now on we know that the starting radius is below source AND receiver !
c
c  integrate toroidal equation from starting radius up to rs,re and down from
c  surface to rs,re. nnsol = jwd if a water layer exists
c
	call torint(nls,rs,nle,re,nlstrt,rst,gh,det,nnsol,iprint)
c
c  construct a solution for a source at rs using the source terms
c  first below the source
c  Note: sign in Jörgs Diss is wrong, sign here is correct
c  second index in gh denotes (1) below and (2) above source
c  first index in gh is component (W, T)
c
	do m=1,nsrc
		mm=jump(m)
		do nr=jstu,jsd
			do i=1,2
				gfs(i,nr,mm)=(-zsrc(1,mm)*gh(2,2,jsu)+zsrc(2,mm)*gh(1,2,jsu))*gh(i,1,nr)/det(1)
			enddo
		enddo
c
c  above the source
c
		do nr=jsu,nnsol
			do i=1,2
				gfs(i,nr,mm)=(-zsrc(1,mm)*gh(2,1,jsd)+zsrc(2,mm)*gh(1,1,jsd))*gh(i,2,nr)/det(1)
			enddo
		enddo
	enddo
c
c  Construct two basis solutions for a source at rnod(nr) and receiver at re
c  using the reciprocity relation
c  first index: basis solution, second index: component
c  below receiver
c
	do nr=jstu,jed
		ratio=(rnod(nr)/re)**2
		grse11=-gh(2,2,jeu)/det(2)*gh(1,1,nr)
		grse12=-gh(2,2,jeu)/det(2)*gh(2,1,nr)
		grse21=+gh(1,2,jeu)/det(2)*gh(1,1,nr)
		grse22=+gh(1,2,jeu)/det(2)*gh(2,1,nr)
		gfe(2,2,nr)=-ratio*grse11
		gfe(2,1,nr)=+ratio*grse21
		gfe(1,2,nr)=+ratio*grse12
		gfe(1,1,nr)=-ratio*grse22
	enddo
c
c  above receiver
c
	do nr=jeu,nnsol
		ratio=(rnod(nr)/re)**2
		grse11=-gh(2,1,jed)/det(2)*gh(1,2,nr)
		grse12=-gh(2,1,jed)/det(2)*gh(2,2,nr)
		grse21=+gh(1,1,jed)/det(2)*gh(1,2,nr)
		grse22=+gh(1,1,jed)/det(2)*gh(2,2,nr)
		gfe(2,2,nr)=-ratio*grse11
		gfe(2,1,nr)=+ratio*grse21
		gfe(1,2,nr)=+ratio*grse12
		gfe(1,1,nr)=-ratio*grse22
	enddo
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
c  We return sk(model-node,comp,ipar) containing the hat kernels.
c
c  TODO
c  better split off computation of skd because loop over ipar is not required
c  for them
c
	if(selflag.eq.1) then
		call layer_getnlay(nlay)
		do m=1,nsrc
			do nb=1,npar
				ipar=pars(nb)
				call frechettor(ipar,ncomp,comps,gfe,gfs(1,1,jump(m)),pdev,skd(1,1,jump(m)),iprint)
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
			do i=1,2
				gfsr(i,jump(m))=gfs(i,jed,jump(m))
			enddo
		enddo
	endif
c
	return
	end
