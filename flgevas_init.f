c------------------------------------------------------------------------
c  $Id: flgevas_init.f,v 1.6 2003/06/01 07:44:50 wolle Exp $
c
c  Initialization of flgevas:
c  Read model, create layer information, print it,
c  make radial nodes using given frequency, and print node info.
c
c  cmode:			name of model file (see flnm_readiso)
c  f:			      frequency for which radial nodes are generated
c  zstore:			depth down to which nodes are generated
c  zs:			source depth
c  ze:			receiver depth
c  nnmin:			minimum number of radial nodes per layer
c  iprint:			print flag (if > 0 you get a lot of output for debugging)
c------------------------------------------------------------------------
	subroutine flgevas_init(cmodel,f,zstore,zs,ze,nnmin,ncomp,comps,npar,pars,
     1	                  ncompsph,compssph,ncomptor,compstor,npartor,parstor,
     1				istyp,nsrcsph,jump_spher,nsrctor,jump_tor,iprint)
	include 'flnm.h'
	integer nnmin,iprint,nlay,i,j
	integer ncompsph,compssph(ncomp),ncomptor,compstor(ncomp),
     1	        npartor,parstor(npar),ncomp,comps(ncomp),npar,pars(npar),
     1	        nsrcsph,nsrctor,jump_spher(4),jump_tor(2),istyp
	double precision f,zstore,zs,ze,r1,r2
	character*(*) cmodel
c
	call flnm_read(1,cmodel)
c
c  Note: the earth model is for the reference frequency. If dispersion
c        is desired, elastic parameters should be recomputed for
c        each frequency. This requires an additional routine.
c
	call earthmodel_fromflnm(flnm_fref,0)
	call earthmodel_layer
	call layer_printinfo
	call earthmodel_spline
	call layer_set_stepsize_wlfrac(0.25d0)
	zs=max(0.001d0,zs)*1.d-3
	ze=max(0.001d0,ze)*1.d-3
c
c  zstore is set here to depth of halfspace+1/4 of the thickness of the overburden.
c  The code assumes that zstore is in the halfspace !!! (nlstore=1)
c
	call layer_getrb(1,r1)
	call layer_getnlay(nlay)
	call layer_getrb(nlay,r2)
	zstore=(r2-r1)*1.25
c
c  create Green nodes
c
	call nodes_create(f,4.d0,zstore,zs,ze,nnmin,iprint)
	if(iprint.eq.1) call nodes_print(6,'hallo')
c
c  separate components into spheroidal and toroidal parts
c  to split computations
c
	ncompsph=0
	ncomptor=0
	do i=1,ncomp
		if(comps(i).le.4) then
			ncompsph=ncompsph+1
			compssph(ncompsph)=comps(i)
		else if(comps(i).gt.4) then
			ncomptor=ncomptor+1
			compstor(ncomptor)=comps(i)-4
		endif
	enddo
c
c  separate parameters into spheroidal and toroidal
c
	npartor=0
	do j=1,npar
		if(pars(j).eq.1.or.pars(j).eq.3.or.pars(j).eq.5) then
			npartor=npartor+1
			parstor(npartor)=pars(j)
		endif
	enddo
c
c  specify for which jump vectors solutions are computed
c  initialize flags with zero
c
c  Explosion
c
	if(istyp.eq.1) then
		nsrcsph = 2
		nsrctor = 0
		jump_spher(1) = 1
		jump_spher(2) = 2
c
c  Vertical single force
c
	else if(istyp.eq.2) then
		nsrcsph = 1
		nsrctor = 0
		jump_spher(1) = 1
c
c  Horizontal single force
c
	else if(istyp.eq.3) then
		nsrcsph = 1
		nsrctor = 1
		jump_spher(1) = 2
		jump_tor(1) = 1
c
c  General single force
c
	else if(istyp.eq.4) then
		nsrcsph = 2
		nsrctor = 1
		jump_spher(1) = 1
		jump_spher(2) = 2
		jump_tor(1) = 1
c
c  General moment tensor source
c
	else if(istyp.eq.5) then
		nsrcsph = 4
		do i=1,4
			jump_spher(i) = i
		enddo
		nsrctor = 2
		do i=1,2
			jump_tor(i) = i
		enddo
	endif
c
	return
	end

