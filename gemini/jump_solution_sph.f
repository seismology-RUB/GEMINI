c----------------------------------------------------------------------
c  compute DSV jump solutions for sources at radial nodes and receiver at re
c
	subroutine jump_solution_sph(wn,f,sigma,ze,gfe,nsrc,zsrc,iprint,alternative,derivflag)
	include 'pi.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	double precision wn,f,sigma,ze
	integer iprint,nsrc,derivflag
	double precision om,elp1,re,rst,rearth
	double complex zom
	double complex det(2),minors(6,2,2),mtil(4,4,2,nnd)
	double complex gfe(nnd,6,4),gh(4,nnd),ghe(4,nnd),zsrc(4,4,nnd),zsh(4)
	integer nle,nlstrt,nr,i,j,top,nb
	integer nb4comp(4)
	logical alternative
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
c
c  get layer of re
c
	call layer_getindex(re,nle,top)
c
c  zero gfe
c
	do j=1,nsrc
		do i=1,6
			do nr=1,nnod
				gfe(nr,i,j)=0.d0
			enddo
		enddo
	enddo
c
c  return if wn=0, gfe = 0
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
c  if starting radius is above receiver
c  stop working , because motion is observed for this k,om.
c  The starting radius may not imply a zero solution there if
c  it is in the halfspace. But startr.f ensures that the starting radius
c  is below the deepest node in that case which may also be source or receiver.
c  Thus, if rst happens to be above re it does imply a zero solution there
c  and we can safely skip the integration.

	if (rst.ge.re) return
c
c  From now on we know that the starting radius is below receiver c
c
c  integrate minors from starting radius up to max(re,re) and down from
c  surface to min(re,re)

	call minorint(nle,re,nle,re,nlstrt,rst,minors,mtil,det,iprint)
c
c  Integrate here four basis DSV vectors for a source at re.
c
	if (alternative) then
		do j=1,nsrc
			do i = 1,4
				zsh(i) = zsrc(i,j,1)
			enddo
			call greenint(nle,zsh,minors(1,1,2),mtil,det(2),gh,0,iprint)
			do i = 1,4
				do nr = jstu,nnod
					gfe(nr,i,j) = gh(i,nr)
				enddo
			enddo
			if (derivflag == 1) call sphderivs(gfe(1,1,j))
		enddo
c
c  The reciprocity principle is used here to determine from
c  these 4 basis DSV vectors the desired DSV components at re 
c  for sources at the r-nodes.
c  ghe: first index is basis solution
c
	else
		do i=1,4
			do nb=1,4
				zsh(nb)=0.d0
			enddo
			nb=nb4comp(i)
			zsh(nb)=1.d0
			call greenint(nle,zsh,minors(1,1,2),mtil,det(2),gh,0,iprint)
			do nr=jstu,nnod
				call reciprocity(nb,rnod(nr),re,elp1,gh(1,nr),ghe(1,nr))
			enddo
c
c  linearly combine basis solutions using prescribed jumps at source
c  nb: basis solution, jump vector = (1,0,0,0) etc
c  j:  different jump vectors needed for calculating seismograms for force or moment tensor sources
c  i:  DSV component (U,R,V,S)
c
			do j=1,nsrc
				do nr=jstu,nnod
					gfe(nr,i,j)=0.d0
					do nb=1,4
						gfe(nr,i,j)=gfe(nr,i,j)+zsrc(nb,j,nr)*ghe(nb,nr)
					enddo
				enddo
			enddo
		enddo
	endif
	return
	end

