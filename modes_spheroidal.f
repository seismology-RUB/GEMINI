c------------------------------------------------------------------
c  compute eigenvalues and eigenfunctions for fixed frequency
c
c  fs:         Frequency in Hz
c  cmin:       min phase velocity
c  cmax:       max phase velocity
c  nmin:       min number of depth nodes per layer
c  noeig:      if set eigenfunctions are not computed, only nov and p are
c              filled with meaningful values 
c  iprint:     print flag (0=no output, 1=normal output, 2=debug output, >2= even more
c  nrt:        number of roots found  
c  nov:        true overtone number
c  p:          slowness of mode in s/km = (el+0.5)/(om R)
c  cg:         group velocity of mode
c  qatt:       Q of mode
c  yeig:       eigenfunctions: (U,R,V,S,dU/dr,dV/dr)
c  fk..:       kernels for (density, P, S)
c---------------------------------------------------------------------------
	subroutine modes_spheroidal(f,cmin,cmax,wlratio,nmin,iprint,noeig,
     1	                      nrt,nov,p,cg,qatt,yeig,fkro,fkvp,fkvs,fkd)
	include 'nvmax.h'
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'pi.h'
	include 'laydim.h'
	integer iprint,knt,nval,n,nmin,js,nrt,noeig,nlstrt
	double precision f,cmin,cmax,wlratio,rstart
	double precision p(nrtt),cg(nrtt),qatt(nrtt),raylquo(nrtt)
	double precision yeig(nvx,nnd,nrtt)
	double precision fkro(nnd,nrtt),fkvp(nnd,nrtt),fkvs(nnd,nrtt),fkd(nlayer,nrtt)
	double precision elp1,om,omh,f,x1min,elmin,elmax
	double precision rearth,ystart(nvmax),sqelp1,anorm(nrtt)
	double precision m5root(nrtt),x1s(nrtt),elroot(nrtt)
	double precision mn1(6,nnd),mn2(6,nnd)
	integer nov(nrtt)
	common/omega/om,omh
c
c  get lmin and lmax
c
	om=2.d0*pi*f
	omh=om
	call earthmodel_getrearth(rearth)
	elmin=om*rearth/cmax
	elmax=om*rearth/cmin
	if(iprint.gt.1) print *,'<modesf>: lmin = ',elmin,' lmax = ',elmax
c
c  bracket roots for given omega
c
 	call root_bracket(elmin,elmax,nval,'S',iprint)
	if(iprint.gt.1) call root_print
	if(nval.gt.nrtt) then
		write(6,'(a)') 'There will be more roots than dimensioned'
		print *,'Found intervals: ',nval,' Dimension: ',nrtt
		print *,'See dimension nrtt in flspherdim.h'
		stop
	endif
	if(nval.eq.0) then
		nrt=0
		print *,'<modes_spheroidal>: No root was found, return!'
		return
	endif
c
c  pin down bracketed roots
c  numbering of overtones starts from one, but nov gives true
c  overtone number!
c
	call root_pindown(nrt,nov,elroot,m5root,'S',iprint,x1s)
	if(noeig.eq.1) write(6,'(a8,4a10,a13)') 'Overtone','Slowness','Phase','M5'
	do n=1,nrt
		p(n)=(elroot(n)+0.5d0)/(om*rearth)
		if(noeig.eq.1) write(6,'(i8,2f10.4,d13.3)') nov(n),p(n),1.d0/p(n),m5root(n)
	enddo
c
c  skip eigenfunctions if noeig is set
c
	if(noeig.eq.1) goto 1
c
c  establish radial nodes for computation of eigenfunctions
c  and radial integrals like group velocity and Q
c  use smallest starting radius for first node
c
	x1min=x1s(1)
	do n=2,nrt
		x1min=min(x1s(n),x1min)
	enddo
	call nodes_create(f,wlratio,rearth-x1min,0.d0,0.d0,nmin,iprint)
	if(iprint.gt.1) call nodes_print(6,'junk')
c
c  construct eigenfunctions using decomposition of dyad for each mode
c
	call layer_set_stepsize_wlfrac(0.0625d0)
	do n=1,nrt
		elp1=elroot(n)*(elroot(n)+1.d0)
		sqelp1=dsqrt(elp1)
c
c  integration of minors from bottom to top
c  second last flag in sphmot determines what is done
c  save minors from this step in mn1
c
		call sphmot(elroot(n),nlstrt,rstart,ystart,knt,0.d0,iprint,1,0)
		call nodes_getnodeindex(rstart,js)
	 	call savemin(js,sqelp1,mn1)
c
c  integration of minors from top to bottom
c  second last flag in sphmot determines what is done
c  save minors from this step in mn2
c
		call sphmot(elroot(n),nlstrt,rstart,ystart,knt,0.d0,iprint,2,0)
		call savemin(js,sqelp1,mn2)
c
c  construct eigenfunction at nodes from minors mn1 and mn2
c  up to a common sign at each r
c  and later fix sign by applying a smoothness criterion
c-
		call eigfun_create(elp1,nlstrt,js,mn1,mn2,yeig(1,1,n))
		call eigfun_fixsign(nlstrt,js,yeig(1,1,n))
c
c  compute norm and normalize eigenfunctions by
c  sqrt( int ro*(U^2 + l(l+1) V^2) r^2 dr )
c-
		call eigfun_integrals(elp1,nlstrt,js,yeig(1,1,n),'S',
     1		                anorm(n),cg(n),qatt(n),raylquo(n))
		call eigfun_normalize(js,yeig(1,1,n),anorm(n))
c
c  compute Frechet kernels
c
		call eigfun_frechet_sph(elp1,nlstrt,js,yeig(1,1,n),
     1		                  fkro(1,n),fkvp(1,n),fkvs(1,n),fkd(1,n))
	enddo
	call eigfun_printeigval(nrt,nov,p,cg,qatt,anorm,raylquo,m5root)
c
 1	return
	end
