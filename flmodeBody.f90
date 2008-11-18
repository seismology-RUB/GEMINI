!------------------------------------------------------------
!  body part of flmode information
!
!  subroutine readFlmodeBody(this,lu,header,jrec)
!  subroutine deallocFlmodeBody(this)
!  subroutine evalSphEigfunFlmodeBody(this,r,u,up,v,vp)
!  subroutine evalTorEigfunFlmodeBody(this,r,w,wp)
!  subroutine excitationFlmodeBody(this,rm,rs,psi)
!  function   displacementFlmodeBody(this,csys,comp,dis,phi,azi,psi)
!  subroutine potentialFlmodeBody(this,dis,phi,psi,pot,iswitch)
!  function getEigfunFlmodeBody(this,k)
!  function getDispersionFlmodeBody(this)
!
 module flmodeBody
	use flmodeHeader
	use mathConstants
	use hankelFunction
	implicit none
	interface new
		module procedure readFlmodeBody
	end interface
	interface dealloc
		module procedure deallocFlmodeBody
	end interface
	interface evalEigfunFlmodeBody
		module procedure evalSphEigfunFlmodeBody
		module procedure evalTorEigfunFlmodeBody
	end interface
	type flmode_body
		private
		integer :: nov
		real :: f,p,cg,qatt     !  p in s/km, cg in km/s, f in Hz, qatt dimensionless
		real, dimension(:,:), pointer :: yeig => null()
		type (flmode_header), pointer :: header => null() ! pointer to header information of node
	end type flmode_body
!
 contains
!-------------------------------------------------------------
!  read one mode from direct access file
!  jrec is current position of file pointer, changed on output
!  header is flmodeHeader information belonging to mode
!
	subroutine readFlmodeBody(this,lu,header,jrec)
	type (flmode_body) :: this
	type (flmode_header), target :: header
	integer :: lu, jrec, ios
	integer :: j,n2,n2rest,nnod,i
!
	this%header => header
	nnod = header%nnod
	read(lu,rec=jrec) this%f,this%nov,this%p,this%cg,this%qatt
	if(header%typ == 'S') then
		allocate(this%yeig(6,nnod))
		do j=1,nnod
			read(lu,rec=jrec+j) (this%yeig(i,j),i=1,6)
		enddo
		jrec=jrec+nnod+1
	else if(header%typ == 'T') then
		n2 = nnod/2
		n2rest = mod(nnod,2)
		allocate(this%yeig(3,nnod))
		do j=1,2*n2,2
			jrec=jrec+1
			read(lu,rec = jrec) (this%yeig(i,j),i=1,3),(this%yeig(i,j+1),i=1,3)
		enddo
		jrec=jrec+1
		if(n2rest.gt.0) then
			read(lu,rec = jrec, iostat = ios) (this%yeig(i,nnod),i=1,3)
!
!  The following is a workaround for a possible compiler bug where in some cases the last record in the
!  file can no more be read. But tests show that this record has been written!
!  Since there is a double node at the surface anyway because in nodes_create the receiver is put
!  at the surface we have the correct eigenfunction value already one node before and therefore use it
!  when the error occurs.
!
			if (ios > 0) then
				print *,'Error reading modes'
				this%yeig(1:3,nnod) = this%yeig(1:3,nnod-1)
!				print *,'Error reading modes: ',ios,jrec,nnod,n2,n2rest,this%nov,this%f
!				print *,(this%yeig(i,nnod),i=1,3)
!				print *,'Ignore error'
			endif
			jrec=jrec+1
		endif
	else
		print *,'Unknown mode type: ',header%typ
		stop
	endif
	end subroutine readFlmodeBody
!--------------------------------------------------------------
!  deallocate
!
	subroutine deallocFlmodeBody(this)
	type (flmode_body) :: this
!
	if(associated(this%yeig)) deallocate(this%yeig)
	if(associated(this%header)) nullify(this%header)
	end subroutine deallocFlmodeBody
!-------------------------------------------------------------------------------------------------
!  evaluate spheroidal eigenfunction at given depth using spline interpolation between neighbouring nodes
!
	subroutine evalSphEigfunFlmodeBody(this,r,u,up,v,vp)
	use splineInterpol
	type (flmode_body) :: this
	double precision :: r,rb,ra
	real :: u,up,v,vp
	integer :: id
!
	id = getLowerNodeIndexFlmodeHeader(this%header,r)
	if (id == 0) then; u=0.; up=0.; v=0.; vp=0.; endif
	rb = this%header%rnod(id)
	ra = this%header%rnod(id+1)
	call splineInterpol2P(r,rb,ra,this%yeig(1,id),this%yeig(1,id+1),this%yeig(5,id),this%yeig(5,id+1),u,up)
	call splineInterpol2P(r,rb,ra,this%yeig(3,id),this%yeig(3,id+1),this%yeig(6,id),this%yeig(6,id+1),v,vp)
	end subroutine evalSphEigfunFlmodeBody
!----------------------------------------------------------------------------------------------
!  evaluate toroidal eigenfunction at given depth using spline interpolation between neighbouring nodes
!  return zero if r is out of node range
!
	subroutine evalTorEigfunFlmodeBody(this,r,w,wp)
	use splineInterpol
	type (flmode_body) :: this
	double precision :: r,rb,ra
	real :: w,wp
	integer :: id
!
	id = getLowerNodeIndexFlmodeHeader(this%header,r)
	if (id == 0) then; w=0.; wp=0.; endif
	rb = this%header%rnod(id)
	ra = this%header%rnod(id+1)
	call splineInterpol2P(r,rb,ra,this%yeig(1,id),this%yeig(1,id+1),this%yeig(3,id),this%yeig(3,id+1),w,wp)
	end subroutine evalTorEigfunFlmodeBody
!-----------------------------------------------------------------------------
!  calculate excitation coefficients for given depth and moment tensor
!
	subroutine excitationFlmodeBody(this,rm,rs,psi)
	
	type (flmode_body) :: this
	real, dimension(:) :: rm
	double precision :: rs
	complex, dimension(0:) :: psi
	real :: wu2,vl,vl2,usf,upsf,vsf,vpsf,wsf,wpsf,fs,xs,zs
!
	wu2=sqrt(2.)
	vl = 2.*pi * this%f * this%header%rearth * this%p
	vl2=vl*vl
	if(this%header%typ == 'S') then
		call evalEigfunFlmodeBody(this,rs,usf,upsf,vsf,vpsf)
		fs=(2.*usf-vl2*vsf)/rs
		xs=vpsf+(usf-vsf)/rs
		psi(0)=rm(1)*upsf+(rm(2)+rm(3))*.5*fs
		psi(1)=(-rm(4)+mc_ci*rm(5))*xs/wu2
		psi(2)=(.5*(rm(2)-rm(3))-mc_ci*rm(6))*.5*vsf/rs
	else if(this%header%typ == 'T') then
		call evalEigfunFlmodeBody(this,rs,wsf,wpsf)
		zs=wpsf-wsf/rs
		psi(0)=0.
		psi(1)=(-rm(4)+mc_ci*rm(5))*mc_ci*zs/wu2
		psi(2)=(.5*(rm(2)-rm(3))-mc_ci*rm(6))*.5*mc_ci*wsf/rs
	else
		print *,'excitationFlmodeBody: unknown mode type'
		stop
	endif
	end subroutine excitationFlmodeBody
!-------------------------------------------------------------------------------------
!  add mode contribution to displacement for given receiver at surface or ocean bottom
!
	complex function displacementFlmodeBody(this,csys,comp,dis,phi,azi,psi)
	type (flmode_body) :: this
	complex, dimension(0:) :: psi
	real :: dis,phi
	character (len=*) :: comp,csys
	real :: ur,vr,rr,wr,azi,caz,saz
	complex :: pot,dt_pot,df_pot,zspl,zspt,zspn,zspe
!
!  spheroidal displacement
!
	if(this%header%typ == 'S') then
		ur = this%yeig(1,this%header%nnod)
		vr = this%yeig(3,this%header%nnod)
		rr = this%yeig(2,this%header%nnod)
		if(this%header%jwd .gt. 0) then
			ur = this%yeig(1,this%header%jwd)
			vr = this%yeig(3,this%header%jwd)
			rr = this%yeig(2,this%header%jwd)
		endif
		if(comp.eq.'Z' .or. comp.eq.'R') then
			call potentialFlmodeBody(this,dis,phi,psi,pot,1)
			displacementFlmodeBody = ur*pot
			return
		else if(comp.eq.'H') then
			call potentialFlmodeBody(this,dis,phi,psi,pot,1)
			displacementFlmodeBody = rr*pot
			return
		else if(comp.eq.'L') then
			call potentialFlmodeBody(this,dis,phi,psi,dt_pot,2)
			displacementFlmodeBody = vr*dt_pot
			return
		else if(comp.eq.'T') then
			call potentialFlmodeBody(this,dis,phi,psi,df_pot,3)
			displacementFlmodeBody = vr*df_pot
			return
		else if(comp.eq.'N' .or.comp.eq.'E') then
			call potentialFlmodeBody(this,dis,phi,psi,dt_pot,2)
			call potentialFlmodeBody(this,dis,phi,psi,df_pot,3)
			zspl = vr*dt_pot
			zspt = vr*df_pot
		endif
!
!  toroidal displacement
!
	else if(this%header%typ == 'T') then
		wr = this%yeig(1,this%header%nnod)
		if(this%header%jwd .gt. 0) then
			wr = this%yeig(1,this%header%jwd)
		endif
		if(comp.eq.'Z' .or. comp.eq.'R' .or. comp.eq.'H') then
			displacementFlmodeBody = (0., 0.)
			return
		else if(comp.eq.'L') then
			call potentialFlmodeBody(this,dis,phi,psi,df_pot,3)
			displacementFlmodeBody = -wr*df_pot
			return
		else if(comp.eq.'T') then
			call potentialFlmodeBody(this,dis,phi,psi,dt_pot,2)
			displacementFlmodeBody = wr*dt_pot
			return
		else if(comp.eq.'N' .or.comp.eq.'E') then
			call potentialFlmodeBody(this,dis,phi,psi,dt_pot,2)
			call potentialFlmodeBody(this,dis,phi,psi,df_pot,3)
			zspl = -wr*df_pot
			zspt = +wr*dt_pot
		endif
	endif
!
!  rotation for N or E component
!
	if(comp.eq.'N' .or.comp.eq.'E') then
		caz=cos(azi)
		saz=sin(azi)
		if(csys.eq.'S') then
			zspn=-zspl*caz+zspt*saz
			zspe=+zspl*saz+zspt*caz
		else
			zspn=+zspl*caz-zspt*saz
			zspe=+zspl*saz+zspt*caz
		endif
		if(comp.eq.'N') displacementFlmodeBody = zspn
		if(comp.eq.'E') displacementFlmodeBody = zspe
	endif
!
	end function displacementFlmodeBody
!------------------------------------------------------------------------------
!  calculate potential for mode
!
!  iswitch:  1 = potential itself
!            2 = d/dtheta of potential
!            3 = 1/sin(theta) d/dphi of potential
!-----------------------------------------------------------------------
	subroutine potentialFlmodeBody(this,dis,phi,psi,pot,iswitch)	
	type (flmode_body) :: this
	complex, dimension(0:) :: psi
	complex :: pot
	integer iswitch,k
	real :: dis,phi,wu2,cph,eps,vl,vl2,arg,arg2,expo
	real :: r1psph,r2psph,cosph,sinph,cos2ph,sin2ph,dfr1psph,dfr2psph
	real bessj, bessy, bessj0, bessj1, bessy0, bessy1
	double precision :: delta,sdelta,deltainv
	complex :: han(0:2),feldfac0,feldfac,dhan(0:2)
!
	wu2=sqrt(2.)
!
	cph = 1./(this%header%rearth*this%p)
	eps=(pi*this%f)/(this%qatt*this%cg/this%header%rearth)
	feldfac0=-mc_ci/(4.*cph*this%cg/this%header%rearth)
	vl=2.*pi*this%f/cph
	vl2=vl*vl
!
!  Compute reference wavefield of exciting mode at receivers
!  
	delta=dble(dis)/this%header%rearth
	deltainv=1.d0/delta
	sdelta=dsin(delta)
	feldfac=feldfac0*dsqrt(delta/sdelta)
	arg=dble(vl)*delta
	expo=dble(eps)*delta
!
!  Hankel functions
!
	if(arg.lt.6.*pi) then
		call hankel(2,arg,han)
	else 
		do k=0,2
			arg2=arg-.5*pi*(k+.5)
			han(k)=sqrt(2./(pi*arg))*cmplx(cos(arg2),-sin(arg2))*exp(-expo)
		enddo
	endif
!
	dhan(0)=-vl*han(1)
	dhan(1)=-deltainv*han(1)+vl*han(0)
	dhan(2)=-2.*deltainv*han(2)+vl*han(1)
!
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
	end subroutine potentialFlmodeBody
!---------------------------------------------------------------------------------
!  return a pointer to a certain eigenfunction
!
	function getEigfunFlmodeBody(this,k)
	type (flmode_body) :: this
	integer :: k
	real, dimension(:), pointer :: getEigfunFlmodeBody
!
	if(k < 1 .or. k > size(this%yeig,1)) then
		print *,'<getEigfunFlmodeBody>: invalid eigenfunction index: ',k
		stop
	endif 
	getEigfunFlmodeBody => this%yeig(k,:)
	end function getEigfunFlmodeBody
!---------------------------------------------------------------------------------
!  return a pointer to all eigenfunctions
!
	function getEigenfunctionsFlmodeBody(this)
	type (flmode_body) :: this
	real, dimension(:,:), pointer :: getEigenfunctionsFlmodeBody
!
	getEigenfunctionsFlmodeBody => this%yeig
	end function getEigenfunctionsFlmodeBody
!----------------------------------------------------------------------------------
!  get spheroidal eigenfunction values at given node
!
	subroutine getEigfunNodeSpheroidalFlmodeBody(this,j,u,up,v,vp)
	type (flmode_body) :: this
	integer :: j
	real :: u,up,v,vp
!
	u = this%yeig(1,j); up = this%yeig(5,j); v = this%yeig(3,j); vp = this%yeig(6,j)
	end subroutine getEigfunNodeSpheroidalFlmodeBody
!----------------------------------------------------------------------------------
!  return f,p,c,cg,q
!
	function getDispersionFlmodeBody(this)
	type (flmode_body) :: this
	real, dimension(5) :: getDispersionFlmodeBody
	getDispersionFlmodeBody = (/ this%f, this%p, 1./this%p, this%cg, this%qatt /)
	end function getDispersionFlmodeBody
!----------------------------------------------------------------------------------
!  return harmonic degree
!
	function getHarmonicDegreeFlmodeBody(this)
	type (flmode_body) :: this
	real :: getHarmonicDegreeFlmodeBody
	getHarmonicDegreeFlmodeBody = 2.*pi * this%f * this%header%rearth * this%p
	end function getHarmonicDegreeFlmodeBody
!----------------------------------------------------------------------------------
!  return angular frequency
!
	function getOmegaFlmodeBody(this)
	type (flmode_body) :: this
	real :: getOmegaFlmodeBody
	getOmegaFlmodeBody = 2.*pi * this%f
	end function getOmegaFlmodeBody
!----------------------------------------------------------------------------------
!  return sqrt(4cUv) with angular velocities c and U
!  4cUv = 4/p*U/R^2*v = 4/p*U/R**2*om*R*p = 4*U*om/R
!
	real function sqr4cuvFlmodeBody(this)
	type (flmode_body) :: this
	sqr4cuvFlmodeBody = sqrt(4.*this%cg*2.*pi*this%f/this%header%rearth)
	end function sqr4cuvFlmodeBody
!----------------------------------------------------------------------------------
!  add something to the slowness
!
	subroutine addSlownessFlmodeBody(this,delp)
	type (flmode_body) :: this
	real :: delp
	this%p = this%p + delp
	end subroutine addSlownessFlmodeBody
!
 end module flmodeBody
