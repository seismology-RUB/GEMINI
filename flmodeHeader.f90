!------------------------------------------------------------
!  header part of flmode information
!  for a fixed frequency
!
 module flmodeHeader
	implicit none
	interface new
		module procedure readFlmodeHeader
	end interface
	interface dealloc
		module procedure deallocFlmodeHeader
	end interface
	type flmode_header
		integer :: nnod,jwd,jwu,nrec           ! nrec: number of records for one mode
		double precision :: rearth             ! rearth: earth radius used for calculation
		double precision, dimension(:), pointer :: rnod => null()
		character (len=1) :: typ
	end type flmode_header
!
 contains
!-------------------------------------------------------------
!  read header from direct access file
!  typ gives mode type
!  jrec is current position of file pointer, changed on output
!  znod in meters
!
	subroutine readFlmodeHeader(this,lu,typ,jrec)
	type (flmode_header) :: this
	integer :: lu,jrec
	character (len=*) :: typ
	integer n3,n3rest,k,i
!
	this%typ = typ
	read(lu,rec=jrec) this%rearth,this%nrec
	read(lu,rec=jrec+1) this%nnod,this%jwd,this%jwu
	n3=this%nnod/3
	n3rest=mod(this%nnod,3)
	allocate(this%rnod(this%nnod))
	do k=1,n3
		read(lu,rec=jrec+1+k) (this%rnod((k-1)*3+i),i=1,3)
	enddo
	jrec=jrec+1+n3+1
	if(n3rest.gt.0) then
		read(lu,rec=jrec) (this%rnod(n3*3+i),i=1,n3rest)
		jrec=jrec+1
	endif
!
	end subroutine readFlmodeHeader
!------------------------------------------------------------
!  deallocate header
!
	subroutine deallocFlmodeHeader(this)
	type (flmode_header) :: this
	if(associated(this%rnod)) deallocate(this%rnod)
	end subroutine deallocFlmodeHeader
!------------------------------------------------------------
!  determine index of node below given radius
!  returns zero if r is below lowest node
!
	integer function getLowerNodeIndexFlmodeHeader(this,r)
	type (flmode_header) :: this
	double precision :: r
	integer :: ids
!
	call dlocate(r,this%nnod,this%rnod,ids)
	getLowerNodeIndexFlmodeHeader = ids
	end function getLowerNodeIndexFlmodeHeader
!------------------------------------------------------------
!  determine index of node above given radius
!  returns nnod+1 if r is above topmost node
!
	integer function getUpperNodeIndexFlmodeHeader(this,r) result(ids)
	type (flmode_header) :: this
	double precision :: r
	integer :: j
!
	call dlocate(r,this%nnod,this%rnod,j)
	ids = j+1
	end function getUpperNodeIndexFlmodeHeader
!-----------------------------------------------------------
!  get number of records used by one set of eigenfunctions
!
	integer function getNrecFlmodeHeader(this)
	type (flmode_header) :: this
	getNrecFlmodeHeader = this%nrec
	end function getNrecFlmodeHeader
!---------------------------------------------------------------------------------
!  return a pointer to the radial nodes
!
	function getNodesFlmodeHeader(this)
	type (flmode_header) :: this
	double precision, dimension(:), pointer :: getNodesFlmodeHeader
!
	getNodesFlmodeHeader => this%rnod
	end function getNodesFlmodeHeader
!---------------------------------------------------------------------------------
!  get mode type
!
	function getModeTypeFlmodeHeader(this)
	type (flmode_header) :: this
	character (len=1) :: getModeTypeFlmodeHeader
!
	getModeTypeFlmodeHeader = this%typ
	end function getModeTypeFlmodeHeader
!---------------------------------------------------------------------------------
!  get earth radius
!
	function getEarthRadiusFlmodeHeader(this)
	type (flmode_header) :: this
	double precision :: getEarthRadiusFlmodeHeader
!
	getEarthRadiusFlmodeHeader = this%rearth
	end function getEarthRadiusFlmodeHeader
!--------------------------------------------------------------------
 end module flmodeHeader
