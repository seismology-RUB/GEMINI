!-------------------------------------------------------------
!  color scale module
!-------------------------------------------------------------
 module pgColorScale
	private :: checkColorCapability, createColorScale, applyColorScale
	interface new; module procedure createColorScale; end interface
	interface apply; module procedure applyColorScale; end interface
	interface dealloc; module procedure deallocColorScale; end interface
	type pg_colscale
		private
		integer :: nci,cimin,cimax     ! number of colors, min index, max index
		real, dimension(:), pointer :: red,green,blue
	end type pg_colscale
!
	contains
!-------------------------------------------------------------
!  create a new color scale object
!
	subroutine createColorScale(this,name)
	type (pg_colscale) :: this
	character (len=*) :: name
!
!  branch to different scales
!
	select case (name)
	case ('rwb')
		call redWhiteBlue(this,.false.)
	case ('bwr')
		call redWhiteBlue(this,.true.)
	case ('rygb')
		call redYellowGreenBlue(this,.false.)
	case ('bgyr')
		call redYellowGreenBlue(this,.true.)
!	case ('bw')
!		call blackWhite(this)
!	case ('wb')
!		call whiteBlack(this)
!	case ('red')
!		call red(this)
!	case ('blue')
!		call blue(this)
	case default
		print *,'ColorScale unknown, use red-white-blue'
		call redWhiteBlue(this,.false.)
	end select
	end subroutine createColorScale
!------------------------------------------------------------
!  apply color scale
!
	subroutine applyColorScale(this)
	type (pg_colscale) :: this
	integer :: res
!
	call checkColorCapability(this)
	res = pgscir(this%cimin, this%cimax)
	do i=1,this%nci
		res = pgscr(this%cimin+i-1,this%red(i),this%green(i),this%blue(i))
	enddo
	end subroutine applyColorScale
!-------------------------------------------------------------
!  deallocate color scale
!
	subroutine deallocColorScale(this)
	type (pg_colscale) :: this
!
	if(associated(this%red)) deallocate(this%red)
	if(associated(this%green)) deallocate(this%green)
	if(associated(this%blue)) deallocate(this%blue)
	end subroutine deallocColorScale
!-------------------------------------------------------------
!  check color capability of device
!
	subroutine checkColorCapability(this)
	type (pg_colscale) :: this
	integer :: ci1,ci2
!
	call pgqcol(ci1,ci2)
	if(ci2 == 1) then
		print *,'Device does not have color capability'
		stop
	endif
	if(this%cimax > ci2) then
		print *,'Too few colors available. Use black and white'
		print *,this%cimax,ci2
		stop
	endif
	end subroutine checkColorCapability
!--------------------------------------------------------------
!  red-white-blue
!
	subroutine redWhiteBlue(this,invflag)
	type (pg_colscale) :: this
	logical :: invflag
	integer :: n,i,nval,ci2
!
	ci2 = 99
	this%cimin=16
!
	n = ((ci2 - 15)/2)*2-1            ! use odd number of indices
	this%nci = n
	this%cimax = this%cimin + n - 1
	allocate(this%red(n), this%green(n), this%blue(n))
!
	nval = n-1
	if(.not.invflag) then
		do i=1,(n+1)/2
			this%red(i) = 1; this%green(i) = (2*i-2)/real(nval); this%blue(i) = this%green(i)
			j=n-i+1
			this%red(j) =  this%green(i); this%green(j) = this%green(i); this%blue(j) = this%red(i) 
		enddo
	else
		do i=1,(n+1)/2
			this%blue(i) = 1; this%green(i) = (2*i-2)/real(nval); this%red(i) = this%green(i)
			j=n-i+1
			this%blue(j) =  this%green(i); this%green(j) = this%green(i); this%red(j) = this%blue(i) 
		enddo
	endif
	end subroutine redWhiteBlue
!--------------------------------------------------------------
!  red-yellow-green-blue
!
	subroutine redYellowGreenBlue(this,invflag)
	type (pg_colscale) :: this
	logical :: invflag
	integer :: n,i,nval,nm,ci2
!
	ci2=99
	this%cimin=16
!
	n = ((ci2 - 15)/2)*2           ! use even number of indices
	this%nci = n
	this%cimax = this%cimin + n - 1
	allocate(this%red(n), this%green(n), this%blue(n))
!
	nval = n-1
	nm=n/2
	if(.not.invflag) then
		do i=1,nm
			this%red(i) = 1.; this%green(i) = (i-1)/real(nm-1); this%blue(i) = 0.
		enddo
		do i=nm+1,n
			this%red(i) = 0.; this%blue(i)=(i-nm-1)/real(nm-1); this%green(i) = 1.-this%blue(i)
		enddo
	else
		do i=1,nm
			this%red(i) = 0.; this%green(i) = (i-1)/real(nm-1); this%blue(i) = 1.-this%green(i)
		enddo
		do i=nm+1,n
			this%red(i) = 1.; this%green(i)=1.-(i-nm-1)/real(nm-1); this%blue(i) = 0.
		enddo
	endif
	end subroutine redYellowGreenBlue
!----------------------------------------------------------------
!  evaluate color index from function value and min/max values
!  linear mapping, full-scale if saturated
!
	integer function getIndexColorScale(this,f,fmin,fmax)
	type (pg_colscale) :: this
	real :: f,fmin,fmax
	integer :: ci
	ci = this%cimin+(f-fmin)/(fmax-fmin)*(this%cimax-this%cimin)
	if (ci < this%cimin) ci = this%cimin
	if (ci > this%cimax) ci = this%cimax
	getIndexColorScale = ci
	end function getIndexColorScale
!---------------------------------------------------------------
 end module pgColorScale
