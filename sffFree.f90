!--------------------------------------------------------
!  module to describe a SFF free block
!--------------------------------------------------------
 module sffFree
	implicit none
	private constructSFFFree, destroySFFFree
	interface new
		module procedure constructSFFFree
	end interface
	interface dealloc
		module procedure destroySFFFree
	end interface dealloc
	type sff_free
		private
		integer :: nline,lenmax
		character (len=80), dimension(:), pointer :: lines
	end type sff_free
!
 contains
!--------------------------------------------------------
!  constructor
!
	subroutine constructSFFFree(this,lu)
	type (sff_free) :: this
	integer :: lu,lindim,ierr
!
	lindim=20
	allocate(this%lines(lindim))
	call sff_RFree(lu,this%nline,this%lines,this%lenmax,lindim,ierr)
	end subroutine constructSFFFree
!---------------------------------------------------------
!  destructor
!
	subroutine destroySFFFree(this)
	type (sff_free) :: this
!
	if(associated(this%lines)) deallocate(this%lines)
	end subroutine destroySFFFree
!----------------------------------------------------------
!  return a pointer to the lines of the free block
!
	function linesSFFFree(this)
	type (sff_free) :: this
	character (len=80), dimension(:), pointer :: linesSFFFree
!
	linesSFFFree => this%lines
	end function linesSFFFree
!-----------------------------------------------------------
!  write a free block to file
!
	subroutine writeSFFFree(this,lu)
	type (sff_free) :: this
	integer :: lu
!
	call sff_WFree(lu,this%nline,this%lines)
	end subroutine writeSFFFree
!-----------------------------------------------------------
 end module sffFree
