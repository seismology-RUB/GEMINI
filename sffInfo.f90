!--------------------------------------------------------
!  module to describe a SFF info line
!--------------------------------------------------------
 module sffInfo
	implicit none
	private constructSFFInfo
	interface new
		module procedure constructSFFInfo
		module procedure setSFFInfo
	end interface
	type sff_info
		private
		character (len=1) :: cs        ! coordinate system
		real :: c1, c2, c3             ! location of source
		integer :: nstack              ! number of stacks
	end type sff_info
!
 contains
!-----------------------------------------------------------
!  constructor
!
	subroutine constructSFFInfo(this,lu)
	type (sff_info) :: this
	integer :: lu,ierr
!
	call sff_RInfo(lu,this%cs,this%c1,this%c2,this%c3,this%nstack,ierr)
	end subroutine constructSFFInfo
!---------------------------------------------------------
!  second constructor
!
	subroutine setSFFInfo(this,cs,c1,c2,c3,nstack)
	type (sff_info) :: this
	character (len=1) :: cs
	real :: c1,c2,c3
	integer :: nstack
!
	this = sff_info(cs,c1,c2,c3,nstack)
	end subroutine setSFFInfo
!-----------------------------------------------------------
!  return receiver location
!
	subroutine locationSFFInfo(this,x,y,z)
	type (sff_info) :: this
	real :: x,y,z
	x = this%c1; y = this%c2; z = this%c3
	end subroutine locationSFFInfo
!-----------------------------------------------------------
!  print info line
!
	subroutine printSFFInfo(this)
	type (sff_info) :: this
	print *,this%cs,this%c1,this%c2,this%c3,this%nstack
	end subroutine printSFFInfo
!---------------------------------------------------------
!  write info line to file
!
	subroutine writeSFFInfo(this,lu)
	type (sff_info) :: this
	integer :: lu
!
	call sff_WInfo(lu,this%cs,this%c1,this%c2,this%c3,this%nstack)
	end subroutine writeSFFInfo
!---------------------------------------------------------
 end module sffInfo
