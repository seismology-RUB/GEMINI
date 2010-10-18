!--------------------------------------------------------
!  module to describe a SFF source line
!--------------------------------------------------------
 module sffSource
	implicit none
	private constructSFFSource
	interface new
		module procedure constructSFFSource
		module procedure setSFFSource
	end interface
	type sff_source
		private
		character (len=20) :: typ      ! type of source
		character (len=1) :: cs        ! coordinate system
		real :: c1, c2, c3             ! location of source
		character (len=6) :: date      ! date of source
		character (len=10) :: time     ! time of source
	end type sff_source
!
 contains
!-----------------------------------------------------------
!  constructor
!
	subroutine constructSFFSource(this,lu)
	type (sff_source) :: this
	integer :: lu,ierr
!
	call sff_RSource(lu,this%typ,this%cs,this%c1,this%c2,this%c3,this%date,this%time,ierr)
	end subroutine constructSFFSource
!-----------------------------------------------------------
!  second constructor
!
	subroutine setSFFSource(this,typ,cs,c1,c2,c3,date,time)
	type (sff_source) :: this
	character typ*20,cs*1,time*10,date*6
	real :: c1,c2,c3
!
	this = sff_source(typ,cs,c1,c2,c3,date,time)
	end subroutine setSFFSource
!-----------------------------------------------------------
!  write a source line
!
	subroutine writeSFFSource(this,lu)
	type (sff_source) :: this
	integer :: lu
!
	call sff_WSource(lu,this%typ,this%cs,this%c1,this%c2,this%c3,this%date,this%time)
	end subroutine writeSFFSource
!-----------------------------------------------------------
!  print source line
!
	subroutine printSFFSource(this)
	type (sff_source) :: this
	print *,this%typ,this%cs,this%c1,this%c2,this%c3,this%date,this%time
	end subroutine printSFFSource
!---------------------------------------------------------
!  read out source location
!
	subroutine locationSFFSource(this,c1,c2,c3)
	type (sff_source) :: this
	real :: c1,c2,c3
!
	c1 = this%c1; c2 = this%c2; c3 = this%c3
	end subroutine locationSFFSource
!--------------------------------------------------------
!  read out date and time of source
!
	subroutine timeSFFSource(this,date,time)
	type (sff_source) :: this
	character (len=*) :: date,time
!
	date = this%date; time = this%time
	end subroutine timeSFFSource
!--------------------------------------------------------
!  read out source description
!
	subroutine descriptionSFFSource(this,typ,cs)
	type (sff_source) :: this
	character (len=*) :: typ,cs
!
	typ = this%typ; cs = this%cs
	end subroutine descriptionSFFSource
!---------------------------------------------------------
!  get source time in seconds after midnight
!
	real function getSourceTimeSFFSource(this)
	type (sff_source) :: this
	integer :: hh,mm
	real :: sec
!
	read(this%time,'(i2,i2,f6.3)') hh,mm,sec
	getSourceTimeSFFSource = sec+60.*(real(mm)+60.*hh)
	end function getSourceTimeSFFSource
!---------------------------------------------------------
 end module sffSource
