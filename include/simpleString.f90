!----------------------------------------------------------------------------
!  module for a string
!----------------------------------------------------------------------------
 module simpleString
	implicit none
	interface new; module procedure createSimpleString; end interface
	interface dealloc; module procedure deallocSimpleString; end interface
	interface operator (.length.); module procedure lengthSimpleString; end interface
	type simple_string
		private
		character (len=1), dimension(:), pointer :: zeichen
	end type simple_string
!
 contains
!----------------------------------------------------------------------------
!  create a string object from a given character sequence
!  strip trailing blanks
!
	subroutine createSimpleString(this,string)
	type (simple_string) :: this
	character (len=*) :: string
	integer :: i
!
!  copy string onto zeichen
!
	allocate(this%zeichen(len_trim(string)))
	do i=1,len_trim(string)
		this%zeichen(i) = string(i:i)
	enddo
	end subroutine createSimpleString
!----------------------------------------------------------------------------
!  deallocate
!
	subroutine deallocSimpleString(this)
	type (simple_string) :: this
	if (associated(this%zeichen)) deallocate(this%zeichen)
	end subroutine deallocSimpleString
!----------------------------------------------------------------------------
!  get length of string
!
	integer function lengthSimpleString(this)
	type (simple_string), intent(in) :: this
	lengthSimpleString = size(this%zeichen)
	end function lengthSimpleString
!----------------------------------------------------------------------------
!  convert string back to character sequence
!
	subroutine convertToCharSimpleString(this,string)
	type (simple_string) :: this
	character (len=*) :: string
	integer :: n,i
	n=min(size(this%zeichen),len(string))
	do i=1,n
		string(i:i) = this%zeichen(i)
	enddo
	do i=n+1,len(string)
		string(i:i) = ' '
	enddo
	end subroutine convertToCharSimpleString
!----------------------------------------------------------------------------	
 end module simpleString
