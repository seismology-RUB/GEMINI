!--------------------------------------------------------------------------
!	Copyright 2013 Wolfgang Friederich
!
!	This file is part of Gemini II.
!
!	Gemini II is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 2 of the License, or
!	any later version.
!
!	Gemini II is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
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
