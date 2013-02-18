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
!  module for an arbitrary long string
!  potentially separated into words by a separating character like one or
!  more blanks or dots or slashes and so on
!----------------------------------------------------------------------------
 module longString
	implicit none
	interface new
		module procedure createLongString
	end interface
	interface dealloc
		module procedure deallocLongString
	end interface
	interface operator(.nword.); module procedure numberOfWordsLongString; end interface
	interface operator(.word.);  module procedure wordLongString; end interface
	type long_string
		private
		character (len=1), dimension(:), pointer :: zeichen
		character (len=1) :: separator
		integer, dimension(:), pointer :: wordstart
		integer, dimension(:), pointer :: wordstop
		integer :: nword
	end type
!
 contains
!----------------------------------------------------------------------------
!  create a long string from a given string
!  strip trailing blanks, words separated by one or more separator characters
!  which default to blanks
!
	subroutine createLongString(this,string,separator)
	type (long_string) :: this
	character (len=*) :: string
	integer :: i,ia,k,ka,ke
	character (len=1), optional :: separator
!
	if (present(separator)) then; this%separator = separator; else; this%separator = ' '; endif
!
!  copy string onto zeichen
!
	allocate(this%zeichen(len_trim(string)))
	do i=1,len_trim(string)
		this%zeichen(i) = string(i:i)
	enddo
!
!  find index of first character which is not a separator
!
	ia = 1; do while(this%zeichen(ia) == this%separator); ia = ia+1; enddo
!
!  find number of words by looking for transitions from separating character to non-separating character
!
	k = 1
	do i=ia,size(this%zeichen)-1
		if(this%zeichen(i) == this%separator .and. this%zeichen(i+1) /= this%separator) k=k+1
	enddo
	this%nword = k
!
!  set starting locations of words at transitions from separator to non-separator
!  set ending   locations of words at transitions from non-separartor to separator
!
	allocate(this%wordstart(k))
	allocate(this%wordstop(k))
	this%wordstart(1) = ia
	this%wordstop(k) = size(this%zeichen)
	ka=1
	ke=0
	do i=ia,size(this%zeichen)-1
		if(this%zeichen(i) == this%separator .and. this%zeichen(i+1) /= this%separator) then; ka=ka+1; this%wordstart(ka) = i+1; endif
		if(this%zeichen(i) /= this%separator .and. this%zeichen(i+1) == this%separator) then; ke=ke+1; this%wordstop(ke) = i; endif
	enddo
	end subroutine createLongString
!----------------------------------------------------------------------------
!  deallocate longString
!
	subroutine deallocLongString(this)
	type (long_string) :: this
	if(associated(this%zeichen)) deallocate(this%zeichen)
	if(associated(this%wordstart)) deallocate(this%wordstart)
	if(associated(this%wordstop)) deallocate(this%wordstop)
	end subroutine deallocLongString
!----------------------------------------------------------------------------
!  find out how many blank-separated words are contained in a long string
!
	integer function numberOfWordsLongString(this)
	type (long_string), intent(in) :: this
!
	numberOfWordsLongString = this%nword
	end function numberOfWordsLongString
!----------------------------------------------------------------------------
!  get k-th word from a long string
!  return empty string if there is none
!
	character (len=232) function wordLongString(this,k)
	type (long_string), intent(in) :: this
	integer, intent(in) :: k
	integer :: i,j
	character (len=232) :: string
!
	j=0
	do i=this%wordstart(k),this%wordstop(k)
		j=j+1
		string(j:j) = this%zeichen(i)
	enddo
	wordLongString = string(1:j)
!
	end function wordLongString
!-----------------------------------------------------------------------
!  search words of longString for equality with a given string
!  returns the number of words that are equal to string
!
	integer function searchWordsEqualLongString(this,string)
	type (long_string), intent(in) :: this
	character (len=*) :: string
	integer :: cnt,k
!
	cnt = 0
	do k=1,this%nword
		if (trim(wordLongString(this,k)) == trim(string)) cnt = cnt+1
	enddo
	searchWordsEqualLongString = cnt
	end function searchWordsEqualLongString
!
 end module longString
!
!  TODO
!
!  - appendString
!  - convertToString
!  - printLongString
