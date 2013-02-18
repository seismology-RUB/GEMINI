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
!--------------------------------------------------------------------
!> \brief  Sort an array
!> \par Description
!>  It works by first finding the smallest element using a linear scan 
!!  and swapping it into the first position in the list, then finding 
!!  the second smallest element by scanning the remaining elements, and so on.
!
 module sortArray
	implicit none
	interface sort
		module procedure realSortArray
		module procedure integerSortArray
		module procedure stringSortArray
	end interface
!
 contains
!---------------------------------------------------------
!> \brief  Sort a real array using selection sort
!> \param  a contains the sorted array
!> \param indx the former indices of the now sorted elements
!
	subroutine realSortArray(a,indx)
	real, dimension(:) :: a
	integer, dimension(:) :: indx
	integer :: i,n,j1,itemp,jm
	integer, dimension(1) :: jmin
	real :: temp
!
	n = size(a)
	j1 = 1
	do i=1,n; indx(i)=i; enddo
	do while (j1 < n)
		jmin=minloc(a(j1:n))
		jm = jmin(1)+j1-1
		temp=a(j1); a(j1)=a(jm); a(jm)=temp
		itemp=indx(j1); indx(j1)=indx(jm); indx(jm)=itemp
		j1=j1+1
	enddo
	end subroutine realSortArray
!---------------------------------------------------------
!> \brief  Sort an integer array using selection sort
!> \param  ia contains the sorted array
!> \param indx the former indices of the now sorted elements
!
	subroutine integerSortArray(ia,indx)
	integer, dimension(:) :: ia
	integer, dimension(:) :: indx
	integer :: i,n,j1,itemp,jm
	integer, dimension(1) :: jmin
	real :: temp
!
	j1 = 1
	n=size(ia)
	do i=1,n; indx(i)=i; enddo
	do while (j1 < n)
		jmin=minloc(ia(j1:n))
		jm = jmin(1)+j1-1
		temp=ia(j1); ia(j1)=ia(jm); ia(jm)=temp
		itemp=indx(j1); indx(j1)=indx(jm); indx(jm)=itemp
		j1=j1+1
	enddo
	end subroutine integerSortArray
!---------------------------------------------------------
!> \brief  Sort a string array using selection sort
!> \param  a contains the sorted array
!> \param indx the former indices of the now sorted elements
!
	subroutine stringSortArray(a,indx)
	character (len=*), dimension(:) :: a
	integer, dimension(:) :: indx
	integer :: i,n,j1,itemp,jm
	character (len=132) :: temp,amin
!
	j1 = 1
	n=size(a)
	do i=1,n; indx(i)=i; enddo
	do while (j1 < n)
		amin = a(j1); jm = j1
		do i = j1,n
			if (llt(a(i),amin)) then
				jm = i; amin = a(i)
			endif
		enddo
		temp = a(j1); a(j1) = a(jm); a(jm) = temp
		itemp = indx(j1); indx(j1) = indx(jm); indx(jm) = itemp
		j1 = j1+1
	enddo
	end subroutine stringSortArray
!------------------------------------------------------------------------
!> \brief Reorder associated arrays to the same sorting order
!
	subroutine reorderIntegerSortArray(ia,indx)
	integer, dimension(:) :: ia,indx
	integer :: n
	integer, dimension(:),allocatable :: iap
!
	n = size(ia)
	allocate(iap(n))
	iap = ia
	ia(1:n) = iap(indx(1:n))
	deallocate(iap)
	end subroutine reorderIntegerSortArray
!-----------------------------------------------------------------------
!> \brief Find the indices of the n smallest elements of an array
!> \param a Input array (remains unchanged)
!
	subroutine findNSmallestRealSortArray(a,indx,nsm)
	real, dimension(:) :: a
	integer, dimension(:) :: indx
	real, dimension(:), allocatable :: acopy
	integer, dimension(:), allocatable :: idxcopy
	integer :: nsm
	integer :: i,n,j1,itemp,jm,jmin
	real :: temp
!
	n = size(a)
	if (nsm > n) nsm = n
	allocate(acopy(n),idxcopy(n))
	acopy = a
	do i = 1,n; idxcopy(i) = i; enddo
	j1 = 1
	do while (j1 <= nsm)
		jmin = minloc(acopy(j1:n),1)
		jm = jmin+j1-1
		temp = acopy(j1); acopy(j1) = acopy(jm); acopy(jm) = temp
		itemp = idxcopy(j1); idxcopy(j1) = idxcopy(jm); idxcopy(jm) = itemp
		j1=j1+1
	enddo
	indx = idxcopy(1:nsm)
	deallocate(acopy,idxcopy)
	end subroutine findNSmallestRealSortArray
!
 end module sortArray
