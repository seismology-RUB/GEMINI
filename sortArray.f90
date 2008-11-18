!--------------------------------------------------------------------
!  Sort an array
!  It works by first finding the smallest element using a linear scan 
!  and swapping it into the first position in the list, then finding 
!  the second smallest element by scanning the remaining elements, and so on.
!
 module sortArray
	implicit none
	interface sort
		module procedure realSortArray
	end interface
!
 contains
!---------------------------------------------------------
!  sort a real array using selection sort
!  a contains the sorted array, and indx the former indices 
!  of the now sorted elements
!
	subroutine realSortArray(a,indx)
	real, dimension(:) :: a
	integer, dimension(:) :: indx
	integer :: i,n,j1,itemp,jm
	integer, dimension(1) :: jmin
	real :: temp
!
	j1 = 1
	n=size(a)
	indx = (/ (i, i=1,n) /)
	do, while (j1 < n)
		jmin=minloc(a(j1:n))
		jm = jmin(1)+j1-1
		temp=a(j1); a(j1)=a(jm); a(jm)=temp
		itemp=indx(j1); indx(j1)=indx(jm); indx(jm)=itemp
		j1=j1+1
	enddo
	end subroutine realSortArray
!
 end module sortArray
