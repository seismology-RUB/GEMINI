!---------------------------------------------------------
!  calculate day of year from year, month, day
!  do not use any more, replaced by timeUtils
!---------------------------------------------------------
 module dayOfYear
	implicit none
!
 contains
!----------------------------------------------------------
!  calculate day of year
!
	integer function getDayOfYear(year,month,day) result(doy)
	integer :: year,month,day,i
	integer, dimension(12) :: dm
	logical :: schalt
!
	dm = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
	schalt = (mod(year,4) == 0)
	if (schalt) dm(2) = 29
	doy = 0
	do i=1,month-1
		doy = doy+dm(i)
	enddo
	doy = doy+day
	end function getDayOfYear
!----------------------------------------------------------
!  calculate month and day from day of year
!
	subroutine monthDayFromDayOfYear(year,doy,month,mday)
	integer :: year,doy,month,mday
	integer, dimension(12), target :: daysum, daysumschalt
	integer, dimension(:), pointer :: daysump
	logical :: schalt
	integer :: jl,ju,jm
	daysum = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
	daysumschalt = (/ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 /)
	schalt = (mod(year,4) == 0)
	if (schalt) then; daysump => daysumschalt; else; daysump => daysum; endif
	jl=0
	ju=13
	do while (ju-jl > 1)
		jm=(ju+jl)/2
		if(doy > daysump(jm)) then; jl=jm; else; ju=jm; endif
	enddo
	month = jl
	mday = doy-daysump(jl)
	end subroutine
!
 end module dayOfYear
