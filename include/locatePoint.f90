!--------------------------------------------------------------
!  module with locate functions
!--------------------------------------------------------------------
!  locate index of point xx with x(i) < xx <= x(i+1)
!  returns value in j
!  returns 0 or n if xx is off-scale
!-------------------------------------------------------
 module locatePoint
	implicit none
	interface locate
		module procedure locate_real
		module procedure locate_double_precision
		module procedure locate_int
	end interface locate
 contains
!--------------------------------------------------------------
!  real array
!
	integer function locate_real(xx,n,x)
	real, dimension(:) :: x
	real :: xx
	integer :: n,jl,ju,jm
!
	jl=0
	ju=n+1
	do while (ju-jl > 1)
		jm=(ju+jl)/2
		if((x(n).gt.x(1)).eqv.(xx.gt.x(jm))) then
			jl=jm
		else
			ju=jm
		endif
	enddo
	locate_real = jl
	end function locate_real
!--------------------------------------------------------------------
!  double precision array
!
	integer function locate_double_precision(xx,n,x)
	double precision, dimension(:) :: x
	double precision :: xx
	integer :: n,jl,ju,jm
!
	jl=0
	ju=n+1
	do while (ju-jl > 1)
		jm=(ju+jl)/2
		if((x(n).gt.x(1)).eqv.(xx.gt.x(jm))) then
			jl=jm
		else
			ju=jm
		endif
	enddo
	locate_double_precision = jl
	end function locate_double_precision
!--------------------------------------------------------------
!  integer array
!
	integer function locate_int(xx,n,x)
	integer, dimension(:) :: x
	integer :: xx
	integer :: n,jl,ju,jm
!
	jl=0
	ju=n+1
	do while (ju-jl > 1)
		jm=(ju+jl)/2
		if((x(n).gt.x(1)).eqv.(xx.gt.x(jm))) then
			jl=jm
		else
			ju=jm
		endif
	enddo
	locate_int = jl
	end function locate_int
!
 end module locatePoint
