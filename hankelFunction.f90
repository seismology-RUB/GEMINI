!--------------------------------------------------------------
!  calculate Hankel function to n-th order
!--------------------------------------------------------------
 module hankelFunction
	use mathConstants
	implicit none
	interface hankel
		module procedure hankel_complex
	end interface hankel
 contains
!--------------------------------------------------------------------
	subroutine hankel_complex(n,arg,han)
	integer :: n,i
	real :: arg
	complex, dimension(0:2) :: han
	real :: bessj0,bessj1,bessj,bessy0,bessy1,bessy
!
	han(0)=bessj0(arg)-mc_ci*bessy0(arg)
	han(1)=bessj1(arg)-mc_ci*bessy1(arg)
	do i=2,n
		han(i)=bessj(i,arg)-mc_ci*bessy(i,arg)
	enddo
	end subroutine hankel_complex
 end module hankelFunction
