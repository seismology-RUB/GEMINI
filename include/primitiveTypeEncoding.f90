!----------------------------------------------------------------------
!  conventions for identifying pritive data types
!----------------------------------------------------------------------
 module primitiveTypeEncoding
	implicit none
	integer, parameter :: T_INTEGER = 1
	integer, parameter :: T_REAL = 2
	integer, parameter :: T_DOUBLE = 3
	integer, parameter :: T_COMPLEX = 4
	integer, parameter :: T_DOUBLE_COMPLEX = 5
	integer, parameter :: T_CHAR = 6
	integer, parameter :: T_FLEXIBLE = 7
	integer, parameter :: T_LOGICAL = 8
end module

