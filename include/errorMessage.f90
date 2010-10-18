!----------------------------------------------------------------------
!> \brief  Module to handle errors occuring in functions or subroutines
!> \par Description
!> Functions producing exceptions may return an error_message object to
!! specify the name of the function where the exception happened, to
!! provide a description of the error and to rate its severity
!! ranging from OK (no error), warning to error.
!<--------------------------------------------------------------------- 
 module errorMessage
       use iso_fortran_env
	use realloc
	implicit none
	interface new
		module procedure createErrorMessage
		module procedure createOKErrorMessage
	end interface
	interface dealloc; module procedure deallocErrorMessage; end interface
	interface print; module procedure printErrorMessage; end interface
	interface operator (.level.); module procedure getLevelErrorMessage; end interface
	type error_message
		private
		integer :: level                  !< error level: 0 = success, 1 = warning, 2 = error
		character (len=132) :: message    !< error message
		character (len=132) :: fctname    !< function name where error happened
		character (len=132), dimension(:), pointer :: trace => null()  !< function names through which error was propagated
	end type
!
 contains
!-----------------------------------------------------------------
!> \brief Create an error message
!
	subroutine createErrorMessage(this,level,message,fctname)
	type (error_message) :: this
	integer :: level
	character (len=*) :: message,fctname
	this = error_message(level,message,fctname,null())
	end subroutine createErrorMessage
!------------------------------------------------------------------
!> \brief Create a default OK message
!
	subroutine createOKErrorMessage(this,fctname)
	type (error_message) :: this
	character (len=*) :: fctname
	this = error_message(0,'success',fctname,null())
	end subroutine createOKErrorMessage
!-----------------------------------------------------------------
!> \brief Deallocate error message
!
	subroutine deallocErrorMessage(this)
	type (error_message) :: this
	if (associated(this%trace)) deallocate(this%trace)
	end subroutine deallocErrorMessage
!------------------------------------------------------------------
!> \brief Print error message
!
	subroutine printErrorMessage(this)
	type (error_message), intent(in) :: this
	integer :: j
	write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	write(OUTPUT_UNIT,*) '>>>> Exception in   --> ',trim(this%fctname)
	write(OUTPUT_UNIT,*) '>>>> Description    --> ',trim(this%message)
	select case (this%level)
	case (0); write(OUTPUT_UNIT,*) '>>>> Level          --> SUCCESS'
	case (1); write(OUTPUT_UNIT,*) '>>>> Level          --> WARNING'
	case (2); write(OUTPUT_UNIT,*) '>>>> Level          --> ERROR'
	end select
	if (associated(this%trace)) then
		do j=1,size(this%trace)
			write(OUTPUT_UNIT,*) '>>>> Passed through --> ',trim(this%trace(j))
		enddo
	endif
	write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	end subroutine printErrorMessage
!---------------------------------------------------------------
!> \brief Add a trace to error message
!
	subroutine addTraceErrorMessage(this,fctname)
	type (error_message) :: this
	character (len=*) :: fctname
	integer :: n
	if (associated(this%trace)) then
		n = size(this%trace)
		this%trace => reallocate(this%trace,n+1)
		this%trace(n+1) = trim(fctname)
	else
		allocate(this%trace(1))
		this%trace(1) = trim(fctname)
	endif
	end subroutine addTraceErrorMessage
!---------------------------------------------------------------
!> \brief Get level of error message
!
	integer function getLevelErrorMessage(this)
	type (error_message), intent(in) :: this
	getLevelErrorMessage = this%level
	end function getLevelErrorMessage
!
 end module errorMessage
