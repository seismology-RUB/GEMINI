!----------------------------------------------------------------
!> \brief Support tools for using MPI

!> \author Wolfgang Friederich

!> \par Description
!>  Module to ease MPI parallelization of code 
!!  Creates an object that contains information about the parallel
!!  environment.
!<---------------------------------------------------------------
 module mpiSupport
       use mpi
	implicit none
	interface new; module procedure createMpiSupport; end interface
	interface dealloc; module procedure deallocMpiSupport; end interface
	interface operator (.myrank.); module procedure myrankMpiSupport; end interface
	interface operator (.numtasks.); module procedure numtasksMpiSupport; end interface
	type mpi_support
		private
		integer myrank                           ! Rank of this process
		integer numtasks                         ! Number of parallel tasks
	end type
!
 contains
!-----------------------------------------------------------------
!> \brief  Initialize MPI
!
	subroutine createMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr,rc
!
	call mpi_init(ierr)
	if (ierr .ne. MPI_SUCCESS) then
		print *,'Error starting MPI program. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
	endif
!
	call MPI_COMM_RANK(MPI_COMM_WORLD, this%myrank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, this%numtasks, ierr)
	end subroutine createMpiSupport
!-----------------------------------------------------------------
!> \brief Dealloc support object and finalize MPI
!
	subroutine deallocMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr
	call MPI_FINALIZE(ierr)
	end subroutine deallocMpiSupport
!-----------------------------------------------------------------
!> \brief get my rank
!
	integer function myrankMpiSupport(this)
	type (mpi_support), intent(in) :: this
	myrankMpiSupport = this%myrank
	end function myrankMpiSupport
!-----------------------------------------------------------------
!> \brief get number of parallel taks
!
	integer function numtasksMpiSupport(this)
	type (mpi_support), intent(in) :: this
	numtasksMpiSupport = this%numtasks
	end function numtasksMpiSupport
!
 end module
