!------------------------------------------------------------
!  Module that describes a complete SFF data section
!------------------------------------------------------------
 module sffDataSection
	use sffHeader
	use sffDatablock
	implicit none
	interface new;     module procedure createSFFDataSection;  end interface
	interface dealloc; module procedure deallocSFFDataSection; end interface
	interface operator (.ndbl.); module procedure getNDatablocksSFFDataSection; end interface
	interface operator (.header.); module procedure getHeaderSFFDataSection; end interface
	interface operator (.datablock.); module procedure getDataBlockSFFDataSection; end interface
	type sff_datasection
		private
		integer :: ndbl
		type (sff_header), pointer :: header                          ! header of data section
		type (sff_datablock), dimension(:), pointer :: datablock      ! array of datablocks
	end type
!
 contains
!------------------------------------------------------------
!  read in the data section
!
	subroutine createSFFDataSection(this,lu,filename)
	type (sff_datasection) :: this
	integer :: lu, n
	character (len=*) :: filename
	logical last
!
	allocate(this%header)
	call new(this%header,lu,filename)
	last = .false.
	n=1
	allocate(this%datablock(1))
	do while (.not. last)
		this%datablock => reallocateSFFDatablock(this%datablock,n)
		call new(this%datablock(n),lu,last)
		n=n+1
	enddo
	this%ndbl = n-1
	end subroutine createSFFDataSection
!------------------------------------------------------------
!  deallocate complete section
!
	subroutine deallocSFFDataSection(this)
	type (sff_datasection) :: this
	integer :: i
!
	if(associated(this%header)) then
		call dealloc(this%header)
		deallocate(this%header)
	endif
	if (associated(this%datablock)) then
		do i=1,this%ndbl
			call dealloc(this%datablock(i))
		enddo
		deallocate(this%datablock)
	endif
	end subroutine deallocSFFDataSection
!-------------------------------------------------------------
!  calculate total number of samples
!
	integer function ntotSamplesSFFDataSection(this)
	type (sff_datasection), intent(in) :: this
	integer ntot,i
!
	ntot = 0
	do i=1,this%ndbl
		ntot = ntot+.nsamp.(this.datablock.i)
	enddo
	ntotSamplesSFFDataSection = ntot
	end function ntotSamplesSFFDataSection
!-------------------------------------------------------------
!  calculate maximum number of samples in a trace
!
	integer function nmaxSamplesSFFDataSection(this)
	type (sff_datasection), intent(in) :: this
	integer :: maxsamp,i
!
	maxsamp = 0
	do i=1,this%ndbl
		maxsamp = max(maxsamp,.nsamp.(this.datablock.i))
	enddo
	nmaxSamplesSFFDataSection = maxsamp
	end function nmaxSamplesSFFDataSection
!-------------------------------------------------------------
!  return number of datablocks
!
	integer function getNDatablocksSFFDataSection(this)
	type (sff_datasection), intent(in) :: this
	getNDatablocksSFFDataSection = this%ndbl
	end function getNDatablocksSFFDataSection
!-------------------------------------------------------------
!  return a pointer to the header
!
	function getHeaderSFFDataSection(this)
	type (sff_datasection), intent(in) :: this
	type (sff_header), pointer :: getHeaderSFFDataSection
	getHeaderSFFDataSection => this%header
	end function getHeaderSFFDataSection
!-------------------------------------------------------------
!  return a pointer to the k-th datablock
!
	function getDatablockSFFDataSection(this,k)
	type (sff_datasection), intent(in) :: this
	integer, intent(in) :: k
	type (sff_datablock), pointer :: getDatablockSFFDataSection
	getDatablockSFFDataSection => this%datablock(k)
	end function getDatablockSFFDataSection
!
 end module sffDataSection
