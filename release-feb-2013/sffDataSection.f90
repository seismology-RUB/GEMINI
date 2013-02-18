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
		type (sff_header), pointer :: header => null()                         ! header of data section
		type (sff_datablock), dimension(:), pointer :: datablock => null()     ! array of datablocks
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
!-------------------------------------------------------------
!> \brief create basic sff section
!
	subroutine createBasicSFFDataSection(this,ndbl)
	type (sff_datasection) :: this
	integer :: ndbl
!
	allocate(this%header)
	allocate(this%datablock(ndbl))
	this%ndbl = 0
	end subroutine createBasicSFFDataSection
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
!> \brief Add header to section
!
	subroutine addHeaderSFFDataSection(this,sffhead)
	type (sff_datasection) :: this
	type (sff_header) :: sffhead
	this%header = sffhead
	end subroutine addHeaderSFFDataSection
!-------------------------------------------------------------
!> \brief Add a datablock to sff data section
!
	subroutine addDataBlockSFFDataSection(this,sffdbl)
	type (sff_datasection) :: this
	type (sff_datablock) :: sffdbl
	if (this%ndbl+1 > size(this%datablock)) then
		print *,'Warning: addDataBlockSFFDataSection'
		print *,'Warning: try to add more datablocks than allocated'
		print *,'Warning: datablock not added!'
		return
	endif
	this%ndbl = this%ndbl+1
	this%datablock(this%ndbl) = sffdbl
	end subroutine addDataBlockSFFDataSection
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
!---------------------------------------------------------------------
!> \brief Write sff data section to filename
!
	subroutine writeSFFDataSection(this,lu,filename)
	type (sff_datasection) :: this
	integer :: lu,j
	character (len=*) :: filename
	open(lu,file = filename)
	call writeSFFHeader(this%header,lu,filename)
	do j = 1,this%ndbl
		call writeSFFDataBlock(this%datablock(j),lu,(j == this%ndbl))
	enddo
	end subroutine writeSFFDataSection
!
 end module sffDataSection
