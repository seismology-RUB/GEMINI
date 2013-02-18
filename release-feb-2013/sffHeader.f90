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
!--------------------------------------------------------
!  module to describe a SFF header
!--------------------------------------------------------
 module sffHeader
	use sffFree
	use sffSource
	implicit none
	private readSFFHeader, destroySFFHeader
	interface new
		module procedure readSFFHeader
		module procedure basicSFFHeader
	end interface
	interface dealloc
		module procedure destroySFFHeader
	end interface dealloc
	interface operator (.cs.); module procedure getCoordinateSystemSFFHeader; end interface
	interface operator (.therd.); module procedure getSourceTimeSFFHeader; end interface
	type sff_header
		private
		real :: version                           ! Version of SFF-library 
		character (len=13) :: timestamp           ! Time the file was written
		character (len=10) :: code                ! Code indicating optional blocks
	 	type (sff_free) :: free                   ! structure for free block
		type (sff_source) :: source               ! structure for source line
	end type sff_header
!
 contains
!----------------------------------------------------------
!  read header from SFF file
!
	subroutine readSFFHeader(this,lu,filename)
	type (sff_header) :: this
	integer :: lu,ierr,j
	character (len=*) :: filename
!
	open(lu, file=filename, err=99, status='old')
	call sff_RStatus(lu,this%version,this%timestamp,this%code,ierr)
!
!  read free block and source line if present
! 
	j=1
	do, while(this%code(j:j) /= ' ')
		if(this%code(j:j) == 'F') then
			call new(this%free,lu)
		endif
		if(this%code(j:j) == 'S') then
			call new(this%source,lu)
		endif
		j=j+1
	enddo
	return
 99	print *,'ERROR: opening file: ',filename
	end subroutine readSFFHeader
!----------------------------------------------------------------
!  construct SFF header from scratch without free and source line
!
	subroutine basicSFFHeader(this)
	type (sff_header) :: this
	real sff_libversion
	character (len=8) :: date
	character (len=10) :: time
!
	this%version = sff_libversion()
	call date_and_time(date,time)
	this%timestamp = date(3:8)//'.'//time(1:6)
	this%code = ''
	end subroutine basicSFFHeader
!------------------------------------------------------------
!  destructor
!
	subroutine destroySFFHeader(this)
	type (sff_header) :: this
	integer :: j
!
	j=1
	do, while(this%code(j:j) /= ' ')
		if(this%code(j:j) == 'F') call dealloc(this%free)
		j=j+1
	enddo
	end subroutine destroySFFHeader
!-------------------------------------------------------------
!  read out free block
!
	function freeLinesSFFHeader(this)
	type (sff_header) :: this
	character (len=80), dimension(:), pointer :: freeLinesSFFHeader
!
	freeLinesSFFHeader => linesSFFFree(this%free)
	end function freeLinesSFFHeader
!----------------------------------------------------------------
!  add source line to header
!
	subroutine addSourceSFFHeader(this,typ,cs,xs,ys,zs,date,time)
	type (sff_header) :: this
	character typ*(*),cs*1,date*6,time*10,typ2*20
	real :: xs,ys,zs
	integer :: i
!
	typ2 = trim(typ)
	do i=len_trim(typ)+1,20; typ2(i:i) = ' '; enddo
	if(cs == 'S') call new(this%source,typ2,cs,90.-xs,ys,zs,date,time)
	if(cs == 'C') call new(this%source,typ2,cs,xs,ys,zs,date,time)
	this%code = trim(this%code)//'S'
	end subroutine addSourceSFFHeader
!-------------------------------------------------------------
!  get source location
!  set to (0,0,1) if source line is not available
!
	subroutine getSourceLocationSFFHeader(this,x,y,z)
	type (sff_header) :: this
	real :: x,y,z
!
	if (index(this%code,'S') > 0) then
		call locationSFFSource(this%source,x,y,z)
	else
		x = 0.; y = 0.; z = 1.
	endif
	end subroutine getSourceLocationSFFHeader
!-------------------------------------------------------------
!  get coordinate system
!
	character (len=1) function getCoordinateSystemSFFHeader(this)
	type (sff_header), intent(in) :: this
	character :: cs*1, typ*20
!
	if (index(this%code,'S') > 0) then
		call descriptionSFFSource(this%source,typ,cs)
	else
		cs = 'C'
	endif
	getCoordinateSystemSFFHeader = cs
	end function getCoordinateSystemSFFHeader
!------------------------------------------------------------
!  get source time
!
	real function getSourceTimeSFFHeader(this)
	type (sff_header), intent(in) :: this
	if (index(this%code,'S') > 0) then
		getSourceTimeSFFHeader = getSourceTimeSFFSource(this%source)
	else 
		getSourceTimeSFFHeader = 0.
	endif
	end function getSourceTimeSFFHeader
!-------------------------------------------------------------
!  write a SFF header to file as is
!
	subroutine writeSFFHeader(this,lu,filename)
	type (sff_header) :: this
	integer :: lu,ierr,j
	character (len=*) :: filename
!
	call sff_New(lu,filename,ierr)    ! deletes an existing file of the same name
	open(lu,file=filename)
	call sff_WStatus(lu,this%code)
!
!  write free block and source line if present
! 
	j=1
	do, while(this%code(j:j) /= ' ')
		if(this%code(j:j) == 'F') then
			call writeSFFFree(this%free,lu)
		endif
		if(this%code(j:j) == 'S') then
			call writeSFFSource(this%source,lu)
		endif
		j=j+1
	enddo
	end subroutine writeSFFHeader
!-------------------------------------------------------------
 end module sffHeader
