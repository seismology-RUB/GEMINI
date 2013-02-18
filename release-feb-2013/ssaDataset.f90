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
!--------------------------------------------------------------------
!  module to read and write seismogram sections from and to 
!  stream access files
!  ssa: Seismograms using Stream Access
!--------------------------------------------------------------------
 module ssaDataset
	use streamAccess
	implicit none
!
!  source block
!
	type ssa_source
		character (len=1) :: cs        ! coordinate system
		real :: c1, c2, c3             ! location of source
		character (len=8) :: date      ! date of source (yyyymmdd)
		integer :: tfs                 ! start time to the full second after midnight
		integer :: tns                 ! fractional part of start time in nanoseconds
	end type
!
!  info block
!
	type ssa_info
		character (len=1) :: cs        ! coordinate system
		real :: c1, c2, c3             ! location of receiver
		integer :: nstack              ! number of stacks
	end type
!
!  data block
!
	interface dealloc; module procedure deallocSSADatablock; end interface
	interface operator (.nsamp.); module procedure getNsampSSADatablock; end interface
	interface operator (.dt.); module procedure getDtSSADatablock; end interface
	interface operator (.station.); module procedure getStationSSADatablock; end interface
	interface operator (.channel.); module procedure getChannelSSADatablock; end interface
	interface operator (.tanf.); module procedure getTanfSSADatablock; end interface
	interface operator (.trace.); module procedure traceSSADatablock; end interface
	interface operator (.therd.); module procedure getSourceTimeSSADatablock; end interface
	interface operator (.cs.); module procedure getCoordinateSystemSSADatablock; end interface
	interface operator (.date.); module procedure getDateSSADatablock; end interface
	interface operator (.tfs.); module procedure getTfsSSADatablock; end interface
	interface operator (.tns.); module procedure getTnsSSADatablock; end interface
	interface operator (.source.); module procedure getSourceSSADatablock; end interface
	interface operator (.info.); module procedure getInfoSSADatablock; end interface
	type ssa_datablock
		private
		integer :: n                         ! number of data 
		real :: dt                           ! sampling interval
		character (len=8) :: date            ! date of dirst data point (yyyymmdd)
		integer :: tfs                       ! time of first data point to the full second after midnight
		integer :: tns                       ! fractional part of start time in nanoseconds
		character (len=5) :: station         ! name of measuring site
		character (len=3) :: channel         ! name of channel
		real, dimension(:), pointer :: y => null()    ! pointer to data
		type (ssa_source) :: source          ! source block
		type (ssa_info) :: info              ! info block
		logical :: link                      ! just link to data
	end type
!
!  complete dataset
!
	interface new; module procedure readSSADataset; end interface
	interface dealloc; module procedure deallocSSADataset; end interface
	interface operator (.ndbl.); module procedure getNDatablocksSSADataset; end interface
	interface operator (.datablock.); module procedure getDatablockSSADataset; end interface
	type ssa_dataset
		private
		integer :: ndbl
		type (ssa_datablock), dimension(:), pointer :: datablock => null()
	end type
!
 contains
!------------------------------------------------------------------------------
!  basic constructor for source block, sets default values
!
	subroutine createSSASource(this,cs,c1,c2,c3,date,tfs,tns)
	type (ssa_source) :: this
	character (len=1), optional :: cs
	real, optional :: c1,c2,c3
	integer, optional :: tfs,tns
	character (len=8), optional :: date
	if (present(cs)) then;   this%cs = cs; else; this%cs = 'S'; endif
	if (present(c1)) then;   this%c1 = c1; else; this%c1 = 0.0; endif
	if (present(c2)) then;   this%c2 = c2; else; this%c2 = 0.0; endif
	if (present(c3)) then;   this%c3 = c3; else; this%c3 = 0.0; endif
	if (present(date)) then; this%date = date; else; this%date = '20090101'; endif
	if (present(tfs)) then;  this%tfs = tfs; else; this%tfs = 0; endif
	if (present(tns)) then;  this%tns = tns; else; this%tns = 0; endif
	end subroutine createSSASource
!------------------------------------------------------------------------------
!  read a source block from a SA dataset
!
	subroutine readSSASource(this,dset,fda)
	type (ssa_source) :: this
	type (data_stream_access) :: dset
	type (file_stream_access) :: fda
	type (flexible), dimension(:), pointer :: source
!
	call readDatasetVectorStreamAccess(dset,fda,source)
	this%cs = source(1)
	this%c1 = source(2); this%c2 = source(3); this%c3 = source(4)
	this%date = source(5); this%tfs = source(6); this%tns = source(7)
	deallocate(source)
	end subroutine readSSASource
!---------------------------------------------------------------------------
!  write a source block to a SA file
!
	subroutine writeSSASource(this,group,fda)
	type (ssa_source) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: group
	type (data_stream_access) :: dset
	type (flexible), dimension(:), allocatable :: source
	integer, dimension(1) :: dims
!
	dims(1) = 7
	allocate(source(7))
	source(1) = this%cs
	source(2) = this%c1; source(3) = this%c2; source(4) = this%c3
	source(5) = this%date; source(6) = this%tfs; source(7) = this%tns
	call createDatasetStreamAccess(dset,1,dims,T_FLEXIBLE)
	call writeDatasetVectorStreamAccess(dset,fda,source)
	call addDatasetStreamAccess(group,dset); call dealloc(dset)
	deallocate(source)
	end subroutine writeSSASource
!---------------------------------------------------------
!  read out source location
!
	subroutine locationSSASource(this,c1,c2,c3)
	type (ssa_source) :: this
	real :: c1,c2,c3
	c1 = this%c1; c2 = this%c2; c3 = this%c3
	end subroutine locationSSASource
!-----------------------------------------------------------------------
!  get double precision source time in seconds after some reference time
!
 	double precision function getSourceTimeSSASource(this,tref) result(tanf)
	type (ssa_source), intent(in) :: this
	integer, intent(in) :: tref
	tanf = dble(this%tfs-tref)+this%tns*1.d-9
	end function getSourceTimeSSASource
!-----------------------------------------------------------------------
!  get source time to the full second after midnight
!
	integer function getFullSecondSourceTimeSSASource(this) result(tfs)
	type (ssa_source), intent(in) :: this
	tfs = this%tfs
	end function getFullSecondSourceTimeSSASource
!-----------------------------------------------------------------------
!  get source time nano seconds
!
	integer function getNanoSecondSourceTimeSSASource(this) result(tns)
	type (ssa_source), intent(in) :: this
	tns = this%tns
	end function getNanoSecondSourceTimeSSASource

!-----------------------------------------------------------------------
!  get coordinate system of source
!
	function getCoordinateSystemSSASource(this) result(cs)
	type (ssa_source), intent(in) :: this
	character (len=1) :: cs
	cs = this%cs
	end function getCoordinateSystemSSASource
!-----------------------------------------------------------------------
!  get date of source
!
	function getDateSSASource(this) result(date)
	type (ssa_source), intent(in) :: this
	character (len=8) :: date
	date = this%date
	end function getDateSSASource
!------------------------------------------------------------------------------
!  basic constructor for info block, sets default values
!
	subroutine createSSAInfo(this,cs,c1,c2,c3,nstack)
	type (ssa_info) :: this
	character (len=1), optional :: cs
	real, optional :: c1,c2,c3
	integer, optional :: nstack
	if (present(cs)) then;   this%cs = cs; else; this%cs = 'S'; endif
	if (present(c1)) then;   this%c1 = c1; else; this%c1 = 1.0; endif
	if (present(c2)) then;   this%c2 = c2; else; this%c2 = 0.0; endif
	if (present(c3)) then;   this%c3 = c3; else; this%c3 = 0.0; endif
	if (present(nstack)) then; this%nstack = nstack; else; this%nstack = 1; endif
	end subroutine createSSAInfo
!------------------------------------------------------------------------------
!  read a info block from a SA dataset
!
	subroutine readSSAInfo(this,dset,fda)
	type (ssa_info) :: this
	type (data_stream_access) :: dset
	type (file_stream_access) :: fda
	type (flexible), dimension(:), pointer :: info
!
	call readDatasetVectorStreamAccess(dset,fda,info)
	this%cs = info(1)
	this%c1 = info(2); this%c2 = info(3); this%c3 = info(4); this%nstack = info(5)
	deallocate(info)
	end subroutine readSSAInfo
!---------------------------------------------------------------------------
!  write an info block to a SA file
!
	subroutine writeSSAInfo(this,group,fda)
	type (ssa_info) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: group
	type (data_stream_access) :: dset
	type (flexible), dimension(:), allocatable :: info
	integer, dimension(1) :: dims
!
	dims(1) = 5
	allocate(info(5))
	info(1) = this%cs
	info(2) = this%c1; info(3) = this%c2; info(4) = this%c3; info(5) = this%nstack
	call createDatasetStreamAccess(dset,1,dims,T_FLEXIBLE)
	call writeDatasetVectorStreamAccess(dset,fda,info)
	call addDatasetStreamAccess(group,dset); call dealloc(dset)
	deallocate(info)
	end subroutine writeSSAInfo
!-----------------------------------------------------------
!  return receiver location
!
	subroutine locationSSAInfo(this,x,y,z)
	type (ssa_info) :: this
	real :: x,y,z
	x = this%c1; y = this%c2; z = this%c3
	end subroutine locationSSAInfo
!--------------------------------------------------------------------------------
!  basic constructor for datablock
!  use only pointer to data
!
	subroutine createSSADatablock(this,n,dt,y,date,tfs,tns,station,channel,source,info)
	type (ssa_datablock) :: this
	integer :: n
	real :: dt
	real, dimension(:), target :: y
	character (len=*), optional :: date,station,channel
	integer, optional :: tfs,tns
	type (ssa_source), optional :: source
	type (ssa_info), optional :: info
!
	this%n = n; this%dt = dt; this%y => y; this%link = .true.
	if (present(date)) then; this%date = date; else; this%date = '20090101'; endif
	if (present(tfs)) then; this%tfs = tfs; else; this%tfs = 0; endif
	if (present(tns)) then; this%tns = tns; else; this%tns = 0; endif
	if (present(station)) then; this%station = station; else; this%station = 'NSP'; endif
	if (present(channel)) then; this%channel = channel; else; this%channel = 'XXZ'; endif
	if (present(source)) then; this%source = source; else; call createSSASource(this%source); endif
	if (present(info)) then; this%info = info; else; call createSSAInfo(this%info); endif
	end subroutine createSSADatablock
!------------------------------------------------------------------------------------
!  deallocate ssa datablock
!
	subroutine deallocSSADatablock(this)
	type (ssa_datablock) :: this
	if (associated(this%y)) then
		if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
	endif
	end subroutine deallocSSADatablock
!------------------------------------------------------------------------------------
!  read a datablock from SA file dataset
!  4 subdatasets: source, info, header, data
!  header: n,dt,date,tfs,tns,station,channel
!
	subroutine readSSADatablock(this,group,fda)
	type (ssa_datablock) :: this
	type (group_stream_access) :: group
	type (file_stream_access) :: fda
	type (group_stream_access), pointer :: pgr
	type (data_stream_access), pointer :: dset
	type (flexible), dimension(:), pointer :: fp
!
	call traversePathStreamAccess(group,0,(/ 1 /),pgr,dset)    ! source block
	call readSSASource(this%source,dset,fda)
!
	call traversePathStreamAccess(group,0,(/2/),pgr,dset)    ! info block
	call readSSAInfo(this%info,dset,fda)
!
	call traversePathStreamAccess(group,0,(/3/),pgr,dset)    ! header block
	call readDatasetVectorStreamAccess(dset,fda,fp)
	this%n = fp(1); this%dt = fp(2); this%date = fp(3); this%tfs = fp(4); this%tns = fp(5)
	this%station = fp(6); this%channel = fp(7)
	this%link = .false.
!
	call traversePathStreamAccess(group,0,(/4/),pgr,dset)    ! data block
	call readDatasetVectorStreamAccess(dset,fda,this%y)
!
	deallocate(fp)
	end subroutine readSSADatablock
!------------------------------------------------------------------------------
!  deep copy of a data block
!
	subroutine copySSADatablock(this,copy)
	type (ssa_datablock) :: this,copy
	integer :: n
	n = this%n
	copy = this
	nullify(copy%y)
	allocate(copy%y(n))
	copy%y(1:n) = this%y(1:n)
	copy%link = .false.
	end subroutine copySSADatablock
!---------------------------------------------------------------------------------------------
!  write a datablock to SA file
!  4 subdatasets: source, info, header, data
!  header: n,dt,date,tfs,tns,station,channel
!
	subroutine writeSSADatablock(this,group,fda)
	type (ssa_datablock) :: this
	type (group_stream_access) :: group
	type (file_stream_access) :: fda
	type (data_stream_access) :: dset
	type (flexible), dimension(:), pointer :: fp
	integer, dimension(1) :: dims
!
	call writeSSASource(this%source,group,fda)
	call writeSSAInfo(this%info,group,fda)
!
!  header block
!
	dims(1) = 7
	allocate(fp(7))
	fp(1) = this%n; fp(2) = this%dt; fp(3) = this%date; fp(4) = this%tfs; fp(5) = this%tns
	fp(6) = this%station; fp(7) = this%channel
	call createDatasetStreamAccess(dset,1,dims,T_FLEXIBLE)
	call writeDatasetVectorStreamAccess(dset,fda,fp)
	call addDatasetStreamAccess(group,dset); call dealloc(dset)
	deallocate(fp)
!
!  data
!
	dims(1) = this%n
	call createDatasetStreamAccess(dset,1,dims,T_REAL)
	call writeDatasetVectorStreamAccess(dset,fda,this%y(1:this%n))
	call addDatasetStreamAccess(group,dset); call dealloc(dset)
	end subroutine writeSSADatablock
!------------------------------------------------------------------------------
!  return a pointer to the data
!
	function traceSSADatablock(this)
	type (ssa_datablock), intent(in) :: this
	real, dimension(:), pointer :: traceSSADatablock
!
	traceSSADatablock => this%y
	end function traceSSADatablock
!------------------------------------------------------------------------------
!  get dt of datablock
!
	function getDtSSADatablock(this)
	type (ssa_datablock), intent(in) :: this
	real :: getDtSSADatablock
	getDtSSADataBlock = this%dt
	end function getDtSSADatablock
!------------------------------------------------------------------------------
!  get station name of datablock
!
	character (len=5) function getStationSSADatablock(this)
	type (ssa_datablock), intent(in) :: this
	getStationSSADatablock = this%station
	end function getStationSSADatablock
!------------------------------------------------------------------------------
!  get channel of datablock
!
	character (len=3) function getChannelSSADatablock(this)
	type (ssa_datablock), intent(in) :: this
	getChannelSSADatablock = this%channel
	end function getChannelSSADatablock
!------------------------------------------------------------------------------
!  get start time of data block with respect to some full second reference time
!
	double precision function getTanfSSADatablock(this,tref)
	type (ssa_datablock), intent(in) :: this
	integer, intent(in) :: tref
	getTanfSSADatablock = dble(this%tfs-tref)+this%tns*1.d-9
	end function getTanfSSADatablock
!------------------------------------------------------------------------------
!  get source time of data block with respect to some full second reference time
!
	function getSourceTimeSSADatablock(this,tref) result(tanf)
	type (ssa_datablock), intent(in) :: this
	integer, intent(in) :: tref
	double precision :: tanf
	tanf = getSourceTimeSSASource(this%source,tref)
	end function getSourceTimeSSADatablock
!------------------------------------------------------------------------------
!  get nsamp of datablock 
!
	integer function getNsampSSADatablock(this)
	type (ssa_datablock), intent(in) :: this
	getNsampSSADatablock = this%n
	end function getNsampSSADatablock
!------------------------------------------------------------------------------
!  get source location
!
	subroutine getSourceLocationSSADatablock(this,c1,c2,c3)
	type (ssa_datablock), intent(in) :: this
	real, intent(out) :: c1,c2,c3 
	call locationSSASource(this%source,c1,c2,c3)
	end subroutine getSourceLocationSSADatablock
!------------------------------------------------------------------------------
!  get receiver location
!
	subroutine getReceiverLocationSSADatablock(this,c1,c2,c3)
	type (ssa_datablock), intent(in) :: this
	real, intent(out) :: c1,c2,c3 
	call locationSSAInfo(this%info,c1,c2,c3)
	end subroutine getReceiverLocationSSADatablock
!------------------------------------------------------------------------------
!  get coordinate system of datablock
!
	function getCoordinateSystemSSADatablock(this) result(cs)
	type (ssa_datablock), intent(in) :: this
	character (len=1) :: cs
	cs = getCoordinateSystemSSASource(this%source)
	end function getCoordinateSystemSSADatablock
!------------------------------------------------------------------------------
!  get date of datablock
!
	function getDateSSADatablock(this) result(date)
	type (ssa_datablock), intent(in) :: this
	character (len=8) :: date
	date = this%date
	end function getDateSSADatablock
!------------------------------------------------------------------------------
!  get full starttime of data block to the full second after midnight
!
	function getTfsSSADatablock(this) result(tfs)
	type (ssa_datablock), intent(in) :: this
	integer :: tfs
	tfs = this%tfs
	end function getTfsSSADatablock
!------------------------------------------------------------------------------
!  get start time nano seconds
!
	function getTnsSSADatablock(this) result(tns)
	type (ssa_datablock), intent(in) :: this
	integer :: tns
	tns = this%tns
	end function getTnsSSADatablock
!------------------------------------------------------------------------------
!  get source
!
	function getSourceSSADatablock(this) result(source)
	type (ssa_datablock), intent(in) :: this
	type (ssa_source) :: source
	source = this%source
	end function getSourceSSADatablock
!------------------------------------------------------------------------------
!  get info
!
	function getInfoSSADatablock(this) result(info)
	type (ssa_datablock), intent(in) :: this
	type (ssa_info) :: info
	info = this%info
	end function getInfoSSADatablock
!-------------------------------------------------------------------------------
!  create basic SSA dataset
!
	subroutine createBasicSSADataset(this,ndbl)
	type (ssa_dataset) :: this
	integer :: ndbl
	allocate(this%datablock(ndbl))
	this%ndbl = 0
	end subroutine createBasicSSADataset
!------------------------------------------------------------------------------
!  add a SSA datablock to SSA dataset
!
	subroutine addDatablockSSADataset(this,dbl)
	type (ssa_dataset) :: this
	type (ssa_datablock) :: dbl
	integer :: n
	n = this%ndbl+1
	if (n > size(this%datablock)) then
		print *,'no space left in ssa dataset'
		stop
	endif
	call copySSADatablock(dbl,this%datablock(n))     ! true copy of datablock
	this%ndbl = n
	end subroutine addDatablockSSADataset
!-------------------------------------------------------------------------------
!  read a complete ssa dataset from stream access file
!
	subroutine readSSADataset(this,lu,filename,vsflag)
	type (ssa_dataset) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: root
	type (group_stream_access), pointer :: group
	type (data_stream_access), pointer :: dset
	integer :: lu
	character (len=*) filename
	integer :: ierr,ndbl,j
	logical, optional :: vsflag
!
	if (present(vsflag)) then; verboseStreamAccess = vsflag; else; verboseStreamAccess = .false.; endif
	ierr = openFileStreamAccess(fda,lu,filename)
	if (ierr /= 0) then
		print *,'Stream Access file: ',trim(filename),' could not be opened!'
		stop
	endif
!
	call readGroupStreamAccess(root,fda)
	ndbl = .nsubgroup.root
	print *,'readSSADataset: Number of datablocks in dataset: ',ndbl
	if (ndbl == 0) stop
	allocate(this%datablock(ndbl))
	do j=1,ndbl
		call traversePathStreamAccess(root,1,(/j,0/),group,dset)
		call readSSADatablock(this%datablock(j),group,fda)
	enddo
	this%ndbl = ndbl
	call clearGroupTree(root); call dealloc(fda)
	end subroutine readSSADataset
!-------------------------------------------------------------------------------
!> \brief Get a SSA dataset from a streamAccess file
!> \param this ssa_dataset to be created
!> \param fda stream access file
!> \param group Group within which ssa_dataset has been written
!
	subroutine getSSADataset(this,fda,group)
	type (ssa_dataset) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: group
	type (group_stream_access), pointer :: g
	type (data_stream_access), pointer :: dset
	integer :: ndbl,j
!
	ndbl = .nsubgroup.group
	if (ndbl == 0) then
		print *,'no datablocks found in ',trim(.groupname.group),.grouptag.group
		stop
	endif
	allocate(this%datablock(ndbl))
	do j=1,ndbl
		call traversePathStreamAccess(group,1,(/ j /),g,dset)
		call readSSADatablock(this%datablock(j),g,fda)
	enddo
	this%ndbl = ndbl
	end subroutine getSSADataset
!-------------------------------------------------------------------
!  deallocate a ssa dataset
!
	subroutine deallocSSADataset(this)
	type (ssa_dataset) :: this
	integer :: j
	if (associated(this%datablock)) then
		do j=1,size(this%datablock)
			call dealloc(this%datablock(j))
		enddo
	endif
	end subroutine deallocSSADataset
!-------------------------------------------------------------------------------
!  get number of datablocks in dataset
!
	integer function getNDatablocksSSADataset(this)
	type (ssa_dataset), intent(in) :: this
	getNDatablocksSSADataset = this%ndbl
	end function getNDatablocksSSADataset
!-------------------------------------------------------------
!  return a pointer to the k-th datablock
!
	function getDatablockSSADataset(this,k)
	type (ssa_dataset), intent(in) :: this
	integer, intent(in) :: k
	type (ssa_datablock), pointer :: getDatablockSSADataset
	getDatablockSSADataset => this%datablock(k)
	end function getDatablockSSADataset
!------------------------------------------------------------------
!  add a ssa dataset to a stream access file
!  under group
!
	subroutine addSSADataset(this,fda,group)
	type (ssa_dataset) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: group
	type (group_stream_access) :: dblgroup
	integer :: j
!
	do j=1,this%ndbl
		call createGroupStreamAccess(dblgroup,'Datablock',j)
		call writeSSADatablock(this%datablock(j),dblgroup,fda)
		call addSubgroupStreamAccess(group,dblgroup); call dealloc(dblgroup)
	enddo
	end subroutine addSSADataset
!-------------------------------------------------------------------
!  write a complete ssa dataset to stream access file
!
	subroutine writeSSADataset(this,lu,filename,vsflag)
	type (ssa_dataset) :: this
	type (file_stream_access) :: fda
	type (group_stream_access) :: root
	type (group_stream_access) :: dblgroup
	integer :: lu
	character (len=*) filename
	integer :: ierr,ndbl,j
	logical, optional :: vsflag
!
	if (present(vsflag)) then; verboseStreamAccess = vsflag; else; verboseStreamAccess = .false.; endif
	if (associated(this%datablock)) then; ndbl = size(this%datablock); else; ndbl = 0; endif
	print *,'writeSSADataset: Number of datablocks in dataset: ',ndbl
	if (ndbl == 0) stop
!
	ierr = createFileStreamAccess(fda,lu,filename)
	if (ierr /= 0) then
		print *,'Stream Access file: ',trim(filename),' could not be created!'
		stop
	endif
!
	call createGroupStreamAccess(root,'Root',0)
	do j=1,ndbl
		call createGroupStreamAccess(dblgroup,'Datablock',j)
		call writeSSADatablock(this%datablock(j),dblgroup,fda)
		call addSubgroupStreamAccess(root,dblgroup); call dealloc(dblgroup)
	enddo
	call writeGroupStreamAccess(root,fda)
	call clearGroupTree(root)
	call dealloc(fda)
	end subroutine writeSSADataset
!
 end module
