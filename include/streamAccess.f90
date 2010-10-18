!----------------------------------------------------------------
!> \brief Handle structured stream input and output

!> \author Wolfgang Friederich

!> \par Description
!>  Support for stream access file reading and writing.
!!  File tree starts with a root group that may contain further groups
!!  and datasets. Each subgroup can contain subgroups and datasets
!!  itself. Idea is to create groups and datasets and to build up a
!!  tree under the root group which contains links to all subgroups
!!  and datasets. Datasets are written to file, file positions are stored
!!  in the group and dataset objects. Finally, the complete tree
!!  with information about the data is also written to file in a way
!!  that by reading the file the complete information and data can
!!  be retrieved.
!<---------------------------------------------------------------
 module streamAccess
	use flexibleType
	use kindDefinitions
	implicit none
	interface dealloc
		module procedure deallocFileStreamAccess
		module procedure nullifyGroupStreamAccess
		module procedure fullyDeallocDatasetStreamAccess
	end interface
	interface clearGroupTree; module procedure recursiveDeallocGroupStreamAccess; end interface
	interface writeDatasetVectorStreamAccess
		module procedure writeDatasetIntegerVectorStreamAccess
		module procedure writeDatasetRealVectorStreamAccess
		module procedure writeDatasetDoubleVectorStreamAccess
		module procedure writeDatasetComplexVectorStreamAccess
		module procedure writeDatasetDoubleComplexVectorStreamAccess
		module procedure writeDatasetCharVectorStreamAccess
		module procedure writeDatasetFlexibleVectorStreamAccess
	end interface
	interface readDatasetVectorStreamAccess
		module procedure readDatasetIntegerVectorStreamAccess
		module procedure readDatasetRealVectorStreamAccess
		module procedure readDatasetDoubleVectorStreamAccess
		module procedure readDatasetComplexVectorStreamAccess
		module procedure readDatasetDoubleComplexVectorStreamAccess
		module procedure readDatasetCharVectorStreamAccess
		module procedure readDatasetFlexibleVectorStreamAccess
	end interface
	interface operator (.ndataset.); module procedure getNdsGroupStreamAccess; end interface
	interface operator (.nsubgroup.); module procedure getNsubGroupStreamAccess; end interface
!
	type file_stream_access
		private
		integer :: lu                            !< Fortran identifier of file
		integer (longint) :: current_file_pos    !< current position of file pointer
		character (len=132) :: filename          !< name of file
	end type
!
	type data_stream_access
!		private
		integer (longint) :: filepos_start                 !< position where data start
		integer :: rank                                    !< dimensionality of data array
		integer, dimension(:), pointer :: dims =>  null()  !< number of data
		integer :: datatype                                !< integer coded primitive type of data
	end type
!
	type group_stream_access
!		private
		integer (longint) :: filepos      !< position where to find subgroup and data information in file
		integer :: tag                    !< positive integer number indexing the group in some way (root has zero!)
		character (len=132) :: name       !< name specifying group properties
		type (group_stream_access), pointer :: parent => null()  !< pointer to parent
		type (group_stream_access), dimension(:), pointer :: subgroup => null()   !< pointer to subgroups
		type (data_stream_access), dimension(:), pointer :: dataset => null()     !< pointer to datasets
	end type
!
	logical :: verboseStreamAccess = .false.
!
 contains
!---------------------------------------------------------------------------
!> \brief Create a new stream access file object. Open file for writing
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
	integer function createFileStreamAccess(this,lu,filename) result(ios)
	type (file_stream_access) :: this
	integer :: lu
	character (len=*) :: filename
	logical :: exflag
!
	inquire(file=filename, exist = exflag)   ! delete file if it exists
	if (exflag) then
		open(lu,file=filename,access = 'stream', form='unformatted')
		close(lu,status = 'delete')
	endif
	open(lu,file=filename,form='unformatted',access='stream',status='new',iostat=ios)
	if (ios /= 0) return
	this%lu = lu
	this%current_file_pos = 1+bit_size(this%current_file_pos)/8      !  leave space for a long integer
	this%filename = filename
	end function createFileStreamAccess
!-----------------------------------------------------------------
!> \brief Open an existing stream access file for reading
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
	integer function openFileStreamAccess(this,lu,filename) result(ios)
	type (file_stream_access) :: this
	integer :: lu
	character (len=*) :: filename
	open(lu,file=filename,form='unformatted',access='stream',status='old',iostat=ios)
	if (ios /= 0) return
	this%lu = lu
	read(lu,pos = 1) this%current_file_pos
	end function openFileStreamAccess
!-----------------------------------------------------------------
!> \brief Deallocate a direct access file object and close file
!> \param this file_stream_access object
!
	subroutine deallocFileStreamAccess(this)
	type (file_stream_access) :: this
	close(this%lu)
	end subroutine deallocFileStreamAccess
!----------------------------------------------------------------
!> \brief Get the logical unit of the file
!> \param this file_stream_access object
!
	integer function getFileUnitStreamAccess(this)
	type (file_stream_access) :: this
	getFileUnitStreamAccess = this%lu
	end function getFileUnitStreamAccess
!----------------------------------------------------------------
!> \brief Get the current file position
!> \param this file_stream_access object
!
	function getCurrentFilePositionStreamAccess(this) result(filepos)
	type (file_stream_access) :: this
	integer (longint) :: filepos
	filepos = this%current_file_pos
	end function getCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Increment current file position
!> \param this file_stream_access object
!> \param n number of bytes to move
!
	function IncrementCurrentFilePositionStreamAccess(this,n) result(filepos)
	type (file_stream_access) :: this
	integer (longint) :: filepos
	integer :: n
	this%current_file_pos = this%current_file_pos + n
	filepos = this%current_file_pos
	end function IncrementCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Set current file position
!> \param this file_stream_access object
!> \param filepos desired position of file pointer
!
	subroutine setCurrentFilePositionStreamAccess(this,filepos)
	type (file_stream_access) :: this
	integer (longint) :: filepos
	this%current_file_pos = filepos
	end subroutine setCurrentFilePositionStreamAccess
!---------------------------------------------------------------------
!> \brief Create a group. Contents to be filled in later.
!> \param this group_stream_access object
!> \param name name of the group
!> \param tag an identifier of the group
!
	subroutine createGroupStreamAccess(this,name,tag)
	type (group_stream_access) :: this
	character (len=*) :: name
	integer :: tag
!
	this%name = name
	this%tag = tag
	if (verboseStreamAccess) then
		print *,'create group named ',trim(this%name),' Tag = ',this%tag
	endif
!
	end subroutine createGroupStreamAccess
!---------------------------------------------------------------------
!> \brief Recursively deallocate group including subgroups and datasets
!> \param this group_stream_access object
!
	recursive subroutine recursiveDeallocGroupStreamAccess(this)
	type (group_stream_access) :: this
	integer :: j
	if (associated(this%subgroup)) then
		do j=1,size(this%subgroup)
			call recursiveDeallocGroupStreamAccess(this%subgroup(j))
		enddo
		deallocate(this%subgroup)
	endif
	if (associated(this%dataset)) then
		do j=1,size(this%dataset)
			call fullyDeallocDatasetStreamAccess(this%dataset(j))
		enddo
		deallocate(this%dataset)
	endif
	end subroutine recursiveDeallocGroupStreamAccess
!-----------------------------------------------------------------------
!> \brief Only nullify pointers in group object. Do not deallocate.
!> \param this group_stream_access object
!
	subroutine nullifyGroupStreamAccess(this)
	type (group_stream_access) :: this
	if (associated(this%parent)) nullify(this%parent)
	if (associated(this%subgroup)) nullify(this%subgroup)
	if (associated(this%dataset)) nullify(this%dataset)
	end subroutine nullifyGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a subgroup to this group.
!> \param this group_stream_access object
!> \param group subgroup object to be incorporated into the tree
!> \par
!> The subgroup must be completely filled with contents.
!! Only a SHALLOW copy of the group object into the tree is performed.
!! Avoid copying the same group into different trees, because
!! clearing one tree leads to a deallocation of memory still associated
!! with pointers in other trees. Clearing these trees will then fail. 
!! After calling this routine \a group may be deallocated using dealloc. 
!<
	subroutine addSubgroupStreamAccess(this,group)
	type (group_stream_access), target :: this
	type (group_stream_access) :: group
	integer :: n,j
!
	if (associated(this%subgroup)) then
		n = size(this%subgroup)
	else
		n = 0
	endif
	this%subgroup => reallocateGroupStreamAccess(this%subgroup,n+1)
	this%subgroup(n+1) = group
	if (verboseStreamAccess) then
		print *,'add group named ',trim(this%subgroup(n+1)%name),' as ',n+1,' th subgroup of group ', &
		      & trim(this%name),' with tag = ',this%tag
	endif
	this%subgroup(n+1)%parent => this
!
!  redirect the parent member of group's subgroups from group to this%subgroup(n+1)
!
	if (associated(group%subgroup)) then
		do j=1,size(group%subgroup)
			group%subgroup(j)%parent => this%subgroup(n+1)
		enddo
	endif
	end subroutine addSubgroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a dataset to some existing group
!> \param this group_stream_access object
!> \param dset dataset object
!> \par 
!> The dataset object must contain all information. This implies
!! that the data have already been written to file
!<
	subroutine addDatasetStreamAccess(this,dset)
	type (group_stream_access) :: this
	type (data_stream_access), target :: dset
	integer :: n
!
	if (associated(this%dataset)) then
		n = size(this%dataset)
	else
		n = 0
	endif
	this%dataset => reallocateDatasetStreamAccess(this%dataset,n+1)
	this%dataset(n+1) = dset
	allocate(this%dataset(n+1)%dims(1))
	this%dataset(n+1)%dims(1) = dset%dims(1)
	if (verboseStreamAccess) then
		print *,'add ',n+1,' th dataset with size ',this%dataset(n+1)%dims,' to group ',trim(this%name),' with tag ',this%tag
	endif
	end subroutine addDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of datasets contained in group
!> \param this group_stream_access object
!
	integer function getNdsGroupStreamAccess(this)
	type (group_stream_access), intent(in) :: this
	if (associated(this%dataset)) then
		getNdsGroupStreamAccess = size(this%dataset)
	else
		getNdsGroupStreamAccess = 0
	endif
	end function getNdsGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of subgroups contains in group
!> \param this group_stream_access object
!
	integer function getNsubGroupStreamAccess(this)
	type (group_stream_access), intent(in) :: this
	if (associated(this%subgroup)) then
		getNsubGroupStreamAccess = size(this%subgroup)
	else
		getNsubGroupStreamAccess = 0
	endif
	end function getNsubGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get name of group
!> \param this group_stream_access object
!
	character (len=132) function getNameGroupStreamAccess(this)
	type (group_stream_access), intent(in) :: this
	getNameGroupStreamAccess = this%name
	end function getNameGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get a pointer to the parent group
!> \param this group_stream_access object
!> \return a pointer to the parent group
!
	function getParentStreamAccess(this) result(parent)
	type (group_stream_access) :: this
	type (group_stream_access), pointer :: parent
	parent => this%parent
	end function getParentStreamAccess
!----------------------------------------------------------------------
!> \brief Get pointer to selected dataset in group
!> \param this group_stream_access object
!> \param k index of selected dataset in group
!> \return a pointer to a data_stream_access object
!
	function getDatasetSelectedStreamAccess(this,k) result(dset)
	type (group_stream_access) :: this
	type (data_stream_access), pointer:: dset
	integer :: k
!
	dset => this%dataset(k)
	end function getDatasetSelectedStreamAccess
!----------------------------------------------------------------------
!> \brief Recursively write the group data to file
!> \param this file_stream_access object
!> \param fda file_stream_access object
!
	recursive subroutine writeGroupStreamAccess(this,fda)
	type (group_stream_access) :: this
	type (file_stream_access) :: fda
	integer :: nsub,nds,j,lu,move
	integer (longint) :: filepos
!
	if (associated(this%subgroup)) then    ! first work through my subgroups
		do j=1,size(this%subgroup)               
			call writeGroupStreamAccess(this%subgroup(j),fda)
		enddo
	endif
!
!  now deal with myself
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
!
	this%filepos = filepos           ! position where infos about my group start
	                                 ! this info will later be written to file by my parent group
!
!  if I am the root group I write now the first position of the file
!  specifiying the position where info about the root group's contents has been written
!
	if (.not. associated(this%parent)) write(lu,pos = 1) filepos
!
	if (associated(this%subgroup)) then; nsub = size(this%subgroup); else; nsub = 0; endif
	if (associated(this%dataset)) then; nds = size(this%dataset); else; nds = 0; endif
!
	write(lu,pos=filepos) this%name,this%tag,nsub,nds   ! write name, tag etc
	if (verboseStreamAccess) then
		print *,'write group: ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
	endif
	filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
!
	if (nsub > 0) then                                  ! write subgroup record info
		write(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
		move = nsub*bit_size(filepos)/8
		filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
	endif
	if (nds > 0) then                                   ! write dataset record info
		write(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
		move = nds*bit_size(filepos)/8
		filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
	endif
	end subroutine writeGroupStreamAccess
!------------------------------------------------------------------------------------
!> \brief Recursively read group data from file into memory starting with root group
!> \param this group_stream_access object
!> \param fda file from which group data read
!
	recursive subroutine readGroupStreamAccess(this,fda)
	type (group_stream_access), target :: this
	type (file_stream_access), target :: fda
	integer :: lu,nsub,nds,j,move
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = getCurrentFilePositionStreamAccess(fda)
	this%filepos = filepos
	read(lu,pos=filepos) this%name,this%tag,nsub,nds      ! read first group record
	if (verboseStreamAccess) then
		print *,'read group ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
	endif
	filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
!
	if (nsub > 0) then                                   ! read subgroup record info
		allocate(this%subgroup(nsub))
		do j=1,nsub
			this%subgroup(j)%parent => this
		enddo
		read(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
		move = nsub*bit_size(filepos)/8
		filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
	else
		this%subgroup => null()
	endif
	if (nds > 0) then                                    ! read dataset record info
		allocate(this%dataset(nds))
		read(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
		move = nds*bit_size(filepos)/8
		filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
	else
		this%dataset => null()
	endif
	if (nsub > 0) then                                   !  read in subgroups
		do j=1,nsub
			call setCurrentFilePositionStreamAccess(fda,this%subgroup(j)%filepos)
			call readGroupStreamAccess(this%subgroup(j),fda)   ! read in subgroups
		enddo
	endif
	end subroutine readGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a group array
!> \param p pointer to an array of groups
!> \param n number of group objects to be allocated
!> \return newp pointer to the new array of groups
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
	function reallocateGroupStreamAccess(p, n) result(newp)
	type (group_stream_access), pointer, dimension(:) :: p, newp
	integer, intent(in) :: n
	integer :: nold, ierr
	allocate(newp(1:n), stat=ierr)
	if(ierr /= 0) stop "allocate error"
	if(.not. associated(p)) return
	nold = min(size(p), n)
	newp(1:nold) = p(1:nold)
	deallocate(p)
	end function reallocateGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Traverse group tree along some given path down to some dataset.
!> \param this Path refers to this group
!> \param depth gives level of the group containing the dataset as seen from \a this
!> \param path an index array of indices specifiying the path through the group tree, 
!!  the last index identifies the data set in the group.
!<
!> \param group output pointer to group the index array points to
!> \param dset a pointer to the dataset specified by the path
!> \par
!!  If \a group does not exist, return a null pointer.
!!  If \a dset does not exist, return a null pointer.
!!  The group pointer is returned to allow fast access to other datasets in the group.
!<
	recursive subroutine traversePathStreamAccess(this,depth,path,group,dset)
	type (group_stream_access), target :: this
	integer :: depth
	integer, dimension(:) :: path
	integer, dimension(:), allocatable :: path_copy
	type (group_stream_access), pointer :: group
	type (data_stream_access), pointer :: dset
	integer :: j,depth2
!
	allocate(path_copy(size(path)))
	path_copy = path
	if (depth > 0 ) then
		if (associated(this%subgroup)) then
			j = path_copy(1)                         ! take index from path
			path_copy = cshift(path_copy,1)          ! put first index at end
			depth2 = depth-1                         ! decrease depth
			call traversePathStreamAccess(this%subgroup(j),depth2,path_copy,group,dset)
		else
			group => null()       ! depth can't be reached
			dset => null()
		endif
	else
		j = path_copy(1)
		group => this
		if (associated(this%dataset) .and. j > 0) then
			if (size(this%dataset) >= j) then
				dset => this%dataset(j)
			else
				dset => null()
			endif
		else
			dset => null()
		endif
	endif
	deallocate(path_copy)
	end subroutine traversePathStreamAccess
!----------------------------------------------------------------------
!> \brief Create a dataset
!> \param this a data_stream_access object
!> \param rank rank of data array
!> \param dims array of integers specifiying size in each rank
!> \param datatype type of data in array using conventions in \link primitiveTypeEncoding.f90 primitiveTypeEncoding \endlink
!> \par
!> Specify rank, dims and datatype but nothing else
!! which is done later on writing to file
!<
	subroutine createDatasetStreamAccess(this,rank,dims,datatype)
	type (data_stream_access) :: this
	integer :: rank
	integer, dimension(:) :: dims
	integer :: datatype
!
	this%rank = rank
	this%datatype = datatype
	allocate(this%dims(rank))
	this%dims = dims
	if (verboseStreamAccess) then
		print *,'create dataset of size ',this%dims
	endif
!
	end subroutine createDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Fully deallocate a dataset
!> \param this data_stream_access object
!
	subroutine fullyDeallocDatasetStreamAccess(this)
	type (data_stream_access) :: this
	if (associated(this%dims)) deallocate(this%dims)
	end subroutine fullyDeallocDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Nullify the pointers contained in a dataset object
!> \param this data_stream_access object
!
	subroutine nullifyDatasetStreamAccess(this)
	type (data_stream_access) :: this
	if (associated(this%dims)) nullify(this%dims)
	end subroutine nullifyDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a dataset pointer array
!> \param p pointer to an array of datasets
!> \param n number of dataset objects to be allocated
!> \return newp pointer to the new array of datasets
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
	function reallocateDatasetStreamAccess(p, n) result(newp)
	type (data_stream_access), pointer, dimension(:) :: p, newp
	integer, intent(in) :: n
	integer :: nold, ierr
	allocate(newp(1:n), stat=ierr)
	if(ierr /= 0) stop "allocate error"
	if(.not. associated(p)) return
	nold = min(size(p), n)
	newp(1:nold) = p(1:nold)
	deallocate(p)
	end function reallocateDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief  Write an integer vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d integer data array
!
	subroutine writeDatasetIntegerVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	integer, dimension(:) :: d
	integer :: lu
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1+size(d))*kind(1))
	end subroutine writeDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a real vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d real data array
!
	subroutine writeDatasetRealVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	real, dimension(:) :: d
	integer :: lu
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.0))
	end subroutine writeDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data array
!
	subroutine writeDatasetDoubleVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	double precision, dimension(:) :: d
	integer :: lu
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.d0))
	end subroutine writeDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d complex data array
!
	subroutine writeDatasetComplexVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	complex, dimension(:) :: d
	integer :: lu
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.0))
	end subroutine writeDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double complex data array
!
	subroutine writeDatasetDoubleComplexVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	double complex, dimension(:) :: d
	integer :: lu
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.d0))
	end subroutine writeDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a character vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d character data array
!
	subroutine writeDatasetCharVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	character(len=*), dimension(:) :: d
	integer :: clen
	integer :: lu
	integer (longint) :: filepos
!
	clen = len(d(1))
	if (clen > 80) then
		print *,'writeDatasetCharVectorStreamAccess: string length greater 80 !'
		stop
	endif
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,this%datatype,clen,d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1)+size(d)*clen)
	end subroutine writeDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a flexible type vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d flexible data array
!
	subroutine writeDatasetFlexibleVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	type (flexible), dimension(:) :: d
	integer :: lu,nbytes,j
	integer (longint) :: filepos
!
	filepos = getCurrentFilePositionStreamAccess(fda)
	lu = getFileUnitStreamAccess(fda)
	this%filepos_start = filepos
	write(lu,pos=filepos) this%rank,this%dims,T_FLEXIBLE
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	do j=1,size(d)
		call writeSAFlexibleType(d(j),lu,filepos,nbytes)
		filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
	enddo
	end subroutine writeDatasetFlexibleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read an integer vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d integer data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetIntegerVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	integer, dimension(:), pointer :: d
	integer :: lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1))
	end subroutine readDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a real vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d real data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetRealVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	real, dimension(:), pointer :: d
	integer :: lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.0))
	end subroutine readDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetDoubleVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	double precision, dimension(:), pointer :: d
	integer :: lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.d0))
	end subroutine readDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d complex data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetComplexVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	complex, dimension(:), pointer :: d
	integer :: lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.0))
	end subroutine readDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double complex data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetDoubleComplexVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	double complex, dimension(:), pointer :: d
	integer :: lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) d
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.d0))
	end subroutine readDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a character vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d character data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetCharVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	character (len=80), dimension(:), pointer :: d
	integer :: clen,i,lu
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype,clen
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1))
	allocate(d(this%dims(1)))
	read(lu,pos=filepos) (d(i)(1:clen),i=1,size(d))
	filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*clen)
	end subroutine readDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a flexible type vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d flexible data array pointer (output, allocated in the subroutine)
!
	subroutine readDatasetFlexibleVectorStreamAccess(this,fda,d)
	type (data_stream_access) :: this
	type (file_stream_access), target :: fda
	type (flexible), dimension(:), pointer :: d
	integer :: lu,nbytes,j
	integer (longint) :: filepos
!
	lu = getFileUnitStreamAccess(fda)
	filepos = this%filepos_start
	call setCurrentFilePositionStreamAccess(fda,filepos)
	allocate(this%dims(1))
	read(lu,pos=filepos) this%rank,this%dims,this%datatype
	filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
	allocate(d(this%dims(1)))
	do j=1,size(d)
		call readSAFlexibleType(d(j),lu,filepos,nbytes)
		filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
	enddo
	end subroutine readDatasetFlexibleVectorStreamAccess
!
 end module
