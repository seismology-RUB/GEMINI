!--------------------------------------------------------
!  module to describe a SFF data block
!--------------------------------------------------------
 module sffDatablock
	use sffFree
	use sffInfo
	use timeSeries
	implicit none
	private readSFFDatablock, destroySFFDatablock, traceSFFDatablock
	interface new
		module procedure readSFFDatablock
		module procedure basicSFFDatablock
	end interface
	interface dealloc
		module procedure destroySFFDatablock
	end interface dealloc
	interface trace
		module procedure traceSFFDatablock
	end interface trace
	interface operator (.nsamp.); module procedure getNSFFDatablock; end interface
	interface operator (.dt.); module procedure getDtSFFDatablock; end interface
	interface operator (.station.); module procedure getStationSFFDatablock; end interface
	interface operator (.channel.); module procedure getChannelSFFDatablock; end interface
	interface operator (.tanf.); module procedure getTanfSFFDatablock; end interface
	interface operator (.year.); module procedure getYearSFFDatablock; end interface
	interface operator (.doy.); module procedure getDoySFFDatablock; end interface
	type sff_datablock
		private
		integer :: nchar;                    ! Number of charcaters in ascii code 
		real :: ampfac;                      ! Scaling factor from float to integer 
		character :: code*10                 ! Code for following optional blocks 
		character :: wid*132                 ! Identification string 
		integer :: n                         ! number of data 
		real :: dt                           ! sampling interval 
		real :: tstart                       ! x-value of first data point 
		real, dimension(:), pointer :: y     ! pointer to data 
		type (sff_free) :: free              ! free block 
		type (sff_info) :: info              ! info line 
	end type sff_datablock
!
 contains
!--------------------------------------------------------------
!  constructor
!
	subroutine readSFFDatablock(this,lu,last)
	type (sff_datablock) :: this
	logical :: last
	integer :: lu,ierr,i
	integer, dimension(:), allocatable :: idata
	character (len=1), dimension(:), allocatable :: cbuf
!
	call sff_RWid2(lu,this%wid,this%n,this%tstart,this%dt,this%nchar,this%ampfac,this%code,ierr)
	if(ierr /= 0) then
		print *,'Error while reading WID2 line: '
		stop
	endif
	allocate(idata(this%n), cbuf(this%nchar), this%y(this%n))
	call sff_RData2(lu,this%n,idata,cbuf,this%nchar,ierr)
	if(ierr /= 0) then
		print *,'Error while reading data block: ', this%wid(1:68)
		stop
	endif
	call sff_i2f(idata, this%y, this%n, this%ampfac)
	deallocate(cbuf, idata)
!
!  read optional free block and info line
!
	i=1
	last=.true.
	do, while (this%code(i:i) /= ' ')
		if (this%code(i:i).eq.'F') then
			call new(this%free,lu)
		endif
		if (this%code(i:i).eq.'I') then
			call new(this%info,lu)
		endif
		if (this%code(i:i).eq.'D') last=.false.
		i=i+1
	enddo
	if(last) close(lu)
	end subroutine readSFFDatablock
!-----------------------------------------------------------------------------
!  create basic SFF datablock
!
	subroutine basicSFFDatablock(this,n,dt,y,net,station,comp,last)
	type (sff_datablock) :: this
	integer :: n,ierr
	logical last
	real :: dt
	real, dimension(:), target :: y
	character (len=*) :: net,station,comp
!
	this%n = n
	this%dt = dt
	this%y => y
!
	call sff_PrepWid2(n, 1./dt, station, -1,-1,-1,-1,-1,comp,net,'NSP',-1.,-1.,-1.,-1.,-1.,this%wid,ierr)
	this%code='D'
	if(last) this%code = ' '
	end subroutine basicSFFDatablock
!-----------------------------------------------------------------------------
!  destructor
!
	subroutine destroySFFDatablock(this)
	type (sff_datablock) :: this
	integer :: i
!
	i=1
	do, while (this%code(i:i) /= ' ')
		if (this%code(i:i).eq.'F') call dealloc(this%free)
		i=i+1
	enddo
	if(associated(this%y)) deallocate(this%y)
	end subroutine destroySFFDatablock
!------------------------------------------------------------------------------
!  return a pointer to the data
!
	function traceSFFDatablock(this)
	type (sff_datablock), intent(in) :: this
	real, dimension(:), pointer :: traceSFFDatablock
!
	traceSFFDatablock => this%y
	end function traceSFFDatablock
!------------------------------------------------------------------------------
!  clear the data array
!
	subroutine clearDataSFFDatablock(this)
	type (sff_datablock) :: this
!
	if(associated(this%y)) deallocate(this%y)
	end subroutine clearDataSFFDatablock
!------------------------------------------------------------------------------
!  nullify the data array pointer
!
	subroutine nullifyDataSFFDatablock(this)
	type (sff_datablock) :: this
!
	if(associated(this%y)) nullify(this%y)
	end subroutine nullifyDataSFFDatablock
!------------------------------------------------------------------------------
!  get dt of datablock
!
	function getDtSFFDatablock(this)
	type (sff_datablock), intent(in) :: this
	real :: getDtSFFDatablock
	getDtSFFDataBlock = this%dt
	end function getDtSFFDatablock
!------------------------------------------------------------------------------
!  get station name of datablock
!
	character (len=5) function getStationSFFDatablock(this)
	type (sff_datablock), intent(in) :: this
	character (len=5) :: station 
	read(this%wid(30:34),'(a5)') station
	getStationSFFDatablock = station
	end function getStationSFFDatablock
!------------------------------------------------------------------------------
!  get channel of datablock
!
	character (len=3) function getChannelSFFDatablock(this)
	type (sff_datablock), intent(in) :: this
	read(this%wid(36:38),'(a3)') getChannelSFFDatablock
	end function getChannelSFFDatablock
!------------------------------------------------------------------------------
!  get start time of data block
!
	real function getTanfSFFDataBlock(this)
	type (sff_datablock), intent(in) :: this
	getTanfSFFDataBlock = this%tstart
	end function getTanfSFFDataBlock
!------------------------------------------------------------------------------
!  get year of data block
!
	integer function getYearSFFDataBlock(this)
	type (sff_datablock), intent(in) :: this
	read(this%wid(6:9),'(i4)') getYearSFFDataBlock
	end function getYearSFFDataBlock
!------------------------------------------------------------------------------
!  get day of year of data block
!
	integer function getDoySFFDataBlock(this)
	type (sff_datablock), intent(in) :: this
	integer year,month,day,sff_TimeGetDOY
	read(this%wid(6:9),'(i4)') year
	read(this%wid(11:12),'(i2)') month
	read(this%wid(14:15),'(i2)') day
	getDoySFFDataBlock = sff_TimeGetDOY(year,month,day)
	end function getDoySFFDataBlock
!------------------------------------------------------------------------------
!  get nsamp of datablock 
!
	integer function getNSFFDatablock(this)
	type (sff_datablock), intent(in) :: this
	getNSFFDatablock = this%n
	end function getNSFFDatablock
!------------------------------------------------------------------------------
!  print wid2 line to unit
!
	subroutine printWid2SFFDatablock(this,lu)
	type (sff_datablock) :: this
	integer :: lu
	write(lu,'(a)') this%wid(1:80)
	end subroutine printWid2SFFDatablock
!------------------------------------------------------------------------------
!  get receiver location
!  return (1,0,0) if info line is not available
!
	subroutine getReceiverLocationSFFDatablock(this,xe,ye,ze)
	type (sff_datablock) :: this
	real :: xe,ye,ze
	if (index(this%code,'I') > 0) then
		call locationSFFInfo(this%info,xe,ye,ze)
	else
		xe = 1.; ye = 0.; ze = 0.
	endif
	end subroutine getReceiverLocationSFFDatablock
!------------------------------------------------------------------------------
!  modify tstart of datablock, also change wid2line 
!
	subroutine modifyTanfSFFDatablock(this,tstart)
	type (sff_datablock) :: this
	real :: tstart,second
	integer :: day,hour,minute
!
	this%tstart = tstart
	call sff_TimeSplit(tstart, day, hour, minute, second)
	call sff_ModWid2time(this%wid, hour, minute, second)
	end subroutine modifyTanfSFFDatablock
!------------------------------------------------------------------------------
!  modify dt of datablock, also change wid2line 
!
	subroutine modifyDtSFFDatablock(this,dt)
	type (sff_datablock) :: this
	real :: dt
!
	this%dt = dt
	call sff_ModWid2samprat(this%wid,1./dt)
	end subroutine modifyDtSFFDatablock
!------------------------------------------------------------------------------
!  modify n of datablock, also change wid2line 
!
	subroutine modifyNSFFDatablock(this,n)
	type (sff_datablock) :: this
	integer :: n
!
	this%n = n
	call sff_ModWid2samps(this%wid,n)
	end subroutine modifyNSFFDatablock
!------------------------------------------------------------------------------
!  modify date (year month day) in wid2line 
!
	subroutine modifyDateSFFDatablock(this,year,month,day)
	type (sff_datablock) :: this
	integer :: year,month,day
!
	call sff_ModWid2date(this%wid, year, month, day)
	end subroutine modifyDateSFFDatablock
!------------------------------------------------------------------------------
!  add info line to datablock, update code
!
	subroutine addInfoSFFDatablock(this,cs,x,y,z,nstack)
	type (sff_datablock) :: this
	character cs*1
	real :: x,y,z
	integer :: nstack
!
	if(cs == 'S') call new(this%info,cs,90.-x,y,z,nstack)
	if(cs == 'C') call new(this%info,cs,x,y,z,nstack)
	this%code = trim(this%code)//'I'
	end subroutine addInfoSFFDatablock
!------------------------------------------------------------------------------
!  write a data block to file id lu as is; file is assumed to be open
!  close file if last = .true.
!  if d is present use d as data instead of this%y
!
	subroutine writeSFFDataBlock(this,lu,last,d)
	type (sff_datablock) :: this
	integer :: lu,j
	logical :: last
	real :: ampfac
	integer, dimension(:), allocatable :: idata
	integer :: lc
	character (len=10) :: code
	real, dimension(:), optional, target :: d
	real, dimension(:), pointer :: temp
!
!  if d is present, redirect data pointer y to d
!  but save y before redirection
!
	if(present(d)) then
		if(associated(this%y)) temp => this%y
		this%y => d
	endif
!
!  use a copy of this%code, but
!  if D is in the code of this datablock, replace it by a blank if last=.true.
!  if D is not in the code of this datablock, add a D if last=.false.
!
	code = this%code
	lc = len_trim(code)+1
	j = index(code,'D')
	if(j > 0 .and. last) code(j:j) = ' '
	if(j == 0 .and. .not.last) code(lc:lc) = 'D'
!
	allocate(idata(this%n))
	call sff_f2i(idata, this%y, this%n, ampfac)
	call sff_WData(lu, this%wid, this%n, idata, ampfac, code)
!
!  write free block and info line if present
! 
	j=1
	do, while(code(j:j) /= ' ')
		if(code(j:j) == 'F') then
			call writeSFFFree(this%free,lu)
		endif
		if(code(j:j) == 'I') then
			call writeSFFInfo(this%info,lu)
		endif
		j=j+1
	enddo
	if(last) close(lu)
	deallocate(idata)
!
!  if d is present, redirect this%y to its former target
!
	if(present(d)) then
		if(associated(temp)) this%y => temp
	endif
	end subroutine writeSFFDataBlock
!------------------------------------------------------------------------------
!  reallocate SFF datablock type
!
	function reallocateSFFDataBlock(p, n)    ! reallocate sff_datablock
	type (sff_datablock), pointer, dimension(:) :: p, reallocateSFFDataBlock
	integer, intent(in) :: n
	integer :: nold, ierr
	allocate(reallocateSFFDataBlock(1:n), stat=ierr)
	if(ierr /= 0) stop "SFFDataBlock reallocate error"
	if(.not. associated(p)) return
	nold = min(size(p), n)
	reallocateSFFDataBlock(1:nold) = p(1:nold)
	deallocate(p)
	end function reallocateSFFDataBlock
!------------------------------------------------------------------------------
!
 end module sffDatablock
