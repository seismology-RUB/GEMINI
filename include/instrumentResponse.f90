!--------------------------------------------------------------------
!> \brief Instrument response module
!!
!!  station:		name of station
!!  net:		network name
!!  channel:		channel name
!!  location:		location id
!!  startyear:	year when operation began
!!  startdoy:		day of year when operation began
!!  endyear:		year when operation stopped
!!  enddoy:		day of year when operation stopped
!!  nzeros:		number of zeros in Laplace domain
!!  npoles:		number of poles in Laplace domain
!!  sz(i):		complex zeros in Laplace s-plane
!!  sp(i):		complex poles in Laplace s-plane
!!  gain_inst:	instrument gain in V/physical unit (e.g. m/s, Pa)
!!  gain_adc:		gain of AD converter in counts/V
!!  sensitivity:	total sensitivity of instrument in counts/physical unit
!!  fs:		frequency at which sensitivity was measured
!<--------------------------------------------------------------------
 module instrumentResponse
	use mathConstants
	implicit none
	private :: normalization
	interface new
		module procedure readInstrumentResponse
		module procedure createInstrumentResponse
	end interface
	interface operator (.sensitivity.); module procedure getSensitivityInstrumentResponse; end interface
	interface operator (.eigenperiod.); module procedure getEigenperiodInstrumentResponse; end interface
	interface operator (.damping.); module procedure getDampingConstantInstrumentResponse; end interface
	interface dealloc
		module procedure deallocInstrumentResponse
	end interface
	type instrument_response
		private
		character :: station*5, net*6, channel*3, location*2
		integer :: nzeros, npoles, startyear, startdoy, endyear, enddoy
		complex, dimension(:), pointer :: sz => null()
		complex, dimension(:), pointer :: sp => null()
		real :: gain_inst, gain_adc, sensitivity, fs
	end type instrument_response
!
	contains
!-------------------------------------------------------------------------------
	subroutine readInstrumentResponse(this,lu,filename)
	type (instrument_response) :: this
	integer :: lu,i
	character (len=*) :: filename
	real :: a,b
!
	open(lu,file=filename, status='old')
	read(lu,'(1x/1x/1x)')
	read(lu,'(a)') this%station
	read(lu,'(a)') this%net
	read(lu,'(a)') this%channel
	read(lu,'(a)') this%location
	read(lu,*) this%startyear
	read(lu,*) this%startdoy
	read(lu,*) this%endyear
	read(lu,*) this%enddoy
	read(lu,*) this%nzeros
	read(lu,*) this%npoles
	allocate(this%sz(this%nzeros),this%sp(this%npoles))
	do i=1,this%nzeros
		read(lu,*) a,b
		this%sz(i) = cmplx(a,b)
	enddo
	do i=1,this%npoles
		read(lu,*) a,b
		this%sp(i) = cmplx(a,b)
	enddo
	read(lu,*) this%gain_inst
	read(lu,*) this%gain_adc
	read(lu,*) this%sensitivity
	read(lu,*) this%fs
	close(lu)
	end subroutine readInstrumentResponse
!-------------------------------------------------------------------------------
!  create instrument response object using arguments
!
	subroutine createInstrumentResponse(this,sta,net,channel,location,ya,da,ye,de,nz,np,sz,sp,fref,pg,gc,lsb)
	type (instrument_response) :: this
	character (len=*) :: sta,net,channel,location
	integer :: ya,ye,da,de,nz,np
	complex, dimension(:) :: sz,sp
	real :: fref,pg,gc,lsb
	this%station = sta
	this%net = net
	this%channel = channel
	this%location = location
	this%startyear = ya
	this%startdoy = da
	this%endyear = ye
	this%enddoy = de
	this%nzeros = nz
	this%npoles = np
	allocate(this%sz(nz),this%sp(np))
	this%sz(1:nz) = sz(1:nz)
	this%sp(1:np) = sp(1:np)
	this%fs = fref
	this%gain_inst = gc
	this%gain_adc = pg/lsb
	this%sensitivity = this%gain_inst*this%gain_adc
	end subroutine createInstrumentResponse
!-------------------------------------------------------------------------------
!  deep copy of instrument response
!
	subroutine deepCopyInstrumentResponse(this,copy)
	type (instrument_response) :: this,copy
	copy%station = this%station
	copy%net = this%net
	copy%channel = this%channel
	copy%location = this%location
	copy%startyear = this%startyear
	copy%startdoy = this%startdoy
	copy%endyear = this%endyear
	copy%enddoy = this%enddoy
	copy%nzeros = this%nzeros
	copy%npoles = this%npoles
	allocate(copy%sz(this%nzeros),copy%sp(this%npoles))
	copy%sz = this%sz
	copy%sp = this%sp
	copy%fs = this%fs
	copy%gain_inst = this%gain_inst
	copy%gain_adc = this%gain_adc
	copy%sensitivity = this%sensitivity
	end subroutine deepCopyInstrumentResponse
!-------------------------------------------------------------------------------
	subroutine deallocInstrumentResponse(this)
	type (instrument_response) :: this
	if(associated(this%sz)) deallocate(this%sz)
	if(associated(this%sp)) deallocate(this%sp)
	end subroutine deallocInstrumentResponse
!-------------------------------------------------------------------------------
	subroutine unlinkDataInstrumentResponse(this)
	type (instrument_response) :: this
	if(associated(this%sz)) nullify(this%sz)
	if(associated(this%sp)) nullify(this%sp)
	end subroutine unlinkDataInstrumentResponse
!------------------------------------------------------------------------------
!  abs value of response at frequency of sensitivity
!
	real function normalization(this)
	type (instrument_response) :: this
	complex :: resp
	integer :: i
	real :: om
!
	resp = (1.,0.)
	om = 2.*pi*this%fs
	do i=1,this%nzeros
		resp=resp*(mc_ci*om-this%sz(i))
	enddo
	do i=1,this%npoles
		resp=resp/(mc_ci*om-this%sp(i))
	enddo
	normalization = cabs(resp)
	end function normalization
!------------------------------------------------------------------------------
!  calculate instrument transfer function
!
	subroutine getInstrumentResponse(this,year,doy,nf1,nf2,df,resp)
	type (instrument_response) :: this
	integer :: nf1,nf2,year,doy,j,i,jf
	real :: df,om
	complex, dimension(:) :: resp
!
	if (.not. checkTimeInstrumentResponse(this,year,doy)) then
		print *,'Response is invalid for year and doy of data block'
		print *,'Year: ',year,' start year: ',this%startyear,' end year: ',this%endyear
		print *,'Doy: ',doy,' start doy: ',this%startdoy,' end doy: ',this%enddoy
		stop
	endif
!
	do j = nf1,nf2
		om = 2.*pi*(j-1)*df
		jf = j-nf1+1
		resp(jf) = (1.,0.)
		do i=1,this%nzeros
			resp(jf)=resp(jf)*(mc_ci*om-this%sz(i))
		enddo
		do i=1,this%npoles
			resp(jf)=resp(jf)/(mc_ci*om-this%sp(i))
		enddo
		resp(jf) = resp(jf)/normalization(this)*this%sensitivity
	enddo
	end subroutine getInstrumentResponse
!------------------------------------------------------------------------------
!  calculate instrument response at one selected frequency
!
	function getSingleFrequencyInstrumentResponse(this,year,doy,f) result(resp)
	type (instrument_response) :: this
	integer :: year,doy,i
	real :: om,f
	complex :: resp
!
	if (.not. checkTimeInstrumentResponse(this,year,doy)) then
		print *,'Response is invalid for year and doy of data block'
		print *,'Year: ',year,' start year: ',this%startyear,' end year: ',this%endyear
		print *,'Doy: ',doy,' start doy: ',this%startdoy,' end doy: ',this%enddoy
		stop
	endif
!
	om = 2.*pi*f
	resp = (1.,0.)
	do i=1,this%nzeros
		resp = resp*(mc_ci*om-this%sz(i))
	enddo
	do i=1,this%npoles
		resp = resp/(mc_ci*om-this%sp(i))
	enddo
	resp = resp/normalization(this)*this%sensitivity
	end function getSingleFrequencyInstrumentResponse
!---------------------------------------------------------------------------
!  check validity of response for given year and doy
!
	logical function checkTimeInstrumentResponse(this,year,doy)
	type (instrument_response) :: this
	integer :: year,doy
!
	checkTimeInstrumentResponse = .true.
	if(year < this%startyear .or. year > this%endyear) checkTimeInstrumentResponse = .false.
	if(year == this%startyear) then
		if(doy < this%startdoy) checkTimeInstrumentResponse = .false.
	endif
	if(year == this%endyear) then
		if(doy > this%enddoy) checkTimeInstrumentResponse = .false.
	endif
	end function checkTimeInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Get sensitivity
!
	function getSensitivityInstrumentResponse(this) result(res)
	type (instrument_response), intent(in) :: this
	real :: res
	res = this%sensitivity
	end function getSensitivityInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Get eigenperiod
!
	function getEigenperiodInstrumentResponse(this) result(res)
	type (instrument_response), intent(in) :: this
	real :: res
	res = -1.
	select case (this%npoles)
	case (2); res = 2.*mc_pi/sqrt(abs(this%sp(1)*this%sp(2)))
	case (1); res = 2.*mc_pi/abs(this%sp(1))
	end select
	end function getEigenperiodInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Get damping constant
!
	function getDampingConstantInstrumentResponse(this) result(res)
	type (instrument_response), intent(in) :: this
	real :: res
	res = -1.
	select case (this%npoles)
	case (2); res = 0.5*abs(this%sp(1)+this%sp(2))/sqrt(abs(this%sp(1)*this%sp(2)))
	end select
	end function getDampingConstantInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Extend a pointer array of instrument_response objects
!
	function extendArrayInstrumentResponse(array,n) result(newarray)
	type (instrument_response), dimension(:), pointer :: array
	type (instrument_response), dimension(:), pointer :: newarray
	integer :: n,nold,i
!
	allocate(newarray(n))
	if (.not. associated(array)) return
	nold = min(size(array),n)
	newarray(1:nold) = array(1:nold)
	do i = 1,nold
		call unlinkDataInstrumentResponse(array(i))
	enddo
	do i = nold+1,size(array)
		call dealloc(array(i))
	enddo
	deallocate(array)
	end function extendArrayInstrumentResponse
!-------------------------------------------------------------------------------
  end module instrumentResponse
