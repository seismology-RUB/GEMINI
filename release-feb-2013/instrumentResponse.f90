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
!> \brief Instrument response module
!!
!!  station:		name of station
!!  channel:		channel name
!!  nzeros:		number of zeros in Laplace domain
!!  npoles:		number of poles in Laplace domain
!!  sz(i):		complex zeros in Laplace s-plane
!!  sp(i):		complex poles in Laplace s-plane
!!  sensitivity:	total sensitivity of instrument in counts/physical unit
!!  fs:			frequency at which sensitivity was measured
!!
!!  Assume here that response refers to ground velocity
!!  Let R(om) be the Fourier transform of the recorded signal
!!  Let VEL(om) be the Fourier transform of the ground velocity
!!  Let HV(om) the velocity transfer function of the instrument 
!!  as calculated from poles and zeros:
!!  HV(om) = product_k (ci*om-sz(k))/product_k (ci*om-sp(k)) 
!!  then: R(om) = S*HV(om)/|HV(omref)|*VEL(om)
!!  where S is the sensitivity constant
!!
!!  Note in seismicTrace we use:
!!  R(om) = C*H(om)*V(om) and hence S = C*|H(omref)|
!<--------------------------------------------------------------------
 module instrumentResponse
	use mathConstants
	use flexibleType
	use errorMessage
	implicit none
	interface new
		module procedure readInstrumentResponse
		module procedure createInstrumentResponse
	end interface
	interface operator (.sensitivity.); module procedure getSensitivityInstrumentResponse; end interface
	interface operator (.eigenperiod.); module procedure getEigenperiodInstrumentResponse; end interface
	interface operator (.damping.); module procedure getDampingConstantInstrumentResponse; end interface
	interface operator (.nz.); module procedure getNzerosInstrumentResponse; end interface
	interface operator (.np.); module procedure getNpolesInstrumentResponse; end interface
	interface dealloc
		module procedure deallocInstrumentResponse
	end interface
	type instrument_response
		private
		character :: station*5, channel*3
		integer :: nzeros, npoles
		complex, dimension(:), pointer :: sz => null()
		complex, dimension(:), pointer :: sp => null()
		real :: sensitivity, fs
		logical :: link
	end type instrument_response
!
	contains
!-------------------------------------------------------------------------------
!> \brief  read response file
!
	subroutine readInstrumentResponse(this,lu,filename)
	type (instrument_response) :: this
	integer :: lu,i
	character (len=*) :: filename
	real :: a,b
!
	open(lu,file=filename, status='old')
	read(lu,'(1x/1x/1x)')
	read(lu,'(a)') this%station
	read(lu,'(1x)')
	read(lu,'(a)') this%channel
	read(lu,'(1x/1x/1x/1x/1x)')
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
	read(lu,'(1x/1x)')
	read(lu,*) this%sensitivity
	read(lu,*) this%fs
	close(lu)
	this%link = .false.
	end subroutine readInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief create instrument response object using arguments
!
	subroutine createInstrumentResponse(this,sta,channel,nz,np,sz,sp,fref,sens)
	type (instrument_response) :: this
	character (len=*) :: sta,channel
	integer :: nz,np
	complex, dimension(:) :: sz,sp
	real :: fref,sens
	this%station = sta
	this%channel = channel
	this%nzeros = nz
	this%npoles = np
	allocate(this%sz(nz),this%sp(np))
	this%sz(1:nz) = sz(1:nz)
	this%sp(1:np) = sp(1:np)
	this%fs = fref
	this%sensitivity = sens
	this%link = .false.
	end subroutine createInstrumentResponse
!--------------------------------------------------------------------------------
!> \brief create response using CONSTANT of SAC poles and zeros files
!! Here R(om) = CONSTANT*HD(om)*DIS(om) = CONSTANT*HD/(i*om)*(i*om*DIS(om)) = CONSTANT*HD/(i*om)*VEL(om)
!! with HD(om) = product_k (ci*om-sz(k))/product_k (ci*om-sp(k))
!! and DIS(om) is displacement
!! Because R(om) = S*HV(om)/|HV(omref)|*VEL(om) = CONSTANT*HD/(i*om)*VEL(om) and HV = HD/(i*om)
!! we find: S/|HV(omref)| = CONSTANT or S = CONSTANT*|HV(omref)|
!> \param fref Reference frequency used to calculate |HV(omref)|
!
	subroutine createFromSACPolesZerosInstrumentResponse(this,sta,channel,nz,np,sz,sp,fref,constant)
	type (instrument_response) :: this
	character (len=*) :: sta,channel
	integer :: nz,np
	complex, dimension(:) :: sz,sp
	real :: constant,norm,fref
!
	this%station = sta
	this%channel = channel
	this%nzeros = nz-1                             ! skip last zero to get HV
	this%npoles = np
	allocate(this%sz(nz-1),this%sp(np))
	this%sz(1:nz-1) = sz(1:nz-1)                   ! skip last zero to get HV
	this%sp(1:np) = sp(1:np)
	this%fs = fref
	norm = normalizationInstrumentResponse(this)   ! |HV(omref)|
	this%sensitivity = constant*norm
	this%link = .false.
	end subroutine createFromSACPolesZerosInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief create from flexible array
!
	subroutine createFromFlexibleArrayInstrumentResponse(this,ft)
	type (instrument_response) :: this
	type (flexible), dimension(:) :: ft
	integer :: j
!
	this%station = ft(1); this%channel = ft(2)
	this%nzeros = ft(3); this%npoles = ft(4)
	this%fs = ft(5); this%sensitivity = ft(6)
	allocate(this%sz(this%nzeros),this%sp(this%npoles))
	do j = 1,this%nzeros
		this%sz(j) = ft(6+j)
	enddo
	do j = 1,this%npoles
		this%sp(j) = ft(6+j+this%nzeros)
	enddo
	this%link = .false.
	end subroutine createFromFlexibleArrayInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief read in SAC poles and zeros
!
	function readSACPolesZerosInstrumentResponse(this,lu,filename,fref) result(errmsg)
	type (instrument_response) :: this
	character (len=*) :: filename
	integer :: lu
	type (error_message) :: errmsg
	character (len=132) :: line
	complex, dimension(:), allocatable :: sz,sp
	integer :: nz,np,poscol,j,cnt,ios
	real :: fref,constant,spr,spi
	character :: sta*5,channel*3
	character (len=35) :: myname = 'readSACPolesZerosInstrumentResponse'
!
	call new(errmsg,myname)
	open(lu,file = filename,status = 'old', iostat = ios)
	if (ios /= 0) then
		call new(errmsg,2,'Response file: '//trim(filename)//' can not be opened',myname)
		return
	endif
	read(lu,'(1x)')
	read(lu,'(a)') line
	poscol = index(line,':')
	do while (index(line,'Station') == 0); read(lu,'(a)') line; enddo
	read(line(poscol+1:),*) sta
	do while (index(line,'Component') == 0); read(lu,'(a)') line; enddo
	read(line(poscol+1:),*) channel
!
!  zeros
!
	do while (index(line,'ZEROS') == 0); read(lu,'(a)') line; enddo
	read(line(6:),*) nz
	allocate(sz(nz))
	read(lu,'(a)') line
	cnt = 0
	do while (index(line,'POLES') == 0)
		cnt = cnt+1
		read(line,*) spr,spi
		sz(cnt) = cmplx(spr,spi)
		read(lu,'(a)') line
	enddo
	sz(cnt+1:nz) = 0.0
	if (cabs(sz(nz)) > 1.e-6) then
		call new(errmsg,2,'No zeros at zero',myname)
		return
	endif
!
!  poles
!
	read(line(6:),*) np
	allocate(sp(np))
	do j = 1,np
		read(lu,*) spr,spi
		sp(j) = cmplx(spr,spi)
	enddo
	read(lu,'(a)') line; read(line(9:),*) constant
	close(lu)
	call createFromSACPolesZerosInstrumentResponse(this,sta,channel,nz,np,sz,sp,fref,constant)
	deallocate(sz,sp)
	end function readSACPolesZerosInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief  deep copy of instrument response
!
	subroutine deepCopyInstrumentResponse(this,copy)
	type (instrument_response) :: this,copy
	copy%station = this%station
	copy%channel = this%channel
	copy%nzeros = this%nzeros
	copy%npoles = this%npoles
	allocate(copy%sz(this%nzeros),copy%sp(this%npoles))
	copy%sz = this%sz
	copy%sp = this%sp
	copy%fs = this%fs
	copy%sensitivity = this%sensitivity
	copy%link = .false.
	end subroutine deepCopyInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief  deallocate instrument response
!
	subroutine deallocInstrumentResponse(this)
	type (instrument_response) :: this
	if (.not. (this%link)) then
		if(associated(this%sz)) deallocate(this%sz)
		if(associated(this%sp)) deallocate(this%sp)
	else
		if(associated(this%sz)) nullify(this%sz)
		if(associated(this%sp)) nullify(this%sp)
	endif
	end subroutine deallocInstrumentResponse
!-------------------------------------------------------------------------------
!> \brief Unlink data instrument response
!
	subroutine unlinkDataInstrumentResponse(this)
	type (instrument_response) :: this
	if(associated(this%sz)) nullify(this%sz)
	if(associated(this%sp)) nullify(this%sp)
	end subroutine unlinkDataInstrumentResponse
!------------------------------------------------------------------------------
!> \brief  abs value of response at frequency of sensitivity
!
	function normalizationInstrumentResponse(this) result(res)
	type (instrument_response) :: this
	complex :: resp
	real :: res,om
	integer :: i
!
	resp = (1.,0.)
	om = 2.*mc_pi*this%fs
	do i=1,this%nzeros
		resp=resp*(mc_ci*om-this%sz(i))
	enddo
	do i=1,this%npoles
		resp=resp/(mc_ci*om-this%sp(i))
	enddo
	res = cabs(resp)
	end function normalizationInstrumentResponse
!------------------------------------------------------------------------------
!> \brief  calculate instrument transfer function
!
	function getInstrumentResponse(this,nf1,nf2,df,sigma) result(resp)
	type (instrument_response) :: this
	integer :: nf1,nf2,j,i,jf
	real :: df,sigma
	complex :: zom
	complex, dimension(:), pointer :: resp
!
	allocate(resp(nf2-nf1+1))
	do j = nf1,nf2
		zom = 2.*mc_pi*(j-1)*df-mc_ci*sigma
		jf = j-nf1+1
		resp(jf) = (1.,0.)
		do i=1,this%nzeros
			resp(jf)=resp(jf)*(mc_ci*zom-this%sz(i))
		enddo
		do i=1,this%npoles
			resp(jf)=resp(jf)/(mc_ci*zom-this%sp(i))
		enddo
		resp(jf) = resp(jf)/normalizationInstrumentResponse(this)*this%sensitivity
	enddo
	end function getInstrumentResponse
!------------------------------------------------------------------------------
!> \brief calculate instrument response at one selected frequency
!
	function getSingleFrequencyInstrumentResponse(this,f,sigma) result(resp)
	type (instrument_response) :: this
	integer :: i
	real :: f,sigma
	complex :: zom
	complex :: resp
!
	zom = 2.*mc_pi*f-mc_ci*sigma
	resp = (1.,0.)
	do i=1,this%nzeros
		resp = resp*(mc_ci*zom-this%sz(i))
	enddo
	do i=1,this%npoles
		resp = resp/(mc_ci*zom-this%sp(i))
	enddo
	resp = resp/normalizationInstrumentResponse(this)*this%sensitivity
	end function getSingleFrequencyInstrumentResponse
!--------------------------------------------------------------------------
!> \brief Get mebers of instrument response in a flexible array
!
	function getMembersInstrumentResponse(this) result(ft)
	type (instrument_response) :: this
	type (flexible), dimension(:), pointer :: ft
	integer :: j
!
	allocate(ft(6+this%nzeros+this%npoles))
	ft(1) = this%station; ft(2) = this%channel
	ft(3) = this%nzeros; ft(4) = this%npoles
	ft(5) = this%fs; ft(6) = this%sensitivity
	do j = 1,this%nzeros
		ft(6+j) = this%sz(j)
	enddo
	do j = 1,this%npoles
		ft(6+j+this%nzeros) = this%sp(j)
	enddo
	end function getMembersInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Get nzeros
!
	function getNzerosInstrumentResponse(this) result(res)
	type (instrument_response), intent(in) :: this
	integer :: res
	res = this%nzeros
	end function getNzerosInstrumentResponse
!---------------------------------------------------------------------------
!> \brief Get npoles
!
	function getNpolesInstrumentResponse(this) result(res)
	type (instrument_response), intent(in) :: this
	integer :: res
	res = this%npoles
	end function getNpolesInstrumentResponse
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
!> \brief shallow copy of instrument response
!
	subroutine shallowCopyInstrumentResponse(this,copy)
	type (instrument_response) :: this,copy
	copy%station = this%station
	copy%channel = this%channel
	copy%nzeros = this%nzeros
	copy%npoles = this%npoles
	copy%sz => this%sz
	copy%sp => this%sp
	copy%fs = this%fs
	copy%sensitivity = this%sensitivity
	copy%link = .true.
	end subroutine shallowCopyInstrumentResponse
!-------------------------------------------------------------------------------
  end module instrumentResponse