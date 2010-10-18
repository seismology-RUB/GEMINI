!-------------------------------------------------------------
!> \brief Routines that operate with Fourier spectra
!-------------------------------------------------------------
 module fourierSpectrum
	use fourierTransform
	use timeSeries
	implicit none
	interface dealloc; module procedure deallocFourierSpectrum; end interface
	interface operator (.fmin.); module procedure fminFourierSpectrum; end interface
	interface operator (.fmax.); module procedure fmaxFourierSpectrum; end interface
	interface operator (.df.); module procedure dfFourierSpectrum; end interface
	interface operator (.nf.); module procedure nfFourierSpectrum; end interface
	interface operator (.amp.); module procedure amplitudeFourierSpectrum; end interface
	interface operator (.phase.); module procedure phaseFourierSpectrum; end interface
	type fourier_spectrum
		private
		integer :: nf1,nf2,nf                        !< Index of first and last frequency and number of frequencies
 		real :: df                                   !< Frequency interval
		logical :: link                              !< If true data is just a pointer
		double complex, dimension(:), pointer :: y => null()   !< Pointer to data
	end type
!
 contains
!----------------------------------------------------------------
!> \brief Create spectrum with a link to the data
!
	subroutine createLinkFourierSpectrum(this,nf1,nf2,df,s)
	type (fourier_spectrum) :: this
	integer :: nf1,nf2
	real :: df
	double complex, dimension(:), target :: s
!
	this%nf1 = nf1; this%nf2 = nf2; this%df = df; this%nf = nf2-nf1+1
	this%y => s
	this%link = .true.
	end subroutine createLinkFourierSpectrum
!------------------------------------------------------------------
!> \brief Create spectrum with a true copy of the data
!
	subroutine createFromDataFourierSpectrum(this,nf1,nf2,df,s)
	type (fourier_spectrum) :: this
	integer :: nf1,nf2
	real :: df
	double complex, dimension(:) :: s
!
	this%nf = nf2-nf1+1
	this%nf1 = nf1; this%nf2 = nf2; this%df = df
	allocate(this%y(this%nf))
	this%y(1:this%nf) = s(1:this%nf)
	this%link = .true.
	end subroutine createFromDataFourierSpectrum
!-----------------------------------------------------------------
!> \brief Compute Fourier spectrum of a time series
!
	subroutine createFromTimeSeriesFourierSpectrum(this,ts)
	type (fourier_spectrum) :: this
	type (time_series) :: ts
	integer :: nplog,np2,nsamp
	double precision, dimension(:), allocatable :: x,y
!
! extend number of samples to next power of 2
!
	nsamp = .nsamp.ts
	nplog = ceiling(log(real(nsamp))/log(2.))
	np2 = 2**nplog
	this%df = 1./(np2*.dt.ts)
!
! copy data to double precision and append zeros
!
	allocate(x(np2),y(np2))
	x(1:nsamp) = .trace.ts
	x(nsamp+1:np2) = 0.d0
	y = 0.d0
!
! Fourier transform
! 
	call fastFourierTransform(x,y,nplog,-1)
!
! fill fourierTransform Structure
!
	this%nf1 = 1
	this%nf2 = np2/2
	this%nf = np2/2
	this%df = 1./(np2*.dt.ts)
	this%link = .false.
	allocate(this%y(np2/2))
	this%y(1:np2/2) = dcmplx(x(1:np2/2),y(1:np2/2))
!
	deallocate(x,y)
	end subroutine createFromTimeSeriesFourierSpectrum
!------------------------------------------------------------
!> \brief Deallocate fourier spectrum
!
	subroutine deallocFourierSpectrum(this)
	type (fourier_spectrum) :: this
	if (associated(this%y)) then
		if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
	endif
	end subroutine deallocFourierSpectrum
!------------------------------------------------------------
!> \brief Amplitude spectrum
!
	function amplitudeFourierSpectrum(this) result(amp)
	type (fourier_spectrum), intent(in) :: this
	real, dimension(:), pointer :: amp
!
	allocate(amp(this%nf))
	amp = sqrt(zabs(this%y))
	end function amplitudeFourierSpectrum
!------------------------------------------------------------
!> \brief Phase spectrum
!
	function phaseFourierSpectrum(this) result(phase)
	type (fourier_spectrum), intent(in) :: this
	real, dimension(:), pointer :: phase
!
	allocate(phase(this%nf))
	phase = real(datan2(dimag(this%y),dreal(this%y)))
	end function phaseFourierSpectrum
!------------------------------------------------------------
!> \brief Get df
!
	real function dfFourierSpectrum(this)
	type (fourier_spectrum), intent(in) :: this
	dfFourierSpectrum = this%df
	end function dfFourierSpectrum
!------------------------------------------------------------
!> \brief Get nf
!
	integer function nfFourierSpectrum(this)
	type (fourier_spectrum), intent(in) :: this
	nfFourierSpectrum = this%nf
	end function nfFourierSpectrum
!------------------------------------------------------------
!> \brief Get fmin
!
	real function fminFourierSpectrum(this)
	type (fourier_spectrum), intent(in) :: this
	fminFourierSpectrum = (this%nf1-1)*this%df
	end function fminFourierSpectrum
!------------------------------------------------------------
!> \brief Get fmax
!
	real function fmaxFourierSpectrum(this)
	type (fourier_spectrum), intent(in) :: this
	fmaxFourierSpectrum = (this%nf2-1)*this%df
	end function fmaxFourierSpectrum
!
 end module
