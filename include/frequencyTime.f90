!----------------------------------------------------------
!> \brief Routines for transformation from frequency to time domain and vice versa
!----------------------------------------------------------
 module frequencyTime
	use fourierTransform
	use errorMessage
	implicit none
!
 contains
!----------------------------------------------------------------------------
!  calculate nsamp and dt of time series after backtransformation from frequency domain
!  such that number of samples is power of 2
!
	subroutine newNsampDtFrequencyTime(df,dt_proposed,nsamp,dt)
	real :: df, dt_proposed,dt
	integer :: k,nsamp
!
!  T = 1/df. With proposed dtp we find rn = T/dtp = 1/(df*dtp). To test if rn is a power of 2
!  we calculate the base 2 log: lb(rn) = lg(rn)/lg(2) = -lg(df*dtp)/lg(2). This is usually a real
!  number. We take the smallest integer which is greater or equal lb(rn): k = ceiling(lb(rn)) and use
!  N = 2**k. The new and correct dt for time series length T is then: dt = T/N = 1/(df*N).
!  To make sure that dt = dtp, you must have chosen: df = 1/(dtp*2**k) or T = dtp*2**k
!  when calculating the frequency spectrum.
!
	k = ceiling( -log(df*dt_proposed)/log(2.) )
	nsamp = 2**k
	dt = 1./(nsamp*df)
	end subroutine newNsampDtFrequencyTime
!-----------------------------------------------------------------------------
!  transform spectrum into time domain
!
!  nsamp: power of 2 specifying total length of spectrum and time series
!         determined using newNsampDtFrequencyTime
!
	function transformFrequencyTime(nf1,nf2,df,sigma,sp,nsamp,dt,ts) result(errmsg)
	integer :: nf1,nf2,nsamp,i,nslog,ns,ndec
	real :: df,sigma,dt,fac,expo,dts
	complex, dimension(:) :: sp
	real, dimension(:) :: ts
	double precision, dimension(:), allocatable :: x,y
	type (error_message) :: errmsg
!
!  treat special case where deisred nsamp is smaller than required by nf2
!  if nsamp < 2*(next 2-power of nf2) take the latter else nsamp
!
	call new(errmsg,'transformFrequencyTime')
	nslog = ceiling(alog(real(nf2))/alog(2.))+1
	ns = 2**nslog
	if (ns <= nsamp) then
		ns = nsamp
		nslog = ceiling( -log(df*dt)/log(2.) )
		dts = dt
		ndec = 1
	else if (ns > nsamp) then
		call new(errmsg,1,'Downsampling to less than 1./(2*fny) requested','transformFrequencyTime')
		dts = 1./(df*ns)            ! smaller than dt
		ndec = nint(dt/dts)         ! decimation factor to obtain dt
		if (abs(dt-ndec*dts) > 1.e-5) then
			call new(errmsg,2,'Non-integer ratio of desired dt and dt defined from frequency spectrum', &
			       & 'transformFrequencyTime')
			return
		endif 
	endif
!
!  fill spectrum for FFT
!
	allocate(x(ns),y(ns))
	x(1:nf1) = 0.d0; y(1:nf1) = 0.d0
	x(nf1:nf2) = real(sp(nf1:nf2)); y(nf1:nf2) = aimag(sp(nf1:nf2))
	x(nf2+1:ns/2+1) = 0.d0; y(nf2+1:ns/2+1) = 0.d0
!
!  mirroring
!
	do i=ns/2+2,ns
		x(i)= x(ns-i+2)
		y(i)=-y(ns-i+2)
	enddo
!
!  Fourier transform to time domain
!
	call fastFourierTransform(x,y,nslog,+1)
!
!  copy xs into single variable
!  and multiply by df and multiply by exp(+sigma*t)
!
	fac=1.
	expo=exp(sigma*dt)
	do i=1,nsamp
		ts(i)=x(ndec*(i-1)+1)*df*fac
		fac=fac*expo
	enddo
	deallocate(x,y)
!
	end function transformFrequencyTime
!--------------------------------------------------------------------------------------
!> \brief Transform from frequency into time domain at non-equally spaced time points
!> \param iflag Sign of exponent in Fourier transform
!> \param Desired error (1e-13 < eps < 1e-1)
!> \param nf1 Index of first value in spectrum with f=(nf1-1)*df
!> \param nf2 Index of last value in spectrum with f=(nf2-1)*df
!> \param df Frequency spacing
!> \param sigma Spectrum evaluated at w-i*sigma
!> \param sp Complex Fourier spectrum
!> \param t Array with time points (input)
!> \param y Array with transformed values at time points (output)
!
	function transformNonuniformFrequencyTime(iflag,eps,nf1,nf2,df,sigma,sp,t,y) result(errmsg)
	integer :: iflag,nf1,nf2,ier,i,ms,nt
	real :: df,sigma,eps
	complex, dimension(:) :: sp
	real, dimension(:) :: t,y
	type (error_message) :: errmsg
	double complex, dimension(:), allocatable :: fk,zy
	double precision, dimension(:), allocatable :: xi
	real :: tlen,tmax
	character (len=132) :: myname = 'transformNonuniformFrequencyTime'
!
	call new(errmsg,myname)
!
!  check if maxval(abs(t)) > tlen = 1./df
!  then df is too small
!
	tlen = 1./df
	tmax = maxval(abs(t))
	if (tmax > tlen) then
		call new(errmsg,2,'Maximum abs(t) is larger than 1/df',myname)
		return
	endif
!
	nt = size(t)
	allocate(zy(nt),xi(nt))
	xi = 2.*mc_pid/tlen*t              ! normalized t-values
!
!  fill spectrum for NUFFT
!
	ms = 2*nf2
	allocate(fk(-ms/2:(ms-1)/2))
	fk(0:nf1-2) = 0.
	fk(nf1-1:nf2-1) = sp(nf1:nf2)
	fk(-ms/2) = 0.
!
!  mirroring
!
	do i = 1,ms/2-1
		fk(-i) = conjg(fk(i))
	enddo
!
	call nonUniformOutputSamplesFourierTransform(nt,xi,zy,iflag,dble(eps),ms,fk,ier)
	if (ier == 1) then
		call new(errmsg,2,'Error limits out of range',myname)
		return
	endif
	if (sigma > 1e-6) then
		do i = 1,nt
			y(i) = real(zy(i))*df*exp(sigma*t(i))
		enddo
	else
		y = real(zy)*df              ! Take real part
	endif
	deallocate(fk,zy,xi)
	end function transformNonuniformFrequencyTime
!
 end module frequencyTime	
