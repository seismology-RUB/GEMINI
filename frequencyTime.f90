!----------------------------------------------------------
!  routines for transformation from frequency to time domain
!----------------------------------------------------------
 module frequencyTime
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
	subroutine transformFrequencyTime(nf1,nf2,df,sigma,sp,nsamp,dt,ts)
	integer :: nf1,nf2,nsamp,i,nslog
	real :: df,sigma,dt,fac,expo
	complex, dimension(:) :: sp
	real, dimension(:) :: ts
	double precision, dimension(:), allocatable :: x,y
!
!  fill spectrum for FFT
!
	allocate(x(nsamp),y(nsamp))
	x(1:nf1) = 0.d0; y(1:nf1) = 0.d0
	x(nf1:nf2) = real(sp(nf1:nf2)); y(nf1:nf2) = aimag(sp(nf1:nf2))
	x(nf2+1:nsamp/2+1) = 0.d0; y(nf2+1:nsamp/2+1) = 0.d0
!
!  mirroring
!
	do i=nsamp/2+2,nsamp
		x(i)= x(nsamp-i+2)
		y(i)=-y(nsamp-i+2)
	enddo
!
!  Fourier transform to time domain
!
	nslog = ceiling( -log(df*dt)/log(2.) )
	call sft(x,y,nslog,+1)
!
!  copy xs into single variable
!  and multiply by df and multiply by exp(+sigma*t)
!
	fac=1.
	expo=exp(sigma*dt)
	do i=1,nsamp
		ts(i)=x(i)*df*fac
		fac=fac*expo
	enddo
	deallocate(x,y)
!
	end subroutine transformFrequencyTime
!
 end module frequencyTime	
