!--------------------------------------------------------
!  routines that operate with time series
!--------------------------------------------------------
 module timeSeries
	implicit none
!
 contains
!--------------------------------------------------------
!  calculate envelope of a time series
!
	subroutine envelopeTimeSeries(y,env)
	real, dimension(:) :: y,env
	integer :: nlog,n,n1,i
	real, dimension(:), allocatable :: e, work
	real, external :: pythag
!
!  SFT of data
!
	n1=size(y)
	nlog=ceiling(alog(float(n1))/alog(2.))
	n=2**nlog
	allocate(work(n), e(n))
	work(1:n) = 0.
	e(1:n1) = y(1:n1)
	e(n1+1:n) = 0.
	call ssft(e,work,nlog,-1)
!
!  zero second half of spectrum
!
	e(n/2+2:n)=0.
	work(n/2+2:n)=0.
!
!  Back transform and 2* absvalue
!
	call ssft(e,work,nlog,+2)
	do i=1,n1
		env(i)=2.*pythag(e(i),work(i))
	enddo
	deallocate(work,e)
	end subroutine envelopeTimeSeries
!---------------------------------------------------------
!  calculate sigma from envelope of time series
!  using water level and magnification and steepness parameter
!  see envsigma.f in wolib
!
	subroutine sigmaEnvelopeTimeSeries(y,mag,nexp,sigma)
	real, dimension(:) :: y,sigma
	integer :: nexp,i
	real :: mag,emax,w,x
	real, dimension(:), allocatable :: env
!
	allocate(env(size(y)))
	call envelopeTimeSeries(y,env)
	emax=maxval(env)
	w=emax/mag
	do i=1,size(y)
		x=env(i)/w
		sigma(i)=w*( 1. + x**(nexp+1) )/( x**nexp + w/emax )
	enddo
	deallocate(env)
	end subroutine sigmaEnvelopeTimeSeries
!--------------------------------------------------------------
!  do a cross correlation in the time domain
!
!  k >= 0: c(k) = sum_{i=1}^{min(ns,nd-k)} d(k+i)*s(i)*dt
!  k <  0: c(k) = sum_{i=|k|+1}^{min(nd+|k|,ns)} d(k+i)*s(i)*dt
!
!  ncc:    dimension of c in calling program (input)
!  nd:     number of d samples (input)
!  ns:     number of s samples (input)
!  nc:     number of cross correlation samples (output)
!  d:      d-array (input)
!  s:      s-array (input)
!  c:      cross correlation function (ouptut)
!  kmin:   true index of leftmost sample (output)
!  kmax:   true index of rightmost sample (output)
!  kl:     desired index of leftmost sample of c (minimum = -ns+1) (optional) 
!  kr:     desired index of rightmost sample of c (maximum = nd-1) (optional)
!
!  Note: index k=0 corresponds to zero lag
!        maximum number of samples of cross-correlation is nd+ns-1
!        index k=kmin is mapped to first element of c in calling program
!        index k=0 is mapped to (first element-kmin) of c in calling program if kmin < 0
!        else k=0 is not considered
!        index k=kmax is mapped to (first element+|kmax-kmin|) of c in calling program
!-------------------------------------------------------------------
	subroutine crossCorrelateTimeSeries(d,s,c,dt,kmin,kmax,kl,kr)
	integer, optional :: kl, kr
	real, dimension(:) :: c,d,s
	real :: help,dt
	integer ::i,k,kmin,kmax,nd,ns,nc
!
	nd = size(d); ns = size(s)
	if(.not.present(kl)) then; kmin = -ns+1; else; kmin = max(kl,-ns+1); endif
	if(.not.present(kr)) then; kmax = nd-1;  else; kmax = min(kr, nd-1); endif
	nc=kmax-kmin+1
	if(nc > size(c)) then
		print *,'<crocotd>: allocation of c is to small to take cross-correlation'
		print *,'size needed: ',nc,', size available: ',size(c)
		stop
	endif
	do k=kmin,-1
		help=0.
		do i=iabs(k)+1,min(nd+iabs(k),ns)
			help=help+d(k+i)*s(i)
		enddo
		c(k-kmin+1)=help*dt
	enddo
	do k=max(kmin,0),kmax
		help=0.
		do i=1,min(nd-k,ns)
			help=help+d(k+i)*s(i)
		enddo
		c(k-kmin+1)=help*dt
	enddo
!
	end subroutine crossCorrelateTimeSeries
!--------------------------------------------------------------------------------
!  resample time series and do anti alias filtering
!
!  h:			filter coefficients
!  ndec:		decimation factor
!  d:			data array (input)
!  nd:               number of data to be considered
!  f:			filtered and resmapled output series
!  nf:			new number of data
!
	subroutine resampleTimeSeries(h,ndec,d,nd,f,nf)
	integer :: ndec,nf,nd
	real, dimension(:) :: d,f,h
	integer :: nh,i,j,k,skip
!
!  check length of f
!
	if(size(f) < (nd+1)/ndec) then
		print *,'<resampleTimeSeries>: not enough space allocated for output'
		stop
	endif
!
!  perform the convolution, sum from k+1-nh to min(k,nd) to ensure
!  that a) k-j+1 >= 1 and b) k-j+1 <= nh and c) j <= nd
!  throw away the first (nh-1)/2 samples to compensate for phase shift of filter
!  also throw away the last (nh-1)/2 samples at the end
!  only calculate every ndec'th sample
!
	nh=size(h)
	skip = (nh-1)/2
	i=0
	do k=skip+1,nd+skip,ndec
		i=i+1
		f(i)=0.
		do j=max(k+1-nh,1),min(k,nd)
			f(i) = f(i)+d(j)*h(k-j+1)
		enddo
	enddo
	nf=i
!		
	end subroutine resampleTimeSeries
!--------------------------------------------------------------------------------
!  convolve time series with some filter response
!
!  h:			filter response
!  d:			data array (input)
!  nd:               number of data to be considered
!  f:			filtered output series (same length as d)
!
	subroutine convolveTimeSeries(h,d,nd,f)
	integer :: nd
	real, dimension(:) :: d,f,h
	integer :: nh,i,j,k,skip
!
!  check length of f
!
	if(size(f) < nd) then
		print *,'<convolveTimeSeries>: not enough space allocated for output'
		stop
	endif
!
!  perform the convolution, sum from k+1-nh to min(k,nd) to ensure
!  that a) k-j+1 >= 1 and b) k-j+1 <= nh and c) j <= nd
!  throw away the first (nh-1)/2 samples to compensate for phase shift of filter
!  also throw away the last (nh-1)/2 samples at the end
!
	nh=size(h)
	skip = (nh-1)/2
	do k=skip+1,nd+skip
		i = k-skip
		f(i)=0.
		do j=max(k+1-nh,1),min(k,nd)
			f(i) = f(i)+d(j)*h(k-j+1)
		enddo
	enddo
!		
	end subroutine convolveTimeSeries
!--------------------------------------------------------------------------------
!  differentiate time series
!
!  d:              data array (input)
!  dp:             differentiated array (output)
!  dt:             sampling interval
!
!  dp(i) = (d(i+1)-d(i-1))/(2*dt)
!  dp(1) = dp(2)
!  dp(n) = dp(n-1)
!
	subroutine diffTimeSeries(d,dp,dt)
	real, dimension(:) :: d,dp
	real :: dt
	integer :: n,i
!
	n = size(d)
	do i=2,n-1
		dp(i) = (d(i+1)-d(i-1))/(2.*dt)
	enddo
	dp(1) = dp(2)
	dp(n) = dp(n-1)
!
	end subroutine diffTimeSeries	
!
 end module timeSeries 
