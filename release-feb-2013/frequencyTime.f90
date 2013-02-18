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
!
 end module frequencyTime	
