c-------------------------------------------------------
c  Do a causal convolution of the form
c  u(i)=sum_{j=max(1,i-n+1)}^min(i,nq) q(j)*x(i-j+1)
c  i runs from 1 to n+nq-1.
c  Since 1 <= j <= i only the present and past
c  values of x are used to construct u.
c  Here this fact is used to overwrite
c  the array x with the values of u.
c
c  ndim:    dimension of data array
c  n:       number of data
c  nq:      number of filter coefficients
c  ncv:     length of convolution n+nq-1
c  dt:      sampling interval
c  x:       array with data
c  q:       Filter impulse response
c
c  The concolved time series is stored onto x.
c-------------------------------------------------------
	subroutine causfalt(ndim,n,nq,ncv,dt,x,q)
	integer n,nq,i,j,ndim,ncv
	real x(ndim),q(nq),help,dt
c
c  start with the last time sample
c  because this won't be used anymore
c
	if(n+nq-1.gt.ndim) then
		print *,'<causfalt>: dimension of x is to small to take convolution'
		print *,'dimension needed: ',n+nq-1,', dimension available: ',ndim
		stop
	endif
	do i=n+nq-1,1,-1
		help=0.
		do j=max(1,i-n+1),min(i,nq)
			help=help+q(j)*x(i-j+1)
		enddo
		x(i)=help*dt
	enddo
	ncv=n+nq-1
c
	return
	end
