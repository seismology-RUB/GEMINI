c-------------------------------------------------------------
c  seismo_oversampling:
c  
c  set true dt and nsamp
c  df: frequency spacing
c  dt_proposed: desired sampling interval
c-------------------------------------------------------------
	subroutine seismo_oversampling(df,dt_proposed)
	include 'seismodim.h'
	include 'seismo.h'
	real fslog,df,dt_proposed
	integer ceil
c	
	fslog=-alog(2.*df*dt_proposed)/alog(2.)+1.
	seismo_nslog=ceil(fslog)
	seismo_nsamp=2**(seismo_nslog)
	seismo_dt=1./(df*seismo_nsamp)
	if(seismo_nsamp.gt.nsampp) then
		print *,'<seismo_oversampling>: more samples than dimensioned'
		print *,seismo_dt,seismo_nsamp
		stop
	endif
	return
	end	
c-------------------------------------------------------------
c  seismo_f2t:
c
c  Compute seismograms from complex frequency spectrum
c
c  nf1:      fmin=(nf1-1)*df
c  nf2:      fmax=(nf2-1)*df
c  df:       frequency spacing
c  ifwin:    flag to cut seismogram according to
c  tmax:     length of seismogram
c  sigma:    negative imaginary part of frequency
c  zdis:     complex array with frequency spectrum
c  nsampw:   number of samples (output)
c  urs:      array for time series (output)
c------------------------------------------------------------
	subroutine seismo_f2t(nf1,nf2,df,ifwin,tmax,sigma,zdis,nsampw,urs)
	include 'seismodim.h'
	include 'seismo.h'
	integer nf1,nf2,ifwin
	real df,tmax,sigma,urs(nsampp)
	double precision xs(nsampp),ys(nsampp)
	complex zdis(nf2)
	integer nsampw,it1,it2,if,it
	real tw2,fac,expo
c
c  apply time window if set
c
	if(ifwin.eq.1) then
		tw2=tmax
		it1=1
		it2=min(tw2/seismo_dt,seismo_nsamp-1)+1
	else
		it1=1
		it2=seismo_nsamp
	endif
	nsampw=it2-it1+1
c
c  fill spectrum for FFT
c
	do if=1,nf1-1
		xs(if)=0.
		ys(if)=0.
	enddo
	do if=nf1,nf2
		xs(if)=real(zdis(if))
		ys(if)=aimag(zdis(if))
	enddo
	do if=nf2+1,seismo_nsamp/2+1
		xs(if)=0.
		ys(if)=0.
	enddo
c
c  mirroring
c
	do if=seismo_nsamp/2+2,seismo_nsamp
		xs(if)=xs(seismo_nsamp-if+2)
		ys(if)=-ys(seismo_nsamp-if+2)
	enddo
c
c  Fourier transform to time domain
c
	call sft(xs,ys,seismo_nslog,+1)
c
c  copy xs into single variable, decimate
c  and multiply by df and multiply by exp(+sigma*t)
c-
	fac=exp(+sigma*(it1-1)*seismo_dt)
	expo=exp(sigma*seismo_dt)
	do it=it1,it2
		urs(it-it1+1)=xs(it)*df*fac
		fac=fac*expo
	enddo
c
	return
	end
c-------------------------------------------------------------
c  seismo_stffromhdurtriangle
c
c  compute a source time function from a given half duration
c  using a triangle function with ramp of length hdur
c
c  This is a moment rate function. If convolved with displacement
c  Green functions, one obtains a velocity seismogram
c
c  normalize source time function to unit area, i.e division
c  by half duration
c
c  hdur: half duration in seconds
c  nstf: number of samples of source time function
c  stf(nsampp): values of source time function (normalized to 1)
c--------------------------------------------------------------
	subroutine seismo_stffromhdurtriangle(hdur,nstf,stf,newhdur)
	include 'seismodim.h'
	include 'seismo.h'
	real stf(nsampp),hdur,newhdur
	integer nstf,imax,i
c
c  if hdur is set to zero, an impules response is caluclated
c  with hdur=0 we get imax=1, nstf=1 and stf(1)=1/dt
c
	imax=nint(hdur/seismo_dt)
	newhdur=imax*seismo_dt
	print *,'Use new half duration: ',newhdur,' instead of ',hdur
	if(imax.eq.0) then
		nstf=1
		stf(1)=1./seismo_dt
		print *,'Half-duration is essentially zero, use stf(1) = 1/dt'
	else
		do i=1,imax+1
			stf(i)=float(i-1)/float(imax)/newhdur
		enddo
		do i=imax+2,2*imax+1
			stf(i)=(1.-float(i-1-imax)/float(imax))/newhdur
		enddo
		nstf=2*imax+1
	endif
	return
	end
c-------------------------------------------------------------
c  seismo_stffromhdursinq
c
c  compute a source time function from a given half duration
c  using a sin**2 function with length 2*hdur
c
c  This is a moment rate function. If convolved with displacement
c  Green functions, one obtains a velocity seismogram
c
c  normalize source time function to unit area, i.e division
c  by half duration
c
c  hdur: half duration in seconds
c  nstf: number of samples of source time function
c  stf(nsampp): values of source time function (normalized to 1)
c--------------------------------------------------------------
	subroutine seismo_stffromhdursinq(hdur,nstf,stf,newhdur)
	include 'seismodim.h'
	include 'seismo.h'
	include 'pis.h'
	real stf(nsampp),hdur,newhdur
	integer nstf,imax,i
c
c  if hdur is set to 0, an impulse response is calculated.
c  this is achieved by setting imax=1 and stf(1)=1./dt
c
	imax=nint(hdur/seismo_dt)
	newhdur=imax*seismo_dt
	print *,'Use new half duration: ',newhdur,' instead of ',hdur
	if(imax.eq.0) then
		nstf=1
		stf(1)=1./seismo_dt
		print *,'Half-duration is essentially zero, use stf(1) = 1/dt'
	else
		do i=1,2*imax+1
			stf(i)=sin((i-1)*seismo_dt*pi*.5/newhdur)**2/newhdur
		enddo
		nstf=2*imax+1
	endif
	return
	end
