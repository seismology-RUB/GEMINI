!-----------------------------------------------------------------------
!  Fast Fourier transform routines
!     schnelle fouriertransformation der 2**n komplexen werte (x,y)
!     is negativ - in den frequenzbereich / positiv - in den zeitbereich
!     normierung - abs(is) =1 - ohne / 2 - mit 1/ng / 3 - mit 1/sqrt(ng)
!    4 - ohne normierung,ohne kontrollausdruck
!-----------------------------------------------------------------------
 module fourierTransform
	use mathConstants
	implicit none
	interface fastFourierTransform
		module procedure realFastFourierTransform
		module procedure doubleFastFourierTransform
	end interface
!
 contains
!----------------------------------------------------------------------
!  single precision version (ssft.f)
!
	subroutine realFastFourierTransform(x,y,n,is)
	real, dimension(:) :: x,y
	integer :: is,n,m
	integer, dimension(21) :: zh
	integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
	real :: gn,alpha,piz,beta,excos,exsin,zx,zy
!
	piz=2.*mc_pi
!
!  tabelle der zweierpotenzen
!
	zh(1)=1
	do l=1,n
		zh(l+1)=2*zh(l)
	enddo
	ng=zh(n+1)
	gn=1./float(ng)
!
!  kernprogramm, dreifache schleife ueber schritt/index/teilserie
!
	do m=1,n
		nar=zh(m)
		lar=ng/nar
		larh=lar/2
		alpha =  piz/float(isign(lar,is))
		do jr=1,larh
			beta=alpha*float(jr-1)
			excos = cos(beta)
			exsin = sin(beta)
			ja=jr-lar
			do nr=1,nar
				ja=ja+lar
				jb=ja+larh
				zx = x(ja)-x(jb)
				zy = y(ja)-y(jb)
				x(ja) = x(ja)+x(jb)
				y(ja) = y(ja)+y(jb)
				x(jb) = zx*excos-zy*exsin
				y(jb) = zx*exsin+zy*excos
			enddo
		enddo
	enddo
!
!     normierung
!
	if (iabs(is) == 3) gn=sqrt(gn)
	if (iabs(is) == 2 .or. iabs(is) == 3) then
		y = y*gn
		x = x*gn
	endif
!
!     umordnung nach "bitreversed" indizes
!
	 do j=1,ng
		js=j-1
		k=1
		nny=n+1
		do ny=1,n
			nny=nny-1
			if (js.lt.zh(nny)) cycle
			js=js-zh(nny)
			k=k+zh(ny)
		enddo
		if (j-k < 0) then
			zx = x(j)
			zy = y(j)
			x(j) = x(k)
			y(j) = y(k)
			x(k) = zx
			y(k) = zy
		else
			cycle
		endif
	enddo
	end subroutine realFastFourierTransform
!----------------------------------------------------------------------
!  double precision version (sft.f)
!
	subroutine doubleFastFourierTransform(x,y,n,is)
	double precision, dimension(:) :: x,y
	integer :: is,n,m
	integer, dimension(21) :: zh
	integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
	double precision :: gn,alpha,piz,beta,excos,exsin,zx,zy
!
	piz=2.*mc_pid
!
!  tabelle der zweierpotenzen
!
	zh(1)=1
	do l=1,n
		zh(l+1)=2*zh(l)
	enddo
	ng=zh(n+1)
	gn=1./dble(ng)
!
!  kernprogramm, dreifache schleife ueber schritt/index/teilserie
!
	do m=1,n
		nar=zh(m)
		lar=ng/nar
		larh=lar/2
		alpha =  piz/dble(isign(lar,is))
		do jr=1,larh
			beta=alpha*dble(jr-1)
			excos = dcos(beta)
			exsin = dsin(beta)
			ja=jr-lar
			do nr=1,nar
				ja=ja+lar
				jb=ja+larh
				zx = x(ja)-x(jb)
				zy = y(ja)-y(jb)
				x(ja) = x(ja)+x(jb)
				y(ja) = y(ja)+y(jb)
				x(jb) = zx*excos-zy*exsin
				y(jb) = zx*exsin+zy*excos
			enddo
		enddo
	enddo
!
!     normierung
!
	if (iabs(is) == 3) gn=dsqrt(gn)
	if (iabs(is) == 2 .or. iabs(is) == 3) then
		y = y*gn
		x = x*gn
	endif
!
!     umordnung nach "bitreversed" indizes
!
	 do j=1,ng
		js=j-1
		k=1
		nny=n+1
		do ny=1,n
			nny=nny-1
			if (js.lt.zh(nny)) cycle
			js=js-zh(nny)
			k=k+zh(ny)
		enddo
		if (j-k < 0) then
			zx = x(j)
			zy = y(j)
			x(j) = x(k)
			y(j) = y(k)
			x(k) = zx
			y(k) = zy
		else
			cycle
		endif
	enddo
	end subroutine doubleFastFourierTransform
!----------------------------------------------------------------------
!  Fast Fourier transform for non-uniform samples
!
!
!                         1  nj
!   iflag > 0:  fk(k1) = -- SUM cj(j) exp(+i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
!                        nj j=1                            
!
!                         1  nj
!   iflag < 0:  fk(k1) = -- SUM !j(j) exp(-i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
!                        nj j=1                            
!
!     References:
!
!     [DR] Fast Fourier transforms for nonequispaced data,
!          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14, 
!          1368-1383, 1993.
!
!     [GL] Accelerating the Nonuniform Fast Fourier Transform,
!          L. Greengard and J.-Y. Lee, SIAM Review 46, 443-454 (2004).
!
! ----------------------------------------------------------------------
!     INPUT:
!
!     nj     number of sources   (integer)
!     xj     location of sources (double precision)
!            on interval [-pi,pi].
!
!     cj     strengths of sources (double complex)
!     iflag  determines sign of FFT (see above)
!     eps    precision request  (between 1.0d-13 and 1.0d-1)
!     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
!
!     OUTPUT:
!
!     fk     Fourier transform values (double complex)
!     ier    error return code
!   
!            ier = 0  => normal execution.
!            ier = 1  => precision eps requested is out of range.
!
!     The type 1 NUFFT proceeds in three steps (see [GL]).
!
!     1) spread data to oversampled regular mesh using convolution with
!        a Gaussian 
!     2) compute FFT on uniform mesh
!     3) deconvolve each Fourier mode independently
!          (mutiplying by Fourier transform of Gaussian).
!
! ----------------------------------------------------------------------
!
!     The oversampled regular mesh is defined by 
!
!     nf1 = rat*ms  points, where rat is the oversampling ratio.
!       
!     For simplicity, we set  
!
!         rat = 2 for eps > 1.0d-11
!         rat = 3 for eps <= 1.0d-11.
!
!     The Gaussian used for convolution is:
!
!        g(x) = exp(-x^2 / 4tau) 
!
!     It can be shown [DR] that the precision eps is achieved when
!
!     nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
!     and tau is chosen as
!
!     tau = pi*lambda/(ms**2)
!     lambda = nspread/(rat(rat-0.5)).
!
!     Note that the Fourier transform of g(x) is
!
!     G(s) = exp(-s^2 tau) = exp(-pi*lambda s^2/ms^2)
!
! ----------------------------------------------------------------------
!     Fast Gaussian gridding is based on the following observation.
!
!     Let hx = 2*pi/nf1. In gridding data onto a regular mesh with
!     spacing nf1, we shift the source point xj by pi so 
!     that it lies in [0,2*pi] to simplify the calculations.
!     Since we are viewing the function
!     as periodic, this has no effect on the result.
!    
!     For source (xj+pi), let kb*hx denote the closest grid point and
!     let  kx*hx be a regular grid point within the spreading
!     distance. We can write
!
!     (xj+pi) - kx*hx = kb*hx + diff*hx - kx*hx = diff*hx - (kx-kb)*hx
!
!     where diff = (xj+pi)/hx - kb.
!
!     Let t1 = hx*hx/(4 tau) = pi/(nf1*nf1)/lambda*ms*ms
!                            = pi/lambda/(rat*rat)
!
!     exp(-( (xj+pi) -kx*hx)**2 / 4 tau)
!         = exp(-pi/lamb/rat^2 *(diff - (kx-kb))**2)
!         = exp(-t1 *(diff - (kx-kb))**2)
!         = exp(-t1*diff**2) * exp(2*t1*diff)**k * exp(-t1*k**2)
!           where k = kx-kb.
! 
      subroutine nonUniformInputSamplesFourierTransform(nj,xj,cj,iflag,eps,ms,fk,ier)
      integer :: ier,iflag,istart,iw1,iwtot,iwsav
      integer :: j,jb1,jb1u,jb1d,k1,ms,nf1,nj,nspread
      double precision :: cross,cross1,diff1,eps,hx,rat,r2lamb,t1,tau
      double precision :: xc(-47:47),xj(nj)
      double complex :: cj(nj),fk(-ms/2:(ms-1)/2),zz,ccj
      double precision, allocatable, save :: fw(:)
!
!     rat is oversampling parameter
!     nspread is number of neighbors to which Gaussian gridding is
!     carried out.
!
      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = int(-log(eps)/(mc_pid*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235FourierTransform(2d0*nspread) 
      endif 
!
!     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
!     It is more convenient to define r2lamb = rat*rat*lambda
!
      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*mc_pid/nf1
!
!     Compute workspace size and allocate
!
      iw1 = 2*nf1
      iwsav = iw1+nspread+1
      iwtot = iwsav+4*nf1+15
      allocate ( fw(0:iwtot) )
!
!     Precompute spreading constants and initialize fw
!     to hold one term needed for fast Gaussian gridding 
!
      t1 = mc_pid/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo
      call dcffti(nf1,fw(iwsav))
!
!     Initialize fine grid data to zero.
!
      do k1 = 0, 2*nf1-1
         fw(k1) = 0d0
      enddo
!
!     Loop over sources (1,...,nj)
!
!     1. find closest mesh point (with periodic wrapping if necessary)
!     2. spread source data onto nearest nspread grid points
!        using fast Gaussian gridding.
!
!     The following is a little hard to read because it takes
!     advantage of fast gridding and optimized to minimize the 
!     the number of multiplies in the inner loops.
!
      do j = 1, nj
         ccj = cj(j)/dble(nj)

         jb1 = int((xj(j)+mc_pid)/hx)
         diff1 = (xj(j)+mc_pid)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
!
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo
!
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    istart = 2*(jb1+k1+nf1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = -jb1d, jb1u
	    istart = 2*(jb1+k1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = jb1u+1, nspread
	    istart = 2*(jb1+k1-nf1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
      enddo
!
!     Compute 1D FFT and carry out deconvolution.
!     There is a factor of (-1)**k1 needed to account for the 
!     FFT phase shift.
!
      if (iflag .ge. 0) then
         call dcfftb(nf1,fw(0),fw(iwsav))
      else
         call dcfftf(nf1,fw(0),fw(iwsav))
      endif
!
      tau = mc_pid * r2lamb / dble(nf1)**2
      cross1 = 1d0/sqrt(r2lamb)
      zz = dcmplx(fw(0),fw(1))
      fk(0) = cross1*zz
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(tau*dble(k1)**2)
	 zz = dcmplx(fw(2*k1),fw(2*k1+1))
         fk(k1) = cross*zz
	 zz = dcmplx(fw(2*(nf1-k1)),fw(2*(nf1-k1)+1))
         fk(-k1) = cross*zz
      enddo
      if (ms/2*2.eq.ms) then
         cross = -cross1*exp(tau*dble(ms/2)**2)
         zz = dcmplx(fw(2*nf1-ms),fw(2*nf1-ms+1))
         fk(-ms/2) = cross*zz
      endif
      deallocate(fw)
      end subroutine nonUniformInputSamplesFourierTransform
!-----------------------------------------------------------------------
!  Fast Fourier transform for non-uniform output samples and regular
!  input samples
!
!              (ms-1)/2
!     cj(j) =    SUM      fk(k1) exp(+i k1 xj(j))  for j = 1,...,nj
!              k1= -ms/2                            
!
!              (ms-1)/2
!     cj(j) =    SUM      fk(k1) exp(-i k1 xj(j))  for j = 1,...,nj
!              k1= -ms/2                            
!
! ----------------------------------------------------------------------
!     INPUT:
!
!     nj     number of output values   (integer)
!     xj     location of output values (double precision array)
!     iflag  determines sign of FFT (see above)
!     eps    precision request  (between 1.0d-13 and 1.0d-1)
!     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
!     fk     Fourier coefficient values (double complex array)
!
!     OUTPUT:
!
!     cj     output values (double complex array)
!     ier    error return code
!   
!            ier = 0  => normal execution.
!            ier = 1  => precision eps requested is out of range.
!
!
!     The type 2 algorithm proceeds in three steps (see [GL]).
!
!     1) deconvolve (amplify) each Fourier mode first 
!     2) compute inverse FFT on uniform fine grid
!     3) spread data to regular mesh using Gaussian
!
!
!     See subroutine nufft1d1f90(nj,xj,cj,iflag,eps,ms,fk,ier)
!     for more comments on fast gridding and parameter selection.
!
      subroutine nonUniformOutputSamplesFourierTransform(nj,xj,cj,iflag,eps,ms,fk,ier)
      integer :: ier,iflag,iw1,iwsav,iwtot,j,jb1,jb1u,jb1d,k1
      integer :: ms,nf1,nj,nspread
      double precision :: cross,cross1,diff1,eps,hx,rat,r2lamb,t1
      double precision :: xj(nj),xc(-47:47)
      double complex :: cj(nj), fk(-ms/2:(ms-1)/2)
      double complex :: zz
      double precision, allocatable, save :: fw(:)
!
!     Precision dependent parameters
!     rat is oversampling parameter
!     nspread is number of neighbors to which Gaussian gridding is
!     carried out.
!
      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = int(-log(eps)/(mc_pid*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235FourierTransform(2d0*nspread) 
      endif 
!
!     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
!     It is more convenient to define r2lamb = rat*rat*lambda
!
      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*mc_pid/nf1
!
!     Compute workspace size and allocate
!
      iw1 = 2*nf1
      iwsav = iw1 + nspread+1
      iwtot = iwsav + 4*nf1 + 15
      allocate ( fw(0:iwtot))
!
!     Precompute spreading constants and initialize fw
!     to hold one term needed for fast Gaussian gridding 
!
      t1 = mc_pid/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo
      call dcffti(nf1,fw(iwsav))
!
!     Deconvolve and compute inverse 1D FFT
!     (A factor of (-1)**k is needed to shift phase.)
!
      t1 = mc_pid * r2lamb / dble(nf1)**2
      cross1 = 1d0/sqrt(r2lamb)
      zz = cross1*fk(0)
      fw(0) = dreal(zz)
      fw(1) = dimag(zz)
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(t1*dble(k1)**2)
         zz = cross*fk(k1)
         fw(2*k1) = dreal(zz)
         fw(2*k1+1) = dimag(zz)
         zz = cross*fk(-k1)
         fw(2*(nf1-k1)) = dreal(zz)
         fw(2*(nf1-k1)+1) = dimag(zz)
      enddo
      cross = -cross1*exp(t1*dble(ms/2)**2)
      if (ms/2*2.eq.ms) then
	 zz = cross*fk(-ms/2)
         fw(2*nf1-ms) = dreal(zz)
         fw(2*nf1-ms+1) = dimag(zz)
      endif
      do k1 = (ms+1)/2, nf1-ms/2-1
         fw(2*k1) = 0d0
         fw(2*k1+1) = 0d0
      enddo
!
      if (iflag .ge. 0) then
         call dcfftb(nf1,fw(0),fw(iwsav))
      else
         call dcfftf(nf1,fw(0),fw(iwsav))
      endif
!
!     Loop over target points (1,...,nj)
!       1. find closest mesh point (with periodic wrapping if needed)
!       2. get contributions from regular fine grid to target
!          locations using Gaussian convolution.
!
      t1 = mc_pid/r2lamb
      do j = 1, nj
         cj(j) = dcmplx(0d0,0d0)
         jb1 = int((xj(j)+mc_pid)/hx)
         diff1 = (xj(j)+mc_pid)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo
!
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    zz = dcmplx(fw(2*(jb1+k1+nf1)),fw(2*(jb1+k1+nf1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
         do k1 = -jb1d, jb1u
	    zz = dcmplx(fw(2*(jb1+k1)),fw(2*(jb1+k1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
         do k1 = jb1u+1, nspread
	    zz = dcmplx(fw(2*(jb1+k1-nf1)),fw(2*(jb1+k1-nf1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
      enddo
      deallocate(fw)
      end subroutine nonUniformOutputSamplesFourierTransform
!------------------------------------------------------------------------
!  returns next multiple of 2,3 and 5 greater or equal than base
!  next235 = 2^p 3^q 5^r >= base  where p>=1, q>=0, r>=0
!
	integer function next235FourierTransform(base) result(next235)
	integer :: numdiv
	double precision :: base
	integer :: next235
!
	next235 = 2 * int(base/2d0+.9999d0)
	if (next235.le.0) next235 = 2
 100	numdiv = next235
	do while (numdiv/2*2 .eq. numdiv)
		numdiv = numdiv /2
	enddo
	do while (numdiv/3*3 .eq. numdiv)
		numdiv = numdiv /3
	enddo
	do while (numdiv/5*5 .eq. numdiv)
		numdiv = numdiv /5
	enddo
	if (numdiv .eq. 1) return
	next235 = next235 + 2
	goto 100
	end function next235FourierTransform
!
 end module fourierTransform
