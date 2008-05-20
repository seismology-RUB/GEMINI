c-----------------------------------------------------------
c  Routines associated with a
c  Green kf-spectrum
c----------------------------------------------------------
c  Read the header info associated with Green kf-spectrum
c  and leave file open for further reading
c
	subroutine green_readkfhead(iunit,greenfile)
	include 'greendim.h'
	include 'green.h'
	integer iunit
	character*(*) greenfile
c
	open(iunit,file=greenfile,form='unformatted')
	read(iunit) green_zs,green_ze,green_sigma,green_nf1,green_nf2,green_df,green_dwn,
     1	      green_rearth,green_jump,green_comp,green_sourcetype,green_nwnmax
	green_nf=green_nf2-green_nf1+1
c
	return
	end
c----------------------------------------------------------
c  Read a Green kf-spectrum from already open file
c  and close file
c
	subroutine green_readkfsp(iunit,greenkfsp)
	include 'greendim.h'
	include 'green.h'
	include 'nkfmax.h'
	integer iunit,if,iwn,nwe
	complex greenkfsp(nkfmax)
c
	nwe=0
	do if=green_nf1,green_nf2
		read(iunit) green_nwn(if),(greenkfsp(nwe+iwn),iwn=1,green_nwn(if))
		nwe=nwe+green_nwn(if)
	enddo
c
	return
	end
c-----------------------------------------------------------
c  Get number of frequencies
c
	subroutine green_getnf(n)
	include 'greendim.h'
	include 'green.h'
	integer n
	n=green_nf2-green_nf1+1
	return
	end
c-------------------------------------------------------------
c  Compute a cos**2 taper for given frequency
c
	subroutine green_taper(if,tapfrac,greentap)
	include 'greendim.h'
	include 'green.h'
	real greentap(npp)
	real tapfrac,c1,pi,wn,c2
	integer if,iwn1,iwn
c
	pi=4.*atan(1.)
	iwn1=min(nint(green_nwn(if)*(1.-tapfrac)),green_nwn(if))
	iwn1=max(1,iwn1)
	c1=(green_nwn(if)-1)*green_dwn*tapfrac*2./pi
	c2=(1.-tapfrac)/tapfrac*0.5*pi
	do iwn=1,iwn1-1
		greentap(iwn)=1.
	enddo
	wn=(iwn1-1)*green_dwn
	do iwn=iwn1,green_nwn(if)
		greentap(iwn)=cos( wn/c1-c2 )**2
		wn=wn+green_dwn
	enddo
	return
	end
c-----------------------------------------------------------
c  Compute wavenumber integrals for single force source
c
c  kint= 1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
c  kint= 2:       R^2/(2*pi)*int ( dk k (-2/kR) U^2 J_1 )  
c  kint= 3:       R^2/(2*pi)*int ( dk k (-kR) V^1 J_1 )  
c  kint= 4:       R^2/(2*pi)*int ( dk k (-2 V^2)(J_0-J_1/(kx)) )  
c  kint= 5:       R^2/(2*pi)*int ( dk k (-2 W^1) J_1/(kx) )  
c  kint= 6:       R^2/(2*pi)*int ( dk k (+2 W^1)(J_0-J_1/(kx)) )  
c  kint= 7:       R^2/(2*pi)*int ( dk k (+2 V^2) J_1/(kx) )
c
c  kint= 8:       R^2/(2*pi)*int ( dk k U^1_r J_0 )  
c  kint= 9:       R^2/(2*pi)*int ( dk k (-2/kR) U^2_r J_1 )  
c  kint=10:       R^2/(2*pi)*int ( dk k (-kR) V^1_r J_1 )  
c  kint=11:       R^2/(2*pi)*int ( dk k (-2 V^2_r)(J_0-J_1/(kx)) )  
c  kint=12:       R^2/(2*pi)*int ( dk k (-2 W^1_r) J_1/(kx) )  
c  kint=13:       R^2/(2*pi)*int ( dk k (+2 W^1_r)(J_0-J_1/(kx)) )  
c  kint=14:       R^2/(2*pi)*int ( dk k (+2 V^2_r) J_1/(kx) )
c
c  kint=15:       R^2/(2*pi)*int ( dk k U^1 (-kR) J_1 )  
c  kint=16:       R^2/(2*pi)*int ( dk k -2 U^2 (J_0-J_1/(kx)) )  
c  kint=17:       R^2/(2*pi)*int ( dk k -(kR)^2 V^1 (J_0-J_1/(kx)) )  
c  kint=18:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR)(-J_1+J_2/(kx)) )  
c  kint=19:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR) J_2/(kx) )  
c  kint=20:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR)(-J_1+J_2/(kx)) )  
c  kint=21:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR) J_2(kx)/(kx) )
c
c
c  for one receiver at x and ALL frequencies. This saves time
c  because the Bessel functions do not depend on frequency
c  and the wavenumber sampling is independent of frequency.
c  A complex array zdis(f) is returned containing the
c  value of the integral.
c  tapfrac is the fraction of wavenumbers multiplied by a cosine taper.
c
c  Note: Unit of zdis(if) is returned in nanometers/N
c
	subroutine green_wnint_force(kint,x,tapfrac,greenkfsp,besselj,zdis)
	include 'greendim.h'
	include 'green.h'
	include 'nkfmax.h'
	real x,wn3,wn,wnbesj(npp),wn3besy,xdk,wnbesy,greentap(npp),tapfrac,pi,besselj(npp,0:2),by0,by1,by2
	integer kint,j3,iwn,if,j33,j,nwn,nwe
	complex zdis(nff),zsum,greenkfsp(nkfmax)
	real bessy0,bessy1
	external bessy0,bessy1
c
	pi=4.*atan(1.)
	if(x.lt.1.e-6) x=1.e-6
	wn3=3./x
c
c  find index of first wn-point greater than 3
c
	j3=int(wn3/green_dwn)+2
	wn3=(j3-1)*green_dwn
c
c  Bessel terms over whole interval
c  division by zero is unproblematic because we avoid wn=0
c  Bessel terms at k=0 vanish anyway
c  and x shouldn't be zero either
c
	wn=green_dwn
	do iwn=2,green_nwnmax
		if(kint.eq.1.or.kint.eq.8)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.2.or.kint.eq.9)  wnbesj(iwn)=-2*besselj(iwn,1)/green_rearth
		if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wnbesj(iwn)=-wn*wn*besselj(iwn,1)*green_rearth
		if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wnbesj(iwn)=-2.*wn*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.5.or.kint.eq.12) wnbesj(iwn)=-2.*wn*besselj(iwn,1)/(wn*x)
		if(kint.eq.6.or.kint.eq.13) wnbesj(iwn)=+2.*wn*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.7.or.kint.eq.14) wnbesj(iwn)=+2.*wn*besselj(iwn,1)/(wn*x)
		if(kint.eq.17) wnbesj(iwn)=-wn*(wn*green_rearth)**2*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.18) wnbesj(iwn)=-2*wn*wn*green_rearth*( -besselj(iwn,1)+besselj(iwn,2)/(wn*x) )
		if(kint.eq.19) wnbesj(iwn)=+2*wn*wn*green_rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.20) wnbesj(iwn)=+2*wn*wn*green_rearth*( -besselj(iwn,1)+besselj(iwn,2)/(wn*x) )
		if(kint.eq.21) wnbesj(iwn)=-2*wn*wn*green_rearth*besselj(iwn,2)/(wn*x)
		wn=wn+green_dwn
	enddo
c
c  Neumann terms at wn3
c
	by0=bessy0(wn3*x)
	by1=bessy1(wn3*x)
	by2=-by0+2./(wn3*x)*by1
	if(kint.eq.1.or.kint.eq.8) wn3besy=wn3*by0
	if(kint.eq.2.or.kint.eq.9) wn3besy=-2*by1/green_rearth
	if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wn3besy=-wn3*wn3*by1*green_rearth
	if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wn3besy=-2.*wn3*(by0-by1/(wn3*x))
	if(kint.eq.5.or.kint.eq.12) wn3besy=-2.*wn3*by1/(wn3*x)
	if(kint.eq.6.or.kint.eq.13) wn3besy=+2.*wn3*(by0-by1/(wn3*x))
	if(kint.eq.7.or.kint.eq.14) wn3besy=+2.*wn3*by1/(wn3*x)
	if(kint.eq.17) wn3besy=-wn3*(wn3*green_rearth)**2*(by0-by1/(wn3*x))
	if(kint.eq.18) wn3besy=-2*wn3*wn3*green_rearth*(-by1+by2/(wn3*x))
	if(kint.eq.19) wn3besy=+2*wn3*wn3*green_rearth*by2/(wn3*x)
	if(kint.eq.20) wn3besy=+2*wn3*wn3*green_rearth*(-by1+by2/(wn3*x))
	if(kint.eq.21) wn3besy=-2*wn3*wn3*green_rearth*by2/(wn3*x)
c
c  perform integration for each frequency
c  Use trapezoidal rule if wn*x < 3 and Filon else
c
	nwe=0
	do if=green_nf1,green_nf2
		nwn=green_nwn(if)
		call green_taper(if,tapfrac,greentap)
		j33=min(j3,nwn)
		zdis(if)=0.
		do j=2,j33-1
			zdis(if)=zdis(if)+greenkfsp(nwe+j)*wnbesj(j)*greentap(j)
		enddo
		zdis(if)=zdis(if)+0.5*greenkfsp(nwe+j33)*wnbesj(j33)*greentap(j33)
		zdis(if)=zdis(if)*green_dwn
		if(j33.eq.nwn) goto 11
c
c  use the Filon integration for kx > 3
c  Neumann terms at the end of integration interval (k_N)
c
		wn=(green_nwn(if)-1)*green_dwn
		by0=bessy0(wn*x)
		by1=bessy1(wn*x)
		by2=-by0+2./(wn*x)*by1
		if(kint.eq.1.or.kint.eq.8) wnbesy=wn*by0
		if(kint.eq.2.or.kint.eq.9) wnbesy=-2*by1/green_rearth
		if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wnbesy=-wn*wn*by1*green_rearth
		if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wnbesy=-2.*wn*(by0-by1/(wn*x))
		if(kint.eq.5.or.kint.eq.12) wnbesy=-2.*wn*by1/(wn*x)
		if(kint.eq.6.or.kint.eq.13) wnbesy=+2.*wn*(by0-by1/(wn*x))
		if(kint.eq.7.or.kint.eq.14) wnbesy=+2.*wn*by1/(wn*x)
		if(kint.eq.17) wnbesy=-wn*(wn*green_rearth)**2*(by0-by1/(wn*x))
		if(kint.eq.18) wnbesy=-2*wn*wn*green_rearth*(-by1+by2/(wn*x))
		if(kint.eq.19) wnbesy=+2*wn*wn*green_rearth*by2/(wn*x)
		if(kint.eq.20) wnbesy=+2*wn*wn*green_rearth*(-by1+by2/(wn*x))
		if(kint.eq.21) wnbesy=-2*wn*wn*green_rearth*by2/(wn*x)
c
		xdk=x*green_dwn
		zsum=greenkfsp(nwe+j3)*wnbesj(j3)*greentap(j3)
		do j=j3+1,nwn-1
			zsum=zsum+2.*greenkfsp(nwe+j)*wnbesj(j)*greentap(j)
		enddo
		zsum=zsum+greenkfsp(nwe+nwn)*wnbesj(nwn)*greentap(nwn)
		zdis(if)=zdis(if)+zsum/x*(1-cos(xdk))/xdk
		zdis(if)=zdis(if)+(greenkfsp(nwe+nwn)*wnbesy*greentap(nwn)
     1		         -greenkfsp(nwe+j3)*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
c
c  multiply by R^2/(2*pi)
c
 11		zdis(if)=zdis(if)*green_rearth**2/(2.*pi)
c
c  convert to nanometers per N
c
		zdis(if)=zdis(if)*1.e-3
c
c  update address of last wavenumber spectrum value
c	
		nwe=nwe+nwn
	enddo
	return
	end
c-----------------------------------------------------------
c  Compute wavenumber integrals for moment tensor source
c
c  kint=1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
c  kint=2:       R^2/(2*pi)*int ( dk k U^2 J_0 )
c  kint=3:       R^2/(2*pi)*int ( dk k U^3 2/(kR) J_1 )
c  kint=4:       R^2/(2*pi)*int ( dk k U^4 2J_2 )
c  kint=5:       R^2/(2*pi)*int ( dk k V^1 -kR J_1 )
c  kint=6:       R^2/(2*pi)*int ( dk k V^2 -kR J_1 )
c  kint=7:       R^2/(2*pi)*int ( dk k V^3 2 (J_0-J_1/(kx)) )
c  kint=8:       R^2/(2*pi)*int ( dk k V^4 2kR (J_1-2J_2/(kx)) )
c  kint=9:       R^2/(2*pi)*int ( dk k W^1 2J_1/(kx) )
c  kint=10:      R^2/(2*pi)*int ( dk k W^2 4kRJ_2/(kx) )
c  kint=11:      R^2/(2*pi)*int ( dk k W^1 2(J_0-J_1/(kx)) )
c  kint=12:      R^2/(2*pi)*int ( dk k W^2 -2kR(J_1-2J_2/(kx)) )
c  kint=13:      R^2/(2*pi)*int ( dk k V^3 2J_1/(kx) )
c  kint=14:      R^2/(2*pi)*int ( dk k V^4 -4kRJ_2/(kx) )
c
c  kint=15:      R^2/(2*pi)*int ( dk k U^1_r J_0 )
c  kint=16:      R^2/(2*pi)*int ( dk k U^2_r J_0 )
c  kint=17:      R^2/(2*pi)*int ( dk k U^3_r 2/(kR) J_1 )
c  kint=18:      R^2/(2*pi)*int ( dk k U^4_r 2J_2 )
c  kint=19:      R^2/(2*pi)*int ( dk k V^1_r -kR J_1 )
c  kint=20:      R^2/(2*pi)*int ( dk k V^2_r -kR J_1 )
c  kint=21:      R^2/(2*pi)*int ( dk k V^3_r 2(J_0-J_1/(kx)) )
c  kint=22:      R^2/(2*pi)*int ( dk k V^4_r 2kR(J_1-2J_2/(kx)) )
c  kint=23:      R^2/(2*pi)*int ( dk k W^1_r 2J_1/(kx) )
c  kint=24:      R^2/(2*pi)*int ( dk k W^2_r 4kRJ_2/(kx) )
c  kint=25:      R^2/(2*pi)*int ( dk k W^1_r 2(J_0-J_1/(kx)) )
c  kint=26:      R^2/(2*pi)*int ( dk k W^2_r -2kR(J_1-2J_2/(kx)) )
c  kint=27:      R^2/(2*pi)*int ( dk k V^3_r 2J_1/(kx) )
c  kint=28:      R^2/(2*pi)*int ( dk k V^4_r -4kRJ_2/(kx) )
c
c  kint=29:      R^2/(2*pi)*int ( dk k U^1 -kR J_1 )  
c  kint=30:      R^2/(2*pi)*int ( dk k U^2 -kR J_1 )
c  kint=31:      R^2/(2*pi)*int ( dk k U^3 2 (J_0-J_1/(kx)) )
c  kint=32:      R^2/(2*pi)*int ( dk k U^4 2kR (J_1-2J_2/(kx)) )
c  kint=33:      R^2/(2*pi)*int ( dk k V^1 -(kR)^2 (J_0-J_1/(kx)) )
c  kint=34:      R^2/(2*pi)*int ( dk k V^2 -(kR)^2 (J_0-J_1/(kx)) )
c  kint=35:      R^2/(2*pi)*int ( dk k V^3 -2kR (J_1-J_2/(kx)) )
c  kint=36:      R^2/(2*pi)*int ( dk k V^4 2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
c  kint=37:      R^2/(2*pi)*int ( dk k W^1 -2kR J_2/(kx) )
c  kint=38:      R^2/(2*pi)*int ( dk k W^2 4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
c  kint=39:      R^2/(2*pi)*int ( dk k W^1 -2kR (J_1-J_2/(kx)) )
c  kint=40:      R^2/(2*pi)*int ( dk k W^2 -2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
c  kint=41:      R^2/(2*pi)*int ( dk k V^3 -2kR J_2/(kx) )
c  kint=42:      R^2/(2*pi)*int ( dk k V^4 -4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
c
c  for one receiver at x and ALL frequencies. This saves time
c  because the Bessel functions do not depend on frequency
c  and the wavenumber sampling is independent of frequency.
c  A complex array zdis(f) is returned containing the
c  value of the integral.
c  tapfrac is the fraction of wavenumbers multiplied by a cosine taper.
c
c  Note: Unit of zdis(if) is returned in nanometers/(Nm)
c
	subroutine green_wnint_mom(kint,x,tapfrac,greenkfsp,besselj,zdis)
	include 'greendim.h'
	include 'green.h'
	include 'nkfmax.h'
	real x,wn3,wn,wnbesj(npp),wn3besy,xdk,wnbesy,greentap(npp),tapfrac,pi,besselj(npp,0:2),by0,by1,by2
	integer kint,j3,iwn,if,j33,j,nwn,nwe
	complex zdis(nff),zsum,greenkfsp(nkfmax)
	real bessy0,bessy1
	external bessy0,bessy1
c
	pi=4.*atan(1.)
	if(x.lt.1.e-6) x=1.e-6
	wn3=3./x
c
c  find index of first wn-point greater than 3
c
	j3=int(wn3/green_dwn)+2
	wn3=(j3-1)*green_dwn
c
c  Bessel terms over whole interval
c  division by zero is unproblematic because we avoid wn=0
c  Bessel terms at k=0 vanish anyway
c  and x shouldn't be zero either
c
	wn=green_dwn
	do iwn=2,green_nwnmax
		if(kint.eq.1.or.kint.eq.15)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.2.or.kint.eq.16)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.3.or.kint.eq.17)  wnbesj(iwn)=2./green_rearth*besselj(iwn,1)
		if(kint.eq.4.or.kint.eq.18)  wnbesj(iwn)=wn*2.*besselj(iwn,2)
		if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wnbesj(iwn)=wn*(-wn*green_rearth)*besselj(iwn,1)
		if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wnbesj(iwn)=wn*(-wn*green_rearth)*besselj(iwn,1)
		if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wnbesj(iwn)=wn*2.*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wnbesj(iwn)=wn*2.*wn*green_rearth*(besselj(iwn,1)-2.*besselj(iwn,2)/(wn*x))
		if(kint.eq.9.or.kint.eq.23)  wnbesj(iwn)=wn*2.*besselj(iwn,1)/(wn*x)
		if(kint.eq.10.or.kint.eq.24) wnbesj(iwn)=wn*4.*wn*green_rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.11.or.kint.eq.25) wnbesj(iwn)=wn*2.*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.12.or.kint.eq.26) wnbesj(iwn)=-wn*2.*wn*green_rearth*(besselj(iwn,1)-2.*besselj(iwn,2)/(wn*x))
		if(kint.eq.13.or.kint.eq.27) wnbesj(iwn)=wn*2.*besselj(iwn,1)/(wn*x)
		if(kint.eq.14.or.kint.eq.28) wnbesj(iwn)=-wn*4.*wn*green_rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.33.or.kint.eq.34) wnbesj(iwn)=-wn*(wn*green_rearth)**2*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.35.or.kint.eq.39) wnbesj(iwn)=-wn*2*wn*green_rearth*(besselj(iwn,1)-besselj(iwn,2)/(wn*x))
		if(kint.eq.36) wnbesj(iwn)=+wn*2*(wn*green_rearth)**2*((6./(wn*x)**2-1.)*besselj(iwn,2)-besselj(iwn,1)/(wn*x))
		if(kint.eq.37.or.kint.eq.41) wnbesj(iwn)=-wn*2*wn*green_rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.38) wnbesj(iwn)=+wn*4.*(wn*green_rearth)**2*(-3.*besselj(iwn,2)/(wn*x)**2+besselj(iwn,1)/(wn*x))
		if(kint.eq.40) wnbesj(iwn)=-wn*2*(wn*green_rearth)**2*((6./(wn*x)**2-1.)*besselj(iwn,2)-besselj(iwn,1)/(wn*x))
		if(kint.eq.42) wnbesj(iwn)=-wn*4.*(wn*green_rearth)**2*(-3.*besselj(iwn,2)/(wn*x)**2+besselj(iwn,1)/(wn*x))
		wn=wn+green_dwn
	enddo
c
c  Neumann terms at wn3
c
	by0=bessy0(wn3*x)
	by1=bessy1(wn3*x)
	by2=-by0+2./(wn3*x)*by1
	if(kint.eq.1.or.kint.eq.15)  wn3besy=wn3*by0
	if(kint.eq.2.or.kint.eq.16)  wn3besy=wn3*by0
	if(kint.eq.3.or.kint.eq.17)  wn3besy=2./green_rearth*by1
	if(kint.eq.4.or.kint.eq.18)  wn3besy=wn3*2.*by2
	if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wn3besy=wn3*(-wn3*green_rearth)*by1
	if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wn3besy=wn3*(-wn3*green_rearth)*by1
	if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wn3besy=wn3*2.*(by0-by1/(wn3*x))
	if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wn3besy=wn3*2.*wn3*green_rearth*(by1-2.*by2/(wn3*x))
	if(kint.eq.9.or.kint.eq.23)  wn3besy=wn3*2.*by1/(wn3*x)
	if(kint.eq.10.or.kint.eq.24) wn3besy=wn3*4.*wn3*green_rearth*by2/(wn3*x)
	if(kint.eq.11.or.kint.eq.25) wn3besy=wn3*2.*(by0-by1/(wn3*x))
	if(kint.eq.12.or.kint.eq.26) wn3besy=-wn3*2.*wn3*green_rearth*(by1-2.*by2/(wn3*x))
	if(kint.eq.13.or.kint.eq.27) wn3besy=wn3*2.*by1/(wn3*x)
	if(kint.eq.14.or.kint.eq.28) wn3besy=-wn3*4.*wn3*green_rearth*by2/(wn3*x)
	if(kint.eq.33.or.kint.eq.34) wn3besy=-wn3*(wn3*green_rearth)**2*(by0-by1/(wn3*x))
	if(kint.eq.35.or.kint.eq.39) wn3besy=-wn3*2*wn3*green_rearth*(by1-by2/(wn3*x))
	if(kint.eq.36) wn3besy=+wn3*2*(wn3*green_rearth)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
	if(kint.eq.37.or.kint.eq.41) wn3besy=-wn3*2*wn3*green_rearth*by2/(wn3*x)
	if(kint.eq.38) wn3besy=+wn3*4.*(wn3*green_rearth)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
	if(kint.eq.40) wn3besy=-wn3*2*(wn3*green_rearth)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
	if(kint.eq.42) wn3besy=-wn3*4.*(wn3*green_rearth)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
c
c  perform integration for each frequency
c  Use trapezoidal rule if wn*x < 3 and Filon else
c
	nwe=0
	do if=green_nf1,green_nf2
		nwn=green_nwn(if)
		call green_taper(if,tapfrac,greentap)
		j33=min(j3,nwn)
		zdis(if)=0.
		do j=2,j33-1
			zdis(if)=zdis(if)+greenkfsp(nwe+j)*wnbesj(j)*greentap(j)
		enddo
		zdis(if)=zdis(if)+0.5*greenkfsp(nwe+j33)*wnbesj(j33)*greentap(j33)
		zdis(if)=zdis(if)*green_dwn
		if(j33.eq.nwn) goto 11
c
c  use the Filon integration for kx > 3
c  Neumann terms at the end of integration interval (k_N)
c
		wn=(green_nwn(if)-1)*green_dwn
		by0=bessy0(wn*x)
		by1=bessy1(wn*x)
		by2=-by0+2./(wn*x)*by1
		if(kint.eq.1.or.kint.eq.15)  wnbesy=wn*by0
		if(kint.eq.2.or.kint.eq.16)  wnbesy=wn*by0
		if(kint.eq.3.or.kint.eq.17)  wnbesy=2./green_rearth*by1
		if(kint.eq.4.or.kint.eq.18)  wnbesy=wn*2.*by2
		if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wnbesy=wn*(-wn*green_rearth)*by1
		if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wnbesy=wn*(-wn*green_rearth)*by1
		if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wnbesy=wn*2.*(by0-by1/(wn*x))
		if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wnbesy=wn*2.*wn*green_rearth*(by1-2.*by2/(wn*x))
		if(kint.eq.9.or.kint.eq.23)  wnbesy=wn*2.*by1/(wn*x)
		if(kint.eq.10.or.kint.eq.24) wnbesy=wn*4.*wn*green_rearth*by2/(wn*x)
		if(kint.eq.11.or.kint.eq.25) wnbesy=wn*2.*(by0-by1/(wn*x))
		if(kint.eq.12.or.kint.eq.26) wnbesy=-wn*2.*wn*green_rearth*(by1-2.*by2/(wn*x))
		if(kint.eq.13.or.kint.eq.27) wnbesy=wn*2.*by1/(wn*x)
		if(kint.eq.14.or.kint.eq.28) wnbesy=-wn*4.*wn*green_rearth*by2/(wn*x)
		if(kint.eq.33.or.kint.eq.34) wnbesy=-wn*(wn*green_rearth)**2*(by0-by1/(wn*x))
		if(kint.eq.35.or.kint.eq.39) wnbesy=-wn*2*wn*green_rearth*(by1-by2/(wn*x))
		if(kint.eq.36) wnbesy=+wn*2*(wn*green_rearth)**2*((6./(wn*x)**2-1.)*by2-by1/(wn*x))
		if(kint.eq.37.or.kint.eq.41) wnbesy=-wn*2*wn*green_rearth*by2/(wn*x)
		if(kint.eq.38) wnbesy=+wn*4.*(wn*green_rearth)**2*(-3.*by2/(wn*x)**2+by1/(wn*x))
		if(kint.eq.40) wnbesy=-wn*2*(wn*green_rearth)**2*((6./(wn*x)**2-1.)*by2-by1/(wn*x))
		if(kint.eq.42) wnbesy=-wn*4.*(wn*green_rearth)**2*(-3.*by2/(wn*x)**2+by1/(wn*x))
		xdk=x*green_dwn
		zsum=greenkfsp(nwe+j3)*wnbesj(j3)*greentap(j3)
		do j=j3+1,nwn-1
			zsum=zsum+2.*greenkfsp(nwe+j)*wnbesj(j)*greentap(j)
		enddo
		zsum=zsum+greenkfsp(nwe+nwn)*wnbesj(nwn)*greentap(nwn)
		zdis(if)=zdis(if)+zsum/x*(1-cos(xdk))/xdk
		zdis(if)=zdis(if)+(greenkfsp(nwe+nwn)*wnbesy*greentap(nwn)
     1		         -greenkfsp(nwe+j3)*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
c
c  multiply by R^2/(2*pi)
c
 11		zdis(if)=zdis(if)*green_rearth**2/(2.*pi)
c
c  convert to nanometer per Nm
c
		zdis(if)=zdis(if)*1.e-6
c
c  update address pf last wavenumber-spectrum value
c
		nwe=nwe+nwn
	enddo
	return
	end
c--------------------------------------------------------------------------
c  Precompute values of Bessel functions at wavenumber points
c
	subroutine green_bessel(x,besselj)
	include 'greendim.h'
	include 'green.h'
	real x,besselj(npp,0:2),wn
	integer iwn
	real bessj0,bessj1
	external bessj0,bessj1
c
	wn=green_dwn
	do iwn=2,green_nwnmax
		besselj(iwn,0)=bessj0(wn*x)
		besselj(iwn,1)=bessj1(wn*x)
		besselj(iwn,2)=-besselj(iwn,0)+2./(wn*x)*besselj(iwn,1)
		wn=wn+green_dwn
	enddo
	return
	end
c--------------------------------------------------------------
c  Transform displacement vector and its gradient from epicentral
c  spherical to cartesian global coorinate system
c
c  dis: epicentral distance of scattering point S in km
c  azi: azimuth (south over east) of scattering point S as viewed from Q or E in radians
c  tm: transformation matrix M (see born-scattering.tex)
c
c  overwrites zsp with cartesian components
c
	subroutine green_trans_epispher_globalcar(r,dis,azi,zsp,csys)
	include 'greendim.h'
	include 'green.h'
	integer j,l,ir,is,if
	character*1 csys
	double precision r
	real dis,azi,tm(3,3),delta,rtandel
	complex zsp(nff,12),ucar(3),graduspher(3,3),work(3,3),graducar(3,3)
c
	if(csys.eq.'S') then
		delta=dis/green_rearth
		rtandel=r*tan(delta)
	endif
c
c  sar_tm_epispher_globalcar knows about csys and
c  provides the appropriate matrix
c
	call sar_tm_epispher_globalcar(delta,azi,tm)
c
	do if=green_nf1,green_nf2
c
c  transformation of displacement: U_j = M_js u_s 
c
		do j=1,3
			ucar(j)=0.
			do is=1,3
				ucar(j)=ucar(j)+tm(j,is)*zsp(if,is)
			enddo
		enddo
c
c  spherical gradient of displacement D_r u_s
c
		if(csys.eq.'S') then
			graduspher(1,1)=zsp(if,4)
			graduspher(1,2)=zsp(if,5)
			graduspher(1,3)=zsp(if,6)
			graduspher(2,1)=zsp(if,7)-zsp(if,2)/r
			graduspher(2,2)=zsp(if,8)+zsp(if,1)/r
			graduspher(2,3)=zsp(if,9)
			graduspher(3,1)=zsp(if,10)-zsp(if,3)/r
			graduspher(3,2)=zsp(if,11)-zsp(if,3)/rtandel
			graduspher(3,3)=zsp(if,12)+zsp(if,1)/r+zsp(if,2)/rtandel
		else
c
c  spherical gradient of displacement D_r u_s in shallow seismic approximation
c  spherical components are essentially treated as cylindrical ones
c  r -> inf, delta -> 0, r sin(delta) -> dis, cot(delta)/r -> 1/dis
c
			graduspher(1,1)=zsp(if,4)
			graduspher(1,2)=zsp(if,5)
			graduspher(1,3)=zsp(if,6)
			graduspher(2,1)=zsp(if,7)
			graduspher(2,2)=zsp(if,8)
			graduspher(2,3)=zsp(if,9)
			graduspher(3,1)=zsp(if,10)
			graduspher(3,2)=zsp(if,11)-zsp(if,3)/dis
			graduspher(3,3)=zsp(if,12)+zsp(if,2)/dis
		endif
c
c  now overwrite displacement components
c
		do j=1,3
			zsp(if,j)=ucar(j)
		enddo
c
c  transformation of gradient d_j u_l = M_jr D_rs M_ls
c  first calculate D_rs M_ls
c
		do l=1,3
			do ir=1,3
				work(ir,l)=0.
				do is=1,3
					work(ir,l)=work(ir,l)+graduspher(ir,is)*tm(l,is)
				enddo
			enddo
		enddo
c
c  now calculate d_j u_l = M_jr * work(ir,l)
c
		do j=1,3
			do l=1,3
				graducar(j,l)=0.
				do ir=1,3
					graducar(j,l)=graducar(j,l)+tm(j,ir)*work(ir,l)
				enddo
			enddo
		enddo
c
c overwrite spherical gradient by cartesian gradient
c e.g. graducar(2,3) = d_y u_z
c
		zsp(if,4)=graducar(1,1)
		zsp(if,5)=graducar(1,2)
		zsp(if,6)=graducar(1,3)
		zsp(if,7)=graducar(2,1)
		zsp(if,8)=graducar(2,2)
		zsp(if,9)=graducar(2,3)
		zsp(if,10)=graducar(3,1)
		zsp(if,11)=graducar(3,2)
		zsp(if,12)=graducar(3,3)
	enddo
c
	return
	end
