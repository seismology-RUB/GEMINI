c--------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
c---------------------------------------------------------------
c  Overview of routines
c--------------------------------------------------------------
c      subroutine earthmodel_fromflnm(f,ir)
c      subroutine earthmodel_spline
c      subroutine earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c      subroutine earthmodel_q(r,qk,qm)
c      subroutine earthmodel_rholiq(r,ro,rop)
c      subroutine earthmodel_isovel(nl,r,alf,bet)
c      subroutine earthmodel_mindr(nl,drmin)
c      subroutine earthmodel_isosample(dr,n,rem,ro,vp,vs,qk,qm)
c      subroutine earthmodel_homlay(dr,n,rem,ro,vp,vs,qk,qm)
c      subroutine earthmodel_getfsflag(flag)
c      subroutine earthmodel_getrearth(r)
c      subroutine earthmodel_getuppernodeindex(nl,r,j)
c----------------------------------------------------------------
c  Create earth model from flnm type
c
c  change to a table of complex elastic constants at nodes
c
c  30/11/04: treat full-space models (earthmodel_fsflag: full-space flag)
c----------------------------------------------------------------
c  m4 type specification for elastic constants
c
	define(m4_elcon_type,`double complex')
c-----------------------------------------------------------------
	subroutine earthmodel_fromflnm(f,ir)
	include 'flnm.h'
	include 'earthmodel.h'
	integer j,ir
	double precision fratio,f
c
	if(flnm_seldamp.eq.1) then 
		fratio=1.d0
	else if(flnm_seldamp.eq.3) then
		fratio=f/flnm_fref
c		print *,f,flnm_fref,fratio
	else 
		print *,' Sorry, invalid or unimplemented value of seldamp'
		stop
	endif
	earthmodel_nk=flnm_nk
	earthmodel_seldamp=flnm_seldamp
	earthmodel_aniflag=flnm_aniflag
	earthmodel_fref=flnm_fref
	earthmodel_rearth=flnm_rearth
	earthmodel_fsflag=flnm_fsflag
c
c  for full-space models subtract 1 from earthmodel_nk
c  but store parameters of upper sphere at index earthmodel_nk+1 !!
c
	if(earthmodel_fsflag.eq.1) earthmodel_nk=flnm_nk-1
	if(ir.eq.0) print *,'nk = ',earthmodel_nk, ' fsflag = ',earthmodel_fsflag
c
	if(ir.eq.0) print *,'    Earthmodel nodes  '
	do j=1,flnm_nk
		earthmodel_rk(j)=flnm_rk(j)
		earthmodel_ro(j)=flnm_ro(j)
		earthmodel_qk(j)=flnm_qk(j)
		earthmodel_qm(j)=flnm_qm(j)
		call zelcon(fratio,flnm_ro(j),flnm_vpv(j),flnm_vph(j),flnm_vsv(j),flnm_vsh(j),
     1		            flnm_eta(j),flnm_qk(j),flnm_qm(j),
     1		            earthmodel_za(j),earthmodel_zc(j),earthmodel_zf(j),
     1		            earthmodel_zl(j),earthmodel_zn(j),
     1		            earthmodel_zkap(j),earthmodel_zmue(j))
		if(ir.eq.0) write(6,'(i5,2f12.3)') j,earthmodel_rk(j),(earthmodel_rearth-earthmodel_rk(j))*1.d3
	enddo
	return
	end
c----------------------------------------------------------------
c  Spline model parameters between discontinuities
c----------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_splco,ifelse(m4_elcon_type,`double complex',`zsplco($*)',`dcubsplco($1,$2,$3,$4,$7)'))
c---------------------------------------------------------------
	subroutine earthmodel_spline
	include 'earthmodel.h'
	include 'layer.h'
	integer nl,j1,j2,n
	double precision y(nkk),y2(nkk),work(nkk,3)
c
c  do not spline halfspace
c
	do nl=2,layer_nlay
		j1=layer_iktop(nl-1)+1
		j2=layer_iktop(nl)
		n=j2-j1+1
c
c  spline elastic parameters
c
		call dcubsplco(earthmodel_rk(j1),earthmodel_ro(j1),n,earthmodel_ro2(j1),work)
		call dcubsplco(earthmodel_rk(j1),earthmodel_qk(j1),n,earthmodel_qk2(j1),work)
		call dcubsplco(earthmodel_rk(j1),earthmodel_qm(j1),n,earthmodel_qm2(j1),work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zkap(j1),n,earthmodel_zkap2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zmue(j1),n,earthmodel_zmue2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_za(j1),n,earthmodel_za2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zc(j1),n,earthmodel_zc2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zf(j1),n,earthmodel_zf2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zl(j1),n,earthmodel_zl2(j1),y,y2,work)
		call m4_splco(earthmodel_rk(j1),earthmodel_zn(j1),n,earthmodel_zn2(j1),y,y2,work)
	enddo
	return
	end
c----------------------------------------------------------------
c  Evaluate model for given radius and layer index
c----------------------------------------------------------------
	subroutine earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
	include 'earthmodel.h'
	include 'layer.h'
	double precision r,ro,a,b,c,d,h,x1,x2
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zkap,zmu
	integer nl,j1,j2,n,j
c
	nl=layer_nlactive
c
c  constant parameters in bottom homomgeneous sphere
c
	if(nl.eq.1) then
		j=layer_iktop(1)
		ro=earthmodel_ro(j)
		zpa=earthmodel_za(j)
		zpc=earthmodel_zc(j)
		zpf=earthmodel_zf(j)
		zpl=earthmodel_zl(j)
		zpn=earthmodel_zn(j)
		zkap=earthmodel_zkap(j)
		zmu=earthmodel_zmue(j)
		return
	endif
c
c  constant parameters in outer homomgeneous sphere (if any)
c  stored at index earthmodel_nk+1 (see earthmodel_fromflnm)
c  layer index of upper sphere is assumed to be nlay+1
c
	if(earthmodel_fsflag.eq.1.and.nl.eq.layer_nlay+1) then
		j=earthmodel_nk+1
		ro=earthmodel_ro(j)
		zpa=earthmodel_za(j)
		zpc=earthmodel_zc(j)
		zpf=earthmodel_zf(j)
		zpl=earthmodel_zl(j)
		zpn=earthmodel_zn(j)
		zkap=earthmodel_zkap(j)
		zmu=earthmodel_zmue(j)
		return
	endif
c
c  above halfspace
c
	j1=layer_iktop(nl-1)+1
	j2=layer_iktop(nl)
	n=j2-j1+1
	call dlocate(r,n,earthmodel_rk(j1),j)
	if(j.eq.0) j=1
	j=j+j1-1
	x1=earthmodel_rk(j)
	x2=earthmodel_rk(j+1)
	h=x2-x1
	a=(x2-r)/h
	b=1.d0-a
	c=h**2*a*(a**2-1.d0)/6.d0
	d=h**2*b*(b**2-1.d0)/6.d0
	ro=a*earthmodel_ro(j)+b*earthmodel_ro(j+1)+c*earthmodel_ro2(j)+d*earthmodel_ro2(j+1)
	zpa=a*earthmodel_za(j)+b*earthmodel_za(j+1)+c*earthmodel_za2(j)+d*earthmodel_za2(j+1)
	zpc=a*earthmodel_zc(j)+b*earthmodel_zc(j+1)+c*earthmodel_zc2(j)+d*earthmodel_zc2(j+1)
	zpf=a*earthmodel_zf(j)+b*earthmodel_zf(j+1)+c*earthmodel_zf2(j)+d*earthmodel_zf2(j+1)
	zpl=a*earthmodel_zl(j)+b*earthmodel_zl(j+1)+c*earthmodel_zl2(j)+d*earthmodel_zl2(j+1)
	zpn=a*earthmodel_zn(j)+b*earthmodel_zn(j+1)+c*earthmodel_zn2(j)+d*earthmodel_zn2(j+1)
	zkap=a*earthmodel_zkap(j)+b*earthmodel_zkap(j+1)+c*earthmodel_zkap2(j)+d*earthmodel_zmue2(j+1)
	zmu= a*earthmodel_zmue(j)+b*earthmodel_zmue(j+1)+c*earthmodel_zmue2(j)+d*earthmodel_zmue2(j+1)
	return
	end
c---------------------------------------------------------------------
c  Get values for qk and qm 
c
	subroutine earthmodel_q(r,qk,qm)
	include 'earthmodel.h'
	include 'layer.h'
	double precision r,qk,qm,a,b,c,d,h,x1,x2
	integer nl,j1,j2,n,j
c
	nl=layer_nlactive
c
c  constant parameters in bottom homomgeneous sphere
c
	if(nl.eq.1) then
		j=layer_iktop(1)
		qk=earthmodel_qk(j)
		qm=earthmodel_qm(j)
		return
	endif
c
c  constant parameters in outer homomgeneous sphere (if any)
c  stored at index earthmodel_nk+1 (see earthmodel_fromflnm)
c  layer index of upper sphere is assumed to be nlay+1
c
	if(earthmodel_fsflag.eq.1.and.nl.eq.layer_nlay+1) then
		j=earthmodel_nk+1
		qk=earthmodel_qk(j)
		qm=earthmodel_qm(j)
		return
	endif
c
c  above halfspace
c
	j1=layer_iktop(nl-1)+1
	j2=layer_iktop(nl)
	n=j2-j1+1
	call dlocate(r,n,earthmodel_rk(j1),j)
	if(j.eq.0) j=1
	j=j+j1-1
	x1=earthmodel_rk(j)
	x2=earthmodel_rk(j+1)
	h=x2-x1
	a=(x2-r)/h
	b=1.d0-a
	c=h**2*a*(a**2-1.d0)/6.d0
	d=h**2*b*(b**2-1.d0)/6.d0
	qk=a*earthmodel_qk(j)+b*earthmodel_qk(j+1)+c*earthmodel_qk2(j)+d*earthmodel_qk2(j+1)
	qm=a*earthmodel_qm(j)+b*earthmodel_qm(j+1)+c*earthmodel_qm2(j)+d*earthmodel_qm2(j+1)
	return
	end
c--------------------------------------------------------------------
c  Get density and derivative in liquid layer
c
	subroutine earthmodel_rholiq(r,ro,rop)
	include 'earthmodel.h'
	include 'layer.h'
	double precision r,ro,rop,a,ap,b,bp,c,cp,d,dp,h,x1,x2
	integer nl,j1,j2,n,j
c
	nl=layer_nlactive
	if(layer_iflso(nl).eq.0) then
		print *,'You called earthmodel_rholiq for a solid layer!!'
		stop
	endif
c
c  constant parameters in bottom homomgeneous sphere
c
	if(nl.eq.1) then
		j=layer_iktop(1)
		ro=earthmodel_ro(j)
		rop = 0.d0
		return
	endif
c
c  constant parameters in outer homomgeneous sphere (if any)
c  stored at index earthmodel_nk+1 (see earthmodel_fromflnm)
c  layer index of upper sphere is assumed to be nlay+1
c
	if(earthmodel_fsflag.eq.1.and.nl.eq.layer_nlay+1) then
		j=earthmodel_nk+1
		ro=earthmodel_ro(j)
		rop=0.d0
		return
	endif		
c
c  above halfspace
c
	j1=layer_iktop(nl-1)+1
	j2=layer_iktop(nl)
	n=j2-j1+1
	call dlocate(r,n,earthmodel_rk(j1),j)
	if(j.eq.0) j=1
	j=j+j1-1
	x1=earthmodel_rk(j)
	x2=earthmodel_rk(j+1)
	h=x2-x1
	a=(x2-r)/h
	ap=-1.d0/h
	b=1.d0-a
	bp=-ap
	c=h**2*a*(a**2-1.d0)/6.d0
	cp=h**2*(3.d0*a**2-1.d0)/6.d0*ap
	d=h**2*b*(b**2-1.d0)/6.d0
	dp=h**2*(3.d0*b**2-1.d0)/6.d0*bp
	ro=a*earthmodel_ro(j)+b*earthmodel_ro(j+1)+c*earthmodel_ro2(j)+d*earthmodel_ro2(j+1)
	rop=ap*earthmodel_ro(j)+bp*earthmodel_ro(j+1)+cp*earthmodel_ro2(j)+dp*earthmodel_ro2(j+1)
	return
	end
c----------------------------------------------------------------
c  Get real isotropic velocity values for given radius in layer nl
c
	subroutine earthmodel_isovel(nl,r,alf,bet)
	include 'earthmodel.h'
	include 'layer.h'
	double precision r,alf,bet,ro,a,b,c,d,h,x1,x2
	m4_elcon_type zkap,zmu
	integer j1,j2,nl,n,j
c
c  constant parameters in halfspace
c
	if(nl.eq.1) then
		j=layer_iktop(1)
		ro=earthmodel_ro(j)
		zkap=earthmodel_zkap(j)
		zmu=earthmodel_zmue(j)
	else
c
c  above halfspace
c
		j1=layer_iktop(nl-1)+1
		j2=layer_iktop(nl)
		n=j2-j1+1
		call dlocate(r,n,earthmodel_rk(j1),j)
		if(j.eq.0) j=1
		j=j+j1-1
		x1=earthmodel_rk(j)
		x2=earthmodel_rk(j+1)
		h=x2-x1
		a=(x2-r)/h
		b=1.d0-a
		c=h**2*a*(a**2-1.d0)/6.d0
		d=h**2*b*(b**2-1.d0)/6.d0
		ro=a*earthmodel_ro(j)+b*earthmodel_ro(j+1)+c*earthmodel_ro2(j)+d*earthmodel_ro2(j+1)
		zkap=a*earthmodel_zkap(j)+b*earthmodel_zkap(j+1)+c*earthmodel_zkap2(j)+d*earthmodel_zmue2(j+1)
		zmu= a*earthmodel_zmue(j)+b*earthmodel_zmue(j+1)+c*earthmodel_zmue2(j)+d*earthmodel_zmue2(j+1)
	endif
c
c  velocities
c
	alf=dsqrt( real(zkap+4.d0*zmu/3.d0)/ro )
	bet=dsqrt( real(zmu)/ro )
c
	return
	end
c----------------------------------------------------------------------
c   Find minimum distance between earthmodel nodes in a layer
c
	subroutine earthmodel_mindr(nl,drmin)
	include 'earthmodel.h'
	include 'layer.h'
	integer nl,i
	double precision drmin
c
	if(nl.eq.1) then
		drmin=1.d12
		return
	endif
	drmin=1.d12
	do i=layer_iktop(nl-1)+2,layer_iktop(nl)
		 drmin=min(drmin,earthmodel_rk(i)-earthmodel_rk(i-1))
	enddo
	return
	end
c------------------------------------------------------------------------
c   Sample earthmodel with defined dr for use with plotting
c   Double precision in- and output
c   qk and qs are true Q's not 1/Q
c
	subroutine earthmodel_isosample(dr,n,rem,ro,vp,vs,qk,qm)
	include 'earthmodel.h'
	include 'layer.h'
	integer n,i,nl,nval
	double precision dr,rem(n),ro(n),vp(n),vs(n),qk(n),qm(n)
	double precision r,drs,h,rho,qkh,qmh
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zkap,zmu
c
	n=1
	rem(n)=earthmodel_rk(1)
	layer_nlactive=1
	call earthmodel_elpar(earthmodel_rk(1),rho,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
	ro(1)=rho
	vp(1)=dsqrt( real(zkap+4.d0*zmu/3.d0)/rho )
	vs(1)=dsqrt( real(zmu)/rho )
	qk(1)=earthmodel_qk(1)
	if (layer_iflso(nl).eq.0) then
		qm(1)=earthmodel_qm(1)
	else
		qm(1)=-1.
	endif
	do nl=2,layer_nlay
		h=layer_rb(nl)-layer_rb(nl-1)
		nval=h/dr
		drs=h/nval
		layer_nlactive=nl
		do i=1,nval+1
			r=layer_rb(nl-1)+(i-1)*drs
			call earthmodel_elpar(r,rho,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			call earthmodel_q(r,qkh,qmh)
			n=n+1
			rem(n)=r
			ro(n)=rho
			vp(n)=dsqrt( real(zkap+4.d0*zmu/3.d0)/rho )
			vs(n)=dsqrt( real(zmu)/rho )
			qk(n)=qkh
			qm(n)=qmh
			if (layer_iflso(nl).eq.1) qm(n)=-1.
		enddo
	enddo
	return
	end
c-------------------------------------------------------------------------------
c  create an earth model consisting of a stack of homogeneous layers
c  Layer i: bottom radius rem(i-1), top radius rem(i), velocity in layer v(i)
c
	subroutine earthmodel_homlay(dr,n,rem,ro,vp,vs,qk,qm)
	include 'earthmodel.h'
	include 'layer.h'
	integer n,i,nl,nval
	double precision dr,rem(0:n),ro(n),vp(n),vs(n),qk(n),qm(n)
	double precision r,drs,h,rho,qkh,qmh
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zkap,zmu
c
	n=0
	rem(0)=layer_rb(1)
	do nl=2,layer_nlay
		h=layer_rb(nl)-layer_rb(nl-1)
		nval=h/dr
		drs=h/nval
		layer_nlactive=nl
		do i=1,nval
			r=layer_rb(nl-1)+i*drs
			call earthmodel_elpar(r-0.5*drs,rho,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			call earthmodel_q(r-0.5*drs,qkh,qmh)
			n=n+1
			rem(n)=r
			ro(n)=rho
			vp(n)=dsqrt( real(zkap+4.d0*zmu/3.d0)/rho )
			vs(n)=dsqrt( real(zmu)/rho )
			qk(n)=qkh
			qm(n)=qmh
			if (layer_iflso(nl).eq.1) qm(n)=-1.
		enddo
	enddo
	return
	end

c---------------------------------------------------------------------------------
c  return value of earthmodel_fsflag
c
	subroutine earthmodel_getfsflag(flag)
	integer flag
	include 'earthmodel.h'
c
	flag = earthmodel_fsflag
c
	return
	end
c---------------------------------------------------------------------------------
c  return value of earthmodel_rearth
c
	subroutine earthmodel_getrearth(r)
	double precision r
	include 'earthmodel.h'
c
	r = earthmodel_rearth
c
	return
	end
c------------------------------------------------
c  return index of model node above given radius
c  located in layer nl
c
	subroutine earthmodel_getuppernodeindex(nl,r,j)
	include 'earthmodel.h'
	include 'layer.h'
	double precision r
	integer j,nl,j1,j2,n
c
	if(nl.eq.1) then
		j=layer_iktop(1)
	else
		j1=layer_iktop(nl-1)+1
		j2=layer_iktop(nl)
		n=j2-j1+1
		call dlocate(r,n,earthmodel_rk(j1),j)
		if(j.eq.0) j=1
		j=j+j1
	endif
c
	return
	end
c---------------------------------------------------
c  get earthmodel node radius for given index
c
	subroutine earthmodel_getrk_selected(k,r)
	include 'earthmodel.h'
	integer k
	double precision r
	r = earthmodel_rk(k)
	return
	end
c---------------------------------------------------
c  get number of earthmodel nodes
c
	subroutine earthmodel_getnk(k)
	include 'earthmodel.h'
	integer k
	k = earthmodel_nk
	return
	end
