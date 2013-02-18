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
c                 STERMS
c
c  evaluate source terms (jump of Green 6-vector at source)
c  for either a moment tensor or single force source.
c
c  See convention for source terms and expansion coefficients
c  in gemini-dsipl.tex.
c
c  sourcetype: f = single force
c              m = moment tensor    
c  zs:         depth of source
c  zsph:       jump of Green 4-vector at source (output)
c              second index is used in case of several sources
c  ztor:       jump of Green 2-vector at source (output)
c              second index is used in case of several sources
c---------------------------------------------------------------
	subroutine sterms(sourcetype,zs,zsph,ztor)
	double complex zsph(4,4),ztor(2,2),zpa,zpc,zpf,zpl,zpn,zkap,zmu,ci
	integer i,j,nls,top,nlay
	double precision r,rr2,rr3,ro,zs,rearth
	character*1 sourcetype
c
	ci=dcmplx(0.d0,1.d0)
c
	call layer_getnlay(nlay)
	call layer_getrb(nlay,rearth)
	r=rearth-zs
	call layer_getindex(r,nls,top)
c
	rr2=1./(r*r)
	rr3=rr2/r
	do j=1,4
	do i=1,4
		zsph(i,j)=0.d0
	enddo
	enddo
	do j=1,2
	do i=1,2
		ztor(i,j)=0.d0
	enddo
	enddo
c
	call layer_setnlactive(nls)
	call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c
c  general force (fr, ftheta, fphi)
c
	if(sourcetype.eq.'f') then
		zsph(2,1)=-rr2
		zsph(4,2)=0.5*rr2
		ztor(2,1)=0.5*rr2
c
c  general moment tensor
c
	else if(sourcetype.eq.'m') then
c
c  m=0, basis solution for Mrr
c
		zsph(1,1)=rr2/zpc
		zsph(2,1)=rr3*2.d0*zpf/zpc
		zsph(3,1)=0.d0
		zsph(4,1)=-rr3*zpf/zpc
c
c  m=0, basis solution for (Mtt+Mff)
c
		zsph(1,2)=0.d0
		zsph(2,2)=-rr3
		zsph(3,2)=0.d0
		zsph(4,2)=0.5*rr3
c
c  m=+-1, basis solution for +-Mrt + iMrf, sqrt(1/[l(l+1)]) not included!
c
		zsph(1,3)=0.d0
		zsph(2,3)=0.d0
		zsph(3,3)=0.5*rr2/zpl
		zsph(4,3)=0.d0
c
c  m=+-2, basis solution for Mff-Mtt +- 2iMtf, sqrt[(l+2)(l-1)]/[l(l+1)]] is not included
c
		zsph(1,4)=0.d0
		zsph(2,4)=0.d0
		zsph(3,4)=0.d0
		zsph(4,4)=0.25*rr3
c
c  m=+-1, basis solution for +-Mrf + iMrt, sqrt(1/[l(l+1)]) not included!
c
		ztor(1,1)=0.5*rr2/zpl
		ztor(2,1)=0.d0
c
c  m=+-2, basis solution for 2Mtf -+ i(Mff-Mtt), sqrt[(l+2)(l-1)]/[l(l+1)]] is not included
c
		ztor(1,2)=0.d0
		ztor(2,2)=0.25*rr3
	endif
c
	return
	end
