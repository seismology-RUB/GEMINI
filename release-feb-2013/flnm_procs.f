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
c  Procedures for class flnm
c
c  30/11/04: treat full-space models (flnm_fsflag: full-space flag)
c            consequences for flinv1d.f unchecked!
c----------------------------------------------------------
c
c  Read model file and fill common/flnm
c
	subroutine flnm_read(iunit,cmodel)
	include 'flnm.h'
	double precision z(nkk1),ro(nkk1),vp(nkk1),vs(nkk1),qk(nkk1),qm(nkk1),f
	double precision vpv(nkk1),vph(nkk1),vsv(nkk1),vsh(nkk1),eta(nkk1),rearth
	integer i,j,iunit,dampflag,nk,aniflag
	character text*80,cmodel*(*)
c
	open(iunit,file=cmodel,status='old')
c
c  read first line of model file
c
	read(iunit,'(a)') text
	if(text.ne.'FLNODE') then
		print *,'<flnm_read>: Wrong model type. Should be FLNODE !'
		stop
	endif
c
c  reading of node-File
c
	read(iunit,'(a)') text
	read(iunit,*) f,dampflag,aniflag,rearth
	read(iunit,*) nk
	if(nk.gt.nkk1) then
		print *,'<flnm_readiso>: Too many model nodes!'
		print *,'actual: ',nk,' dimensioned: ',nkk1
		stop
	endif
	if(aniflag.eq.0) then
		do i=1,nk
			j=nk+1-i
			read(iunit,*) z(j),ro(j),vp(j),vs(j),qk(j),qm(j)
		enddo
		call flnm_isofromarrays(nk,dampflag,f,rearth,z,ro,vp,vs,qk,qm)
	else
		do i=1,nk
			j=nk+1-i
			read(iunit,*) z(j),ro(j),vpv(j),vph(j),vsv(j),vsh(j),eta(j),qk(j),qm(j)
		enddo
		call flnm_anisofromarrays(nk,dampflag,f,rearth,z,ro,vpv,vph,vsv,vsh,eta,qk,qm)
	endif
	close(iunit)
	return
	end
c-----------------------------------------------------------------
c  fill flnm_array from arrays in the argument list
c  isotropic model
c
c
c  nk:         number of model nodes
c  dampflag:   switch to select damping convention
c  f:          reference frequency fref
c  z:          node depth in km
c  ro:         density
c  vp:         P-velocity
c  vs:         S-Velocity
c  qk:         Qkappa (not 1/Q)
c  qm:         Qmu
c
	subroutine flnm_isofromarrays(nk,dampflag,f,rearth,z,ro,vp,vs,qk,qm)
	include 'flnm.h'
	double precision qk(nk),qm(nk),z(nk),ro(nkk1),vp(nkk1),vs(nkk1),f,rearth
	integer i,nk,dampflag
c
	if(nk.gt.nkk1) then
		print *,'<flnm_isofromarrays>: Too many model nodes!'
		print *,'actual: ',nk,' dimensioned: ',nkk1
		stop
	endif
	flnm_nk=nk
	flnm_fref=f
	flnm_seldamp=dampflag
	flnm_aniflag=0
	flnm_rearth=rearth
	flnm_fsflag=0
c
c  check for full-space (z=0.0 occurs twice!)
c  and set flnm_fsflag
c
	if(dabs(z(nk)-z(nk-1)).lt.1.e-6) then
		flnm_fsflag=1
	endif
	do i=1,nk
		flnm_rk(i)=rearth-z(i)*1.d-3
		flnm_ro(i)=ro(i)
		flnm_vpv(i)=vp(i)
		flnm_vph(i)=vp(i)
		flnm_vsv(i)=vs(i)
		flnm_vsh(i)=vs(i)
		flnm_eta(i)=1.d0
		flnm_qk(i)=1./qk(i)
		if(qm(i).lt.0.) then
			flnm_qm(i)=0.d0
		else
			flnm_qm(i)=1./qm(i)
		endif
	enddo
	return
	end
c-----------------------------------------------------------------
c  fill flnm_array from arrays in the argument list
c  transversely isotropic model
c
c
c  nk:         number of model nodes
c  dampflag:   switch to select damping convention
c  f:          reference frequency fref
c  z:          node depth in km
c  ro:         density
c  vpv,vph:    P-velocities
c  vsv,vsh:    S-Velocities
c  eta:        F/(A-2L)
c  qk:         Qkappa (not 1/Q)
c  qm:         Qmu
c
	subroutine flnm_anisofromarrays(nk,dampflag,f,rearth,z,ro,vpv,vph,vsv,vsh,eta,qk,qm)
	include 'flnm.h'
	double precision qk(nk),qm(nk),z(nk),ro(nkk1),vpv(nkk1),vph(nkk1)
	double precision vsv(nkk1),vsh(nkk1),eta(nkk1),f,rearth
	integer i,nk,dampflag
c
	if(nk.gt.nkk1) then
		print *,'<flnm_isofromarrays>: Too many model nodes!'
		print *,'actual: ',nk,' dimensioned: ',nkk1
		stop
	endif
	flnm_nk=nk
	flnm_fref=f
	flnm_seldamp=dampflag
	flnm_aniflag=1
	flnm_rearth=rearth
	flnm_fsflag=0
c
c  check for full-space (z=0.0 occurs twice!)
c  and set flnm_fsflag
c
	if(dabs(z(nk)-z(nk-1)).lt.1.e-6) then
		flnm_fsflag=1
	endif
	do i=1,nk
		flnm_rk(i)=rearth-z(i)*1.e-3
		flnm_ro(i)=ro(i)
		flnm_vpv(i)=vpv(i)
		flnm_vph(i)=vph(i)
		flnm_vsv(i)=vsv(i)
		flnm_vsh(i)=vsh(i)
		flnm_eta(i)=eta(i)
		flnm_qk(i)=1./qk(i)
		if(qm(i).lt.0.) then
			flnm_qm(i)=0.d0
		else
			flnm_qm(i)=1./qm(i)
		endif
	enddo
	return
	end
c------------------------------------------------------------------
c   Read model file but do not fill common but return parameters
c   via the argument list. Single precision!
c
	subroutine flnm_readiso_sngl(iunit,cmodel,dampflag,fref,rearth,nparr,nkk,nk,zem,parfun)
	integer iunit,nk,nparr,nkk
	integer i,j,dampflag,aniflag,ipar
	real zem(nk),parfun(nkk,nparr),rearth,fref
	character text*80,cmodel*(*)
c
	open(iunit,file=cmodel,status='old')
c
c  read first line of model file
c
	read(iunit,'(a)') text
	if(text.ne.'FLNODE') then
		print *,'<flnm_read>: Wrong model type. Should be FLNODE !'
		stop
	endif
c
c  reading of node-File
c
	read(iunit,'(a)') text
	read(iunit,*) fref,dampflag,aniflag,rearth
	read(iunit,*) nk
	if(nk.gt.nkk) then
		print *,'<flnm_readiso_sngl>: Too many model nodes!'
		print *,'actual: ',nk,' dimensioned: ',nkk
		stop
	endif
	if(aniflag.eq.1) then
		print *,'Not an isotropic model !'
		stop
	endif
	do i=1,nk
		j=nk+1-i
		read(iunit,*) zem(j),(parfun(j,ipar),ipar=1,nparr)
		parfun(j,4)=1./parfun(j,4)
		if(parfun(j,5).lt.0.) then
			parfun(j,5)=0.d0
		else
			parfun(j,5)=1./parfun(j,5)
		endif
	enddo
	close(iunit)
	return
	end
c------------------------------------------------------------------
c   Write model file. Single precision!
c
	subroutine flnm_writeiso_sngl(iunit,cmodel,dampflag,fref,rearth,nparr,
     1	                              nkk,nk,zem,parfun)
	integer iunit,nk,nparr,nkk,aniflag
	integer i,j,dampflag,ipar
	real zem(nk),parfun(nkk,nparr),rearth,fref,parqk,parqm
	character cmodel*(*)
	data aniflag/0/
c
	open(iunit,file=cmodel)
c
c  read first line of model file
c
	write(iunit,'(a)') 'FLNODE'
c
c  reading of node-File
c
	write(iunit,'(a)') 'Inverted model'
	write(iunit,*) fref,dampflag,aniflag,rearth
	write(iunit,*) nk
	if(nk.gt.nkk) then
		print *,'<flnm_writeiso_sngl>: Too many model nodes!'
		print *,'actual: ',nk,' dimensioned: ',nkk
		stop
	endif
	do i=1,nk
		j=nk+1-i
		parqk=1./parfun(j,4)
		if(parfun(j,3).le.0.) then
			parqm=-1.d0
		else
			parqm=1./parfun(j,5)
		endif
		write(iunit,'(f10.7,5f11.6)') zem(j),(parfun(j,ipar),ipar=1,3),parqk,parqm
	enddo
	close(iunit)
	return
	end
c---------------------------------------------------------------
c  get flnm_fref
c
	subroutine flnm_get_fref(fref)
	include 'flnm.h'
	double precision fref
	fref = flnm_fref
	return
	end
