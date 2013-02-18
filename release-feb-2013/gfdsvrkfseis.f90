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
!---------------------------------------------------------------------
!  GFDSVRKFSEIS
!---------------------------------------------------------------------
 program gfdsvrkfseis
	use mathConstants
	use commandLine
	use sourceReceiver
	use greenFrequencyWavenumber
	use greenSection
	use errorMessage
	use fileUnitHandler
!
	implicit none
	type (cmdLine) :: cl
	type (sar) :: info
	type (green_frequency_wavenumber) :: gfk
	type (green_section) :: greensec
	type (error_message) :: errmsg
	type (file_unit_handler) :: fuh
	real, dimension(:), allocatable :: stf
!	
	integer :: i,ifwin,ifqtf,lu
	integer :: nstf,nordh,nordl
	real :: fch,fcl,dt,dtp,dtstf,gv,tcut,zs,epimax,dum
	character :: string*132,seisfile*132,infofile*80
	character :: stfasc*80,dsvbasename*132,outformat*3
	character (len=132) :: myname = 'gfdsvrkfseis'
	logical :: vsflag,parflag,twodflag
	external printhelp
!
!  defaults for input parameters
!
	data seisfile/'greenseis.sff'/nordh/0/fch/10./nordl/0/fcl/100./epimax/100./ifwin/0/ &
	   & gv/2.0/tcut/0./dt/0.001/ifqtf/0/vsflag/.false./parflag/.false./twodflag/.false./ &
	   & outformat/'sff'/
!----------------------------------------------------------
!  Declarations for tf_cmdline
!  Use data to provide optid and hasarg
!
	call new(cl,13,2,'h o f hp lp e gvw dt qtf zs vs p 2d','0 1 1 1 1 1 1 1 1 1 0 0 0',printhelp)
	infofile = clManarg(cl,1)
	dsvbasename = clManarg(cl,2)
	if(clOptset(cl,'o')) then; seisfile = clOptarg(cl,'o'); endif
	if(clOptset(cl,'f')) then; outformat = clOptarg(cl,'f'); endif
	if(clOptset(cl,'hp')) then; string = clOptarg(cl,'hp'); read(string,*) nordh,fch; endif
	if(clOptset(cl,'lp')) then; string = clOptarg(cl,'lp'); read(string,*) nordl,fcl; endif
	if(clOptset(cl,'e')) then; string = clOptarg(cl,'e'); read(string,*) epimax; endif
	if(clOptset(cl,'gvw')) then; ifwin = 1; string = clOptarg(cl,'gvw'); read(string,*) gv,tcut; endif
	if(clOptset(cl,'dt')) then; string = clOptarg(cl,'dt'); read(string,*) dtp; endif
	if(clOptset(cl,'qtf')) then; ifqtf = 1; stfasc = clOptarg(cl,'qtf'); endif
	if(clOptset(cl,'zs')) then; string = clOptarg(cl,'zs'); read(string,*) zs; endif
	if(clOptset(cl,'vs')) then; vsflag = .true.; endif
	if(clOptset(cl,'p')) then; parflag = .true.; endif
!	if(clOptset(cl,'2d')) then; twodflag = .true.; endif
!-------------------------------------------------------
!  go
!
	call new(fuh,20)
!
!  Read info.eventid file
!  and convert station coordinates to epicentral
!  system -> dis(i), phi(i), dis(i) in km
!
	call new(info,1,infofile,epimax)
	if(sar_getNetFlag(info) .lt. 2) then
		print *,'You need an infofile with component specification'
		stop
	endif                                  !< take source depth from info file if not on command line
	if (.not. clOptset(cl,'zs')) then
		if (.csys.info /= 'S') then
			zs = .zs.info
		else
			zs = (.zs.info)*1.e3
		endif
		print *,'Source depth taken from info file: ',zs
	endif
!
!  open frequency-wavenumber spectrum file
!  and read in group tree and header info
!
	if (.not. parflag) then
		errmsg = readHeaderGreenFrequencyWavenumber(gfk,fuh,dsvbasename,vsflag)
	else
		errmsg = readHeaderParallelGreenFrequencyWavenumber(gfk,fuh,dsvbasename,vsflag)
	endif
	if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
	if (.styp.gfk /= .styp.info) then
		call new(errmsg,2,'Source type in frequency wavenumber file different from that in infofile',myname)
		call print(errmsg)
		stop
	endif
!
!  read fk-spectra for depth zs
!
	if (.not. parflag) then
		errmsg = readBodyGreenFrequencyWavenumber(gfk,zs,vsflag)
	else
		errmsg = readParallelBodyGreenFrequencyWavenumber(gfk,zs,vsflag)
	endif
	if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  convert from rad to km if spherical cs
!
	if(.csys.info == 'S') call sar_rad2km(info,getRerdeGreenFrequencyWavenumber(gfk))
!
!  green section
!
!	if (twodflag) then
!		errmsg = create2DGreenSection(greensec,gfk,info,dtp,dt,nordl,fcl,nordh,fch)
!	else
	errmsg = createGreenSection(greensec,gfk,info,dtp,dt,nordl,fcl,nordh,fch)
!	endif
	if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  convolve with source time function if desired
!
	if(ifqtf.eq.1) then
		lu = get(fuh)
		open(lu,file=stfasc)
		read(lu,*) nstf,dtstf
		allocate(stf(nstf))
		do i=1,nstf
			read(lu,*) dum,stf(i)
		enddo
		close(lu); call undo(fuh)
		if(abs(dtstf-dt)/dt.gt.1.e-6) then
			call new(errmsg,2,'dt of syns and stf inconsistent',myname)
			call print(errmsg); stop
		endif
		call convolveSourceTimeFunctionGreenSection(greensec,stf)
	endif
!
!  write to sff or ssa or sgy
!
	if (outformat == 'sff') then
		print *,'write to sff file'
		if (ifwin == 0) then
			call writeSFFGreenSection(greensec,get(fuh),seisfile)
		else
			call writeSFFGreenSection(greensec,get(fuh),seisfile,gv,tcut)
		endif
	else if (outformat == 'ssa') then
		print *,'write to ssa file'
		if (ifwin == 0) then
			call writeSSAGreenSection(greensec,get(fuh),seisfile)
		else
			call writeSSAGreenSection(greensec,get(fuh),seisfile,gv,tcut)
		endif
	else if (outformat == 'sgy') then
		print *,'write to sgy file'
		if (ifwin == 0) then
			call writeSGYGreenSection(greensec,get(fuh),seisfile)
		else
			call writeSGYGreenSection(greensec,get(fuh),seisfile,gv,tcut)
		endif
	else if (outformat == 'asc') then
		print *,'write to ascii file'
		if (ifwin == 0) then
			call writeAscGreenSection(greensec,get(fuh),seisfile)
		else
			call writeAscGreenSection(greensec,get(fuh),seisfile,gv,tcut)
		endif
	else
		call new(errmsg,2,'Unknown output format',myname)
		call print(errmsg); stop
	endif
!
!  deallocation
!
	if (allocated(stf)) deallocate(stf)
	call dealloc(gfk)
	call dealloc(info)
	call dealloc(errmsg)
	call dealloc(fuh)
!
 end program
!-----------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   GFDSVRKFSEIS'
	print *,''
	print *,'Usage: gfdsvrkfseis  [-h] [-o] [-f] [-hp] [-lp] [-e] [-gvw] [-dt] [-qtf] [-zs] [-p]'
	print *,'                     infofile dsvbasename'
	print *,''
	print *,'Compute synthetic seismograms for receivers listed in an <infofile>'
	print *,'using DSV frequency-wavenumber spectra read from files with'
	Print *,'names starting with <dsvbasename>.'
	print *,''
	print *,'<infofile> contains numeric station ids, names, components, network, locations,'
	print *,'source location, force or moment tensor components, and event time.'
	print *,'' 
	print *,'By default, the seismograms are NOT filtered.'
	print *,''
	print *,'Unit for displacement is nanometers/Nm'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-o    : name of output SFF file (def = greenseis.sff)'
	print *,'-f    : format of output: sff or ssa or sgy or asc, (default = sff)'
	print *,'-hp   : order and corner frequency (Hz) of highpass filter'
	print *,'        such as def = <0 10.> '
	print *,'-lp   : order and corner frequency (Hz) of lowpass filter'
	print *,'        such as def = <0 100.> '
	print *,'-e    : maximum distance up to which stations'
	print *,'        are considered (def = 100 km), unit km,'
	print *,'        but if sar_csys in infofile is S, unit must be degrees.'
	print *,'-gvw  : cut length of seismogram according to group velocity <gv> and minimum length <tcut>:'
	print *,'        tmax = tcut + dis/gv'
	print *,'-dt   : desired sampling interval (def = 0.001 s). Actual sampling may be'
	print *,'        different because of th 2**n requirement of the FFT. It is written'
	print *,'        into the header of the seismogram file. Never use the value given here'
	print *,'        when plotting seismograms!'
	print *,'-qtf  : name of source time file (either use -qsp or -qtf or nothing!)'
	print *,'-zs   : desired source depth in meter (def = from info file)'
	print *,'-vs   : verbosity flag for stream access input (def = false)'
	print *,'-p    : FK-spectra were calculated with a parallel run (def = .false.)'
!	print *,'-2d   : assume line source perepndicular to x-direction, use 2D formulae'
!	print *,'        station positions are projected onto x-direction'
	print '(50(1h-))'
	return
	end
