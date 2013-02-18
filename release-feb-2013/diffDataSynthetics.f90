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
 program diffDataSynthetics
	use sffDataSection
	use commandLine
	use timeSeries
	implicit none
!
	type (cmdLine) :: cl
	type (sff_datasection) :: dsec, ssec
	type (sff_datablock), pointer :: dbl,sbl
	type (time_series) :: dat_ts,syn_ts
	character (len=132) :: datafile,synfile,misfitfile,string
	real, dimension(:), allocatable :: d
	real, dimension(:), pointer :: sigma
	logical :: scaleflag, envflag
	integer :: nsamp_total,nsamp,ntrace,ntraced,ntraces
	integer :: j,jd
	real :: misfit,mag,dmax,smax
	external printhelp
!----------------------------------------------------------------------------
!  default values
!
	data mag/1.0/scaleflag/.false./envflag/.false./misfitfile/'misfit'/
!-------------------------------------------------------------------------------------------
	call new(cl,5,2,'h mag sc env o','0 1 0 0 1',printhelp)
	if(clOptset(cl,2)) then; string=clOptarg(cl,2); read(string,*) mag; endif
	if(clOptset(cl,3)) then; scaleflag = .true.; endif
	if(clOptset(cl,4)) then; envflag = .true.; endif
	if(clOptset(cl,5)) then; misfitfile = clOptarg(cl,5); endif
	datafile = clManarg(cl,1)
	synfile = clManarg(cl,2)
	call dealloc(cl)
!
!  read sections with data, synthetics
!
	call new(dsec,1,datafile)
	call new(ssec,2,synfile)
!
!  number of traces
!
	ntraced = .ndbl.dsec; ntraces = .ndbl.ssec
	if(ntraced /= ntraces) then
		print *,'Number of traces in data and synthetics are different'
		stop
	endif
	ntrace = ntraced
	nsamp_total = ntotSamplesSFFDataSection(dsec)
!
!  allocate space for sigma weighted residual data vector d, sigma and bm
!
	allocate(d(nsamp_total))
!
!  fill data residual vector
!
	jd = 0
	do j = 1,ntrace
		dbl => dsec.datablock.j
		nsamp = .nsamp.dbl
		call createLinkTimeSeries(dat_ts,nsamp,.tanf.dbl,.dt.dbl,trace(dbl))
		sbl => ssec.datablock.j
		call createLinkTimeSeries(syn_ts,.nsamp.sbl,.tanf.sbl,.dt.sbl,trace(sbl))
		if(nsamp /= .nsamp.sbl) then
			print *,'Sample number of data and syns differ! ',nsamp,.nsamp.sbl
			stop
		endif
!
!  calculate sigma weighted residual seimograms for this data block
!
		sigma => sigmaEnvelopeTimeSeries(dat_ts,mag,5)
!
!  overwrite dat_ts and syn_ts with envelope if desired
!
		if(envflag) then
			dat_ts = envelopeTimeSeries(dat_ts)
			syn_ts = envelopeTimeSeries(syn_ts)
		endif
		dmax = absmaxTimeSeries(dat_ts)
		smax = absmaxTimeSeries(syn_ts)
		if(.not. scaleflag) then
			d(jd+1:jd+nsamp) = (.trace.dat_ts - .trace.syn_ts)/sigma
		else
			d(jd+1:jd+nsamp) = (.trace.dat_ts-dmax*(.trace.syn_ts)/smax)/sigma
		endif
		call dealloc(dat_ts); call dealloc(syn_ts)
!
!  update sample base index
!
		jd=jd+nsamp
		deallocate(sigma)
	enddo
!
!  check size of d and nsamp_total
!
	if(size(d) /= nsamp_total) then
		print *,size(d),nsamp_total
		print *,'Something wrong with size of data vector'
	endif
!
!  calculate misfit
!
	misfit = dot_product(d,d)/nsamp_total
!
!  write to file
!
	open(3,file=misfitfile)
	write(3,*) misfit
	close(3)
!
	print *,'<diffDataSynthetics>: Misift written to file'
!
	deallocate(d)
	call dealloc(dsec); call dealloc(ssec)
!
 end program
!----------------------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                diffDataSynthetics'
	print *,''
	print *,'Usage: diffDataSynthetics [-h] [-mag] [-sc] [-env] [-o] datafile synfile'
	print *,''
	print *,'Evaluate misfit between data and synthetics in different ways'
	print *,'according to options'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-mag  : max magnification of signal achieved with envelope'
	print *,'        weighting (def=1)'
	print *,'-sc   : scale synthetics to data (def = 0)'
	print *,'-env  : use envelopes instead of time series (def=no)'
	print *,'-o    : name of output misfit file (def = misfit)'
	print '(50(1h-))'
	end subroutine printhelp

