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
!--------------------------------------------------------
!  module for green seismogram section for either force or
!  moment tensor associated with seismogram traces in info
!  file
!--------------------------------------------------------
 module greenSection
	use greenFrequencyWavenumber
	use greenDisplacementSpectrum
	use frequencyTime
	use sourceReceiver
	use timeSeries
	use sffHeader
	use sffDatablock
	use ssaDataset
	use dataSegy
	use mathConstants
	use instrumentResponse
	use errorMessage
	implicit none
	interface dealloc; module procedure deallocGreenSection; end interface
	interface operator (.nts.); module procedure getNTimeSeriesGreenSection; end interface
	interface operator (.ts.); module procedure getSelectedTimeSeriesGreenSection; end interface
	type green_section
		private
		integer :: istyp                                              ! Quelltyp
		type (time_series), dimension(:), pointer :: seis => null()   ! Spuren
		type (sar), pointer :: info => null()                         ! Info zu sar
		type (green_frequency_wavenumber), pointer :: gfk => null()   ! fk-spectra
	end type
!
 contains
!-----------------------------------------------------------
!  create Green seismogram section
!  assume that gfk object is completely filled
!  dtp: desired sampling interval
!  dt:  really chosen sampling interval
!  nl,fcl,nh,fch: order and corner frequency (Hz)
!  of Butterworth low and high pass filters
!  response: (optional) array of instrument response objects ordered according to info file
!
	function createGreenSection(this,gfk,info,dtp,dt,nl,fcl,nh,fch,response,silent) result(errmsg)
	type (green_section) :: this
	type (green_frequency_wavenumber), target :: gfk
	type (sar), target :: info
	type (instrument_response), dimension(:), optional :: response
	type (error_message) :: errmsg
	real :: dtp,dt,fcl,fch
	integer :: nl,nh,nf
	logical, optional :: silent
	logical :: verbose
	real, dimension(:,:), allocatable :: besselj
	complex, dimension(:), allocatable :: zsp,hflow,hfhigh
	complex, dimension(:), pointer :: resp
	complex, dimension(:,:), allocatable :: zdis
	real, dimension(:), allocatable :: urs
	real :: force(3),mt(6)
	integer :: nf1,nf2,i,j,kintmax,nsamp,nwint
	integer, dimension(:), allocatable :: kint
	complex :: zoma
	character (len=132) :: myname = 'createGreenSection'
!
	if (.not. present(silent)) then; verbose = .true.; else; verbose = .not.silent; endif
!
	call new(errmsg,myname)
!
	this%istyp = .styp.gfk
	this%info => info
	this%gfk => gfk
!
	nf1 = .nfa.gfk; nf2 = .nfb.gfk
	nf = nf2-nf1+1
	zoma = getZomaGreenFrequencyWavenumber(gfk)
	kintmax = (1+.styp.gfk)*getMaxWavenumberIntegrals()
!
!  calculate kint needed for synthetics
!
	if (dsvSteppingGreenFrequencyWavenumber(gfk) == 1) then
		if (this%istyp == 1) then
			nwint = 18
			allocate(kint(18))
			kint = (/ (i,i=1,14),(i,i=43,46) /)    ! stresses available
		else
			nwint = 9
			allocate(kint(9))
			kint = (/ (i,i=1,7),(i,i=22,23) /)     ! stresses available
		endif
	else
		if (this%istyp == 1) then
			nwint = 14
			allocate(kint(14))
			kint = (/ (i,i=1,14) /)
		else
			nwint = 7
			allocate(kint(7))
			kint = (/ (i,i=1,7) /)
		endif
	endif
!
!  dt and nsamp
!
	call newNsampDtFrequencyTime(.df.gfk,dtp,nsamp,dt)
	if (verbose) print *,'<createGreenSection>: nsamp,dtp,dt: ',nsamp,dtp,dt
!
!  Filters
!
	allocate(hflow(nf),hfhigh(nf))
	call hbutcw("low ",nl,nf,2.*mc_pi*fcl,zoma,2.*mc_pi*(.df.gfk),hflow)    ! Filtering
	call hbutcw("high",nh,nf,2.*mc_pi*fch,zoma,2.*mc_pi*(.df.gfk),hfhigh)
!
!  Loop over receivers
!
	allocate(besselj(.nwnmax.gfk,0:2))
	allocate(zdis(nf2,kintmax))
	allocate(zsp(nf2))
	allocate(this%seis(.nr.info))
	allocate(urs(nsamp))
	allocate(resp(nf))
	do i=1,.nr.info
		call besselGreenFrequencyWavenumber(gfk,info.dis.i,besselj(:,0:2))
		do j=1,nwint
			if(this%istyp == 1) then
				call wnintMomentGreenFrequencyWavenumber(gfk,kint(j),info.dis.i,0.2,besselj(:,0:2),zdis(:,kint(j)))
			else
				call wnintForceGreenFrequencyWavenumber(gfk,kint(j),info.dis.i,0.2,besselj(:,0:2),zdis(:,kint(j)))
			endif
		enddo
		if(this%istyp == 0) then
			call sar_getForce(info,force)
			call forceGreenDisplacementSpectrum(nf1,nf2,.csys.info,info.comp.i,info.propdir.i, &
			                                  & info.phi.i,force,zdis,zsp)
		else if(this%istyp == 1) then
			call sar_getMomentTensor(info,mt)
			call momentGreenDisplacementSpectrum(nf1,nf2,.csys.info,info.comp.i,info.propdir.i, &
			                                   & info.phi.i,mt,zdis,zsp)
		else
			call new(errmsg,2,'unknown source type',myname)
			return
		endif
		zsp(nf1:nf2) = zsp(nf1:nf2)*hflow(1:nf)*hfhigh(1:nf)
		if (present(response)) then
			resp => getInstrumentResponse(response(i),nf1,nf2,.df.gfk,.sigma.gfk)
			zsp(nf1:nf2) = zsp(nf1:nf2)*resp(1:nf)*1.e-9     !  back to counts: d[c] = H* 1c/(m/s) *u* 1e-9*m/s
		endif
		errmsg = transformFrequencyTime(nf1,nf2,.df.gfk,.sigma.gfk,zsp,nsamp,dt,urs)
		if (.level.errmsg == 2) then
			call addTraceErrorMessage(errmsg,myname)
			return
		else if (.level.errmsg == 1) then
			call addTraceErrorMessage(errmsg,myname)
		endif
		call createFromDataTimeSeries(this%seis(i),nsamp,.therd.info,dt,urs)
		if (verbose) then
			if (present(response)) then
				print *,'Max values: ',trim(info.station.i),info.comp.i,absmaxTimeSeries(this%seis(i)),'counts'
			else
				print *,'Max values: ',trim(info.station.i),info.comp.i,absmaxTimeSeries(this%seis(i)),'nm or mPa'
			endif
		endif
	enddo
!
	deallocate(besselj,zsp,zdis,urs,hflow,hfhigh,resp)
	end function createGreenSection
!-----------------------------------------------------------
!  create 2D Green seismogram section
!  assume that gfk object is completely filled
!  dtp: desired sampling interval
!  dt:  really chosen sampling interval
!  nl,fcl,nh,fch: order and corner frequency (Hz)
!  of Butterworth low and high pass filters
!  response: (optional) array of instrument response objects ordered according to info file
!
!	function create2DGreenSection(this,gfk,info,dtp,dt,nl,fcl,nh,fch,response) result(errmsg)
!	type (green_section) :: this
!	type (green_frequency_wavenumber), target :: gfk
!	type (sar), target :: info
!	type (instrument_response), dimension(:), optional :: response
!	type (error_message) :: errmsg
!	real :: dtp,dt,fcl,fch
!	integer :: nl,nh,nf
!	complex, dimension(:), allocatable :: zsp,hflow,hfhigh
!	complex, dimension(:), pointer :: resp
!	complex, dimension(:,:,:), allocatable :: zdis
!	real, dimension(:), allocatable :: urs,x
!	real :: force(3),mt(6)
!	integer :: nf1,nf2,i,j,kintmax,nsamp,nwint,nr
!	integer, dimension(:), allocatable :: kint
!	complex :: zoma
!	character (len=132) :: myname = 'create2DGreenSection'
!
!	call new(errmsg,myname)
!
!	this%istyp = .styp.gfk
!	this%info => info
!	this%gfk => gfk
!
!	nf1 = .nfa.gfk; nf2 = .nfb.gfk
!	nf = nf2-nf1+1
!	zoma = getZomaGreenFrequencyWavenumber(gfk)
!	kintmax = (1+.styp.gfk)*getMaxWavenumberIntegrals()
!
!  Force or Moment
!
!	if (this%istyp == 0) then
!		call sar_getForce(info,force)
!	else if (this%istyp == 1) then
!		call sar_getMomentTensor(info,mt)
!	else
!		call new(errmsg,2,'unknown source type',myname)
!		return
!	endif
!
!  calculate kint needed for synthetics
!
!	if (dsvSteppingGreenFrequencyWavenumber(gfk) == 1) then
!		if (this%istyp == 1) then
!			nwint = 18
!			allocate(kint(18))
!			kint = (/ (i,i=1,14),(i,i=43,46) /)    ! stresses available
!		else
!			nwint = 9
!			allocate(kint(9))
!			kint = (/ (i,i=1,7),(i,i=22,23) /)     ! stresses available
!		endif
!	else
!		if (this%istyp == 1) then
!			nwint = 14
!			allocate(kint(14))
!			kint = (/ (i,i=1,14) /)
!		else
!			nwint = 7
!			allocate(kint(7))
!			kint = (/ (i,i=1,7) /)
!		endif
!	endif
!
!  dt and nsamp
!
!	call newNsampDtFrequencyTime(.df.gfk,dtp,nsamp,dt)
!	print *,'<create2DGreenSection>: nsamp,dtp,dt: ',nsamp,dtp,dt
!
!  Filters
!
!	allocate(hflow(nf),hfhigh(nf))
!	call hbutcw("low ",nl,nf,2.*mc_pi*fcl,zoma,2.*mc_pi*(.df.gfk),hflow)    ! Filtering
!	call hbutcw("high",nh,nf,2.*mc_pi*fch,zoma,2.*mc_pi*(.df.gfk),hfhigh)
!
!  Loop over receivers
!
!	nr = .nr.info
!	allocate(zdis(nf2,nr,kintmax))
!	allocate(zsp(nf2))
!	allocate(this%seis(nr),x(nr))
!	allocate(urs(nsamp))
!	allocate(resp(nf))
!	do i = 1,nr
!		x(i) = (info.dis.i)*cos(info.phi.i)
!	enddo
!	do j=1,nwint
!		if(this%istyp == 1) then
!			call wnintMomentGreenFrequencyWavenumber(gfk,kint(j),x,0.2,zdis(:,:,kint(j)))
!		else
!			call wnint2DForceGreenFrequencyWavenumber(gfk,kint(j),x,0.2,zdis(:,:,kint(j)))
!		endif
!	enddo
!	do i = 1,nr
!		if (this%istyp == 0) call  force2DGreenDisplacementSpectrum(nf1,nf2,info.comp.i,info.propdir.i,force,zdis(:,i,:),zsp)
!		if (this%istyp == 1) call moment2DGreenDisplacementSpectrum(nf1,nf2,info.comp.i,info.propdir.i,mt,zdis(:,i,:),zsp)
!		zsp(nf1:nf2) = zsp(nf1:nf2)*hflow(1:nf)*hfhigh(1:nf)
!		if (present(response)) then
!			resp => getInstrumentResponse(response(i),nf1,nf2,.df.gfk,.sigma.gfk)
!			zsp(nf1:nf2) = zsp(nf1:nf2)*resp(1:nf)*1.e-9     !  back to meters
!		endif
!		errmsg = transformFrequencyTime(nf1,nf2,.df.gfk,.sigma.gfk,zsp,nsamp,dt,urs)
!		if (.level.errmsg == 2) then
!			call addTraceErrorMessage(errmsg,myname)
!			return
!		else if (.level.errmsg == 1) then
!			call addTraceErrorMessage(errmsg,myname)
!		endif
!		call createFromDataTimeSeries(this%seis(i),nsamp,.therd.info,dt,urs)
!		print *,'Max values: ',trim(info.station.i),info.comp.i,absmaxTimeSeries(this%seis(i)),'nm or milli-Pa'
!	enddo
!
!	deallocate(zsp,zdis,urs,hflow,hfhigh,resp,x)
!	end function create2DGreenSection
!--------------------------------------------------------------------------
!  deallocate green section object
!
	subroutine deallocGreenSection(this)
	type (green_section) :: this
	integer :: i
	if (associated(this%seis)) then
		do i=1,size(this%seis); call dealloc(this%seis(i)); enddo
		deallocate(this%seis)
	endif
	end subroutine deallocGreenSection
!--------------------------------------------------------------------------
!  convolve Green section with a source time function
!
	subroutine convolveSourceTimeFunctionGreenSection(this,stf)
	type (green_section) :: this
	type (time_series) :: cts
	real, dimension(:) :: stf
	integer :: j
!
	do j = 1,size(this%seis)
		cts = this%seis(j)    ! just copy pointer to data
		this%seis(j) = convolveTimeSeries(cts,stf)
		call dealloc(cts)
	enddo
	end subroutine convolveSourceTimeFunctionGreenSection
!---------------------------------------------------------------------------
!  get a pointer to the time series
!
	function getTimeSeriesGreenSection(this) result(p)
	type (green_section), intent(in) :: this
	type (time_series), dimension(:), pointer :: p
	p => this%seis
	end function getTimeSeriesGreenSection
!---------------------------------------------------------------------------
!  get a pointer to the time series with index k
!
	function getSelectedTimeSeriesGreenSection(this,k) result(p)
	type (green_section), intent(in) :: this
	type (time_series), pointer :: p
	integer, intent(in) :: k
	p => this%seis(k)
	end function getSelectedTimeSeriesGreenSection
!---------------------------------------------------------------------------
!  set a selected time series of section to a new one
!
	subroutine setSelectedTimeSeriesGreenSection(this,ts,k)
	type (green_section), intent(in) :: this
	type (time_series) :: ts
	integer, intent(in) :: k
	call copyTimeSeries(ts,this%seis(k))
	end subroutine setSelectedTimeSeriesGreenSection
!---------------------------------------------------------------------------
!  get the number of time series in green section
!
	integer function getNTimeSeriesGreenSection(this)
	type (green_section), intent(in) :: this
	getNTimeSeriesGreenSection = size(this%seis)
	end function getNTimeSeriesGreenSection
!------------------------------------------------------------------------
!  write to sff file
!
!  optional arguments: 	gv, tbase for group velocity window
!
	subroutine writeSFFGreenSection(this,lu,filename,gv,tbase)
	type (green_section) :: this
	integer :: lu
	character (len=*) :: filename
	type (sff_header) :: head
	type (sff_datablock) :: dbl
	type (time_series), pointer :: ts
	logical :: last
	integer :: ntraces,j,nsamp
	real, optional :: gv,tbase
!
	call basicSFFHeader(head)
	call addSourceSFFHeader(head,'auto',.csys.(this%info),.xs.(this%info),.ys.(this%info), &
	                      & .zs.(this%info),.sdate.(this%info),.stime.(this%info))
	call writeSFFHeader(head,lu,filename)
	last = .false.
	ntraces = size(this%seis)
	do j=1,ntraces
		if (j == ntraces) last = .true.
		ts => this%seis(j)
		if (present(gv)) then
			nsamp = min(.nsamp.ts,ceiling( (tbase+((this%info).dis.j)/gv)/.dt.ts) )
		else
			nsamp = .nsamp.ts
		endif
		call basicSFFDatablock(dbl,nsamp,.dt.ts,.trace.ts,(this%info).netid.j, &
		                    & (this%info).station.j,'XU'//((this%info).comp.j),last)
		call modifyDateSFFDataBlock(dbl,.syear.(this%info),.smonth.(this%info),.sday.(this%info))
		call modifyTanfSFFDatablock(dbl,.therd.(this%info))
		call addInfoSFFDatablock(dbl,.csys.(this%info),(this%info).x.j,(this%info).y.j,0.0,1)
		call writeSFFDataBlock(dbl,lu,last)
		call nullifyDataSFFDatablock(dbl)
	enddo
	end subroutine writeSFFGreenSection
!---------------------------------------------------------------------------
!  write as ssa datset
!
	subroutine writeSSAGreenSection(this,lu,filename,gv,tbase)
	type (green_section) :: this
	integer :: lu
	character (len=*) :: filename
	type (ssa_dataset) :: dset
	type (ssa_source) :: source
	type (ssa_info) :: info
	type (ssa_datablock) :: dbl
	type (time_series), pointer :: ts
	integer :: ntraces,j,nsamp
	real, optional :: gv,tbase
!
	ntraces = size(this%seis)
	call createBasicSSADataset(dset,ntraces)
	call createSSASource(source,.csys.(this%info),.xs.(this%info),.ys.(this%info), &
	                     .zs.(this%info),.sdatefy.(this%info),.tfs.(this%info),.tns.(this%info))
	do j = 1,ntraces
		ts => this%seis(j)
		if (present(gv)) then
			nsamp = min(.nsamp.ts,ceiling( (tbase+((this%info).dis.j)/gv)/.dt.ts) )
		else
			nsamp = .nsamp.ts
		endif
		call createSSAInfo(info,.csys.(this%info),(this%info).x.j,(this%info).y.j,0.0,1)
		call createSSADatablock(dbl,nsamp,.dt.ts,.trace.ts,.sdatefy.(this%info),.tfs.(this%info),.tns.(this%info), &
		                       (this%info).station.j,'XU'//((this%info).comp.j),source,info)
		call addDatablockSSADataset(dset,dbl)
		call dealloc(dbl)
	enddo
	call writeSSADataset(dset,lu,filename,.false.)
	call dealloc(dset)
	end subroutine writeSSAGreenSection
!---------------------------------------------------------------------------
!  write as sgy datset
!
	subroutine writeSGYGreenSection(this,lu,filename,gv,tbase)
	type (green_section) :: this
	integer :: lu
	character (len=*) :: filename
	type (time_series), pointer :: ts
	type (data_segy) :: data
	integer :: ntraces,j
	integer(kind=2) :: csi,nsamp,hod,moh,som
	integer(kind=4) :: cscale,tscale
	real(kind=4) :: rcscale,rtscale
	character (len=1) :: csc
	real, optional :: gv,tbase
	character (len=132) :: myname = 'writeSGYGreenSection'
	character (len=132) :: tmpstring
!
	ntraces = size(this%seis)
	csc = .csys.(this%info)
	if (csc == 'C') then
		csi = 1
	else if (csc == 'S') then
		csi = 3
	else
		csi = 1
		write(*,*) trim(myname)//"unknown coordinate system. assuming cartesian system."
	end if
!
	hod = 0;moh = 0;som = 0
	if (real(.tfs.(this%info),4)/3600. .gt. 1.) then
		hod = floor(real(.tfs.(this%info),4)/3600.)
	end if
	if (real(.tfs.(this%info)-3600*hod,4)/60. .gt. 1.) then
		moh = floor(real(.tfs.(this%info)-3600*hod,4)/60.)
	end if
	if (real(.tfs.(this%info)-3600*hod-60*moh,4) .gt. 1.) then
		som = floor(real(.tfs.(this%info)-3600*hod-60*moh,4))
	end if
!
	cscale = -10000
	if (ceiling(log10(2.d0**exponent(.xs.(this%info)))) .gt. cscale) then
		cscale = ceiling(log10(2.d0**exponent(.xs.(this%info))))
	end if
	if (ceiling(log10(2.d0**exponent(.ys.(this%info)))) .gt. cscale) then
		cscale = ceiling(log10(2.d0**exponent(.xs.(this%info))))
	end if
	do j = 1,ntraces
		if (ceiling(log10(2.d0**exponent((this%info).x.j))) .gt. cscale) then
			cscale = ceiling(log10(2.d0**exponent((this%info).x.j)))
		end if
		if (ceiling(log10(2.d0**exponent((this%info).y.j))) .gt. cscale) then
			cscale = ceiling(log10(2.d0**exponent((this%info).y.j)))
		end if
	enddo
!
	if (cscale .gt. 9 .and. cscale .le. 13) then
		rcscale = 1*10**(cscale-9);cscale = 1*10**(cscale-9);rcscale = 1/rcscale
	else if (cscale .ge. 5 .and. cscale .le. 9) then
		rcscale = 1*10**(9-cscale);cscale = -1*10**(9-cscale)
	else if (cscale .lt. 5 .and. cscale .gt. -2) then
		rcscale = 1*10000;cscale = -1*10000
	else
		rcscale = 0;cscale = 0
		write(*,*) trim(myname)//": source/receiver coordinates out of range"
		write(*,*) trim(myname)//": file will be missing source/receiver coordinates"
	end if
!
	do j = 1,ntraces
		ts => this%seis(j)
		if (present(gv)) then
			nsamp = min(.nsamp.ts,ceiling( (tbase+((this%info).dis.j)/gv)/.dt.ts) )
		else
			nsamp = .nsamp.ts
		end if
		if (j .eq. 1) then
			rtscale = 1.
			if (nint(((.dt.ts)*1d6),4) .lt. 1 .or. nint(((.dt.ts)*1d6),4) .gt. 32767) then
				write(*,*) trim(myname)//": warning: sampling interval out of range [1,32767] µs"
				tscale = ceiling(log10(2.d0**exponent((.dt.ts)*1d6)))
				if (tscale .gt. 0) then
					rtscale = 10**(tscale-4);rtscale = 1/rtscale;tscale = 10**(tscale-4)
					write(tmpstring,'(i6)') nint((.dt.ts)*1d6*rtscale,2)
					write(*,*) trim(myname)//": warning: rescaled sampling interval to "//trim(tmpstring)//" µs"
					write(tmpstring,'(f8.0)') 1./rtscale
					write(*,*) trim(myname)//": warning: multiply with "//trim(tmpstring)//" to get real sampling interval"
				else
					rtscale = 10**(4-tscale);tscale = 10**(4-tscale)
					write(tmpstring,'(i6)') nint((.dt.ts)*1d6*rtscale,2)
					write(*,*) trim(myname)//": warning: rescaled sampling interval to "//trim(tmpstring)//" µs"
					write(tmpstring,'(f8.0)') rtscale
					write(*,*) trim(myname)//": warning: divide by "//trim(tmpstring)//" to get real sampling interval"
				end if
			end if
			call buildBinaryDataSegy(data,int(ntraces,2),fixedlen=int(0,2),skipext=int(0,2),&
		                               & gensampint=nint(((.dt.ts)*1d6*rtscale),2),&
                                               & gennumsamp=nsamp,meassys_opt=int(1,2))
		end if
		if (csi .eq. 3) then
		call buildTraceDataSegy(data,j,.trace.ts,trid_opt=int(1,2),&
				    & srcdepth_opt=abs(int(.zs.(this%info),4)),coordpow_opt = int(cscale,2),&
				    & xsrccoord_opt=nint(rcscale*90,4)-int(rcscale*(.ys.(this%info)),4),&
				    & ysrccoord_opt=nint(rcscale*(.xs.(this%info)),4),&
				    & xgrpcoord_opt=nint(rcscale*90,4)-int(rcscale*((this%info).y.j),4),&
                                    & ygrpcoord_opt=nint(rcscale*((this%info).x.j),4),&
				    & coordunits_opt=csi,numsamp_opt=nsamp,&
				    & hod_opt = hod, moh_opt = moh, som_opt = som,&
				    & sampint_opt=nint(((.dt.ts)*1d6*rtscale),2),valunits_opt=int(0,2),srctype_opt=int(0,2))
		else
		call buildTraceDataSegy(data,j,.trace.ts,trid_opt=int(1,2),&
				    & srcdepth_opt=abs(int(.zs.(this%info),4)),coordpow_opt = int(cscale,2),&
				    & xsrccoord_opt=nint(rcscale*(.xs.(this%info)),4),&
				    & ysrccoord_opt=nint(rcscale*(.ys.(this%info)),4),&
				    & xgrpcoord_opt=nint(rcscale*((this%info).x.j),4),&
                                    & ygrpcoord_opt=nint(rcscale*((this%info).y.j),4),&
				    & coordunits_opt=csi,numsamp_opt=nsamp,&
				    & hod_opt = hod, moh_opt = moh, som_opt = som,&
				    & sampint_opt=nint(((.dt.ts)*1d6*rtscale),2),valunits_opt=int(0,2),srctype_opt=int(0,2))
		end if
	enddo
	call writeDataSegy(data,lu,filename)
!
	end subroutine writeSGYGreenSection
!---------------------------------------------------------------------------
!  write as ascii datset
!
	subroutine writeAscGreenSection(this,lu,filename,gv,tbase)
	type (green_section) :: this
	integer :: lu
	character (len=*) :: filename
	integer :: ntraces,j,nsamp,i
	real, optional :: gv,tbase
	real, dimension(:), pointer :: vals
!
	ntraces = size(this%seis)
	open(lu, file = filename)
	do j = 1,ntraces
		vals => .trace.(this%seis(j))
		if (present(gv)) then
			nsamp = min(.nsamp.(this%seis(j)),ceiling( (tbase+((this%info).dis.j)/gv)/.dt.(this%seis(j))) )
		else
			nsamp = .nsamp.(this%seis(j))
		endif
		write(lu,*) (this%info).station.j,'XU'//((this%info).comp.j), &
		          & nsamp,.dt.(this%seis(j)),.therd.(this%info)
		do i = 1,nsamp
			write(lu,*) vals(i)
		enddo
	enddo
	close(lu)
	end subroutine writeAscGreenSection
!
 end module greenSection
