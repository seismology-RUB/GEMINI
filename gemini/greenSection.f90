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
	function createGreenSection(this,gfk,info,dtp,dt,nl,fcl,nh,fch,response) result(errmsg)
	type (green_section) :: this
	type (green_frequency_wavenumber), target :: gfk
	type (sar), target :: info
	type (instrument_response), dimension(:), optional :: response
	type (error_message) :: errmsg
	real :: dtp,dt,fcl,fch
	integer :: nl,nh,nf
	real, dimension(:,:), allocatable :: besselj
	complex, dimension(:), allocatable :: zsp,hflow,hfhigh,resp
	complex, dimension(:,:), allocatable :: zdis
	real, dimension(:), allocatable :: urs
	real :: force(3),mt(6)
	integer :: nf1,nf2,i,j,kintmax,nsamp,nwint
	integer, dimension(:), allocatable :: kint
	complex :: zoma
	character (len=132) :: myname = 'createGreenSection'
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
	print *,'<createGreenSection>: nsamp,dtp,dt: ',nsamp,dtp,dt
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
			call getInstrumentResponse(response(i),.syear.info,.sdoy.info,nf1,nf2,.df.gfk,resp)
			zsp(nf1:nf2) = zsp(nf1:nf2)*resp(1:nf)*1.e-9     !  back to meters
		endif
		errmsg = transformFrequencyTime(nf1,nf2,.df.gfk,.sigma.gfk,zsp,nsamp,dt,urs)
		if (.level.errmsg == 2) then
			call addTraceErrorMessage(errmsg,myname)
			return
		else if (.level.errmsg == 1) then
			call addTraceErrorMessage(errmsg,myname)
		endif
		call createFromDataTimeSeries(this%seis(i),nsamp,.therd.info,dt,urs)
		print *,'Max values: ',trim(info.station.i),info.comp.i,absmaxTimeSeries(this%seis(i)),'nm or milli-Pa'
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
	function create2DGreenSection(this,gfk,info,dtp,dt,nl,fcl,nh,fch,response) result(errmsg)
	type (green_section) :: this
	type (green_frequency_wavenumber), target :: gfk
	type (sar), target :: info
	type (instrument_response), dimension(:), optional :: response
	type (error_message) :: errmsg
	real :: dtp,dt,fcl,fch
	integer :: nl,nh,nf
	complex, dimension(:), allocatable :: zsp,hflow,hfhigh,resp
	complex, dimension(:,:,:), allocatable :: zdis
	real, dimension(:), allocatable :: urs,x
	real :: force(3),mt(6)
	integer :: nf1,nf2,i,j,kintmax,nsamp,nwint,nr
	integer, dimension(:), allocatable :: kint
	complex :: zoma
	character (len=132) :: myname = 'create2DGreenSection'
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
!  Force or Moment
!
	if (this%istyp == 0) then
		call sar_getForce(info,force)
	else if (this%istyp == 1) then
		call sar_getMomentTensor(info,mt)
	else
		call new(errmsg,2,'unknown source type',myname)
		return
	endif
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
	print *,'<create2DGreenSection>: nsamp,dtp,dt: ',nsamp,dtp,dt
!
!  Filters
!
	allocate(hflow(nf),hfhigh(nf))
	call hbutcw("low ",nl,nf,2.*mc_pi*fcl,zoma,2.*mc_pi*(.df.gfk),hflow)    ! Filtering
	call hbutcw("high",nh,nf,2.*mc_pi*fch,zoma,2.*mc_pi*(.df.gfk),hfhigh)
!
!  Loop over receivers
!
	nr = .nr.info
	allocate(zdis(nf2,nr,kintmax))
	allocate(zsp(nf2))
	allocate(this%seis(nr),x(nr))
	allocate(urs(nsamp))
	allocate(resp(nf))
	do i = 1,nr
		x(i) = (info.dis.i)*cos(info.phi.i)
	enddo
	do j=1,nwint
		if(this%istyp == 1) then
!			call wnintMomentGreenFrequencyWavenumber(gfk,kint(j),x,0.2,zdis(:,:,kint(j)))
		else
			call wnint2DForceGreenFrequencyWavenumber(gfk,kint(j),x,0.2,zdis(:,:,kint(j)))
		endif
	enddo
	do i = 1,nr
		if (this%istyp == 0) call  force2DGreenDisplacementSpectrum(nf1,nf2,info.comp.i,info.propdir.i,force,zdis(:,i,:),zsp)
!		if (this%istyp == 1) call moment2DGreenDisplacementSpectrum(nf1,nf2,info.comp.i,info.propdir.i,mt,zdis(:,i,:),zsp)
		zsp(nf1:nf2) = zsp(nf1:nf2)*hflow(1:nf)*hfhigh(1:nf)
		if (present(response)) then
			call getInstrumentResponse(response(i),.syear.info,.sdoy.info,nf1,nf2,.df.gfk,resp)
			zsp(nf1:nf2) = zsp(nf1:nf2)*resp(1:nf)*1.e-9     !  back to meters
		endif
		errmsg = transformFrequencyTime(nf1,nf2,.df.gfk,.sigma.gfk,zsp,nsamp,dt,urs)
		if (.level.errmsg == 2) then
			call addTraceErrorMessage(errmsg,myname)
			return
		else if (.level.errmsg == 1) then
			call addTraceErrorMessage(errmsg,myname)
		endif
		call createFromDataTimeSeries(this%seis(i),nsamp,.therd.info,dt,urs)
		print *,'Max values: ',trim(info.station.i),info.comp.i,absmaxTimeSeries(this%seis(i)),'nm or milli-Pa'
	enddo
!
	deallocate(zsp,zdis,urs,hflow,hfhigh,resp,x)
	end function create2DGreenSection
!--------------------------------------------------------------------------
!  deallocate basisGreenSeis
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
!
 end module greenSection
