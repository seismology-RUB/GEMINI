!--------------------------------------------------------------------------
!	Copyright 2013 Nils Müller, Wolfgang Friederich
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
!
!> \mainpage Module package with segyformat utilities
!>
!> \section purpose
!> - reading segy-files
!> - writing segy-files
!> \section version
!> - last update: 01.06.2010
!> - supports SEG Y rev 1 (http://www.seg.org/SEGportalWEBproject/prod/SEG-Publications/Pub-Technical-Standards/Documents/seg_y_rev1.pdf)
!> - currently the following sampling formats are supported:
!>   i)		2 byte integer
!>   ii)	4 byte integer
!>   iii)	4 byte IEEE floating-point
!> - currently support for the following sampling formats are being tested:
!>   		4 byte IBM floating-point
!> - module expects files in big_endian mode
!-------------------------------------------------------------
!> \brief Module with segyformat utilities
!
module dataSegy
!
	use binarySegy
	use traceSegy
!
	implicit none
!
!-------------------------------------------------------------
!> \brief data_segy object
!
	type data_segy
		private
		type (binary_segy) :: binary_segy
		type (trace_segy), allocatable, dimension(:) :: trace_segy
	end type data_segy
!
!-------------------------------------------------------------
! interfaces & operators
!
	interface new; module procedure readDataSegy; end interface
	interface dealloc; module procedure deallocDataSegy; end interface
!
	interface operator (.numtr.); module procedure getNumTrDataSegy; end interface
	interface operator (.numtren.); module procedure getNumTrEnDataSegy; end interface
	interface operator (.gensampint.); module procedure getGenSampIntDataSegy; end interface
	interface operator (.gennumsamp.); module procedure getGenNumSampDataSegy; end interface
	interface operator (.sampformat.); module procedure getSampFormatDataSegy; end interface
	interface operator (.enfold.); module procedure getEnFoldDataSegy; end interface
	interface operator (.trsort.); module procedure getTrSortDataSegy; end interface
	interface operator (.meassys.); module procedure getMeasSysDataSegy; end interface
	interface operator (.sigpol.); module procedure getSigPolDataSegy; end interface
	interface operator (.revnr.); module procedure getRevNrDataSegy; end interface
	interface operator (.fixedlen.); module procedure getFixedLenDataSegy; end interface
	interface operator (.numext.); module procedure getNumExtDataSegy; end interface
!
	interface operator (.trnr.); module procedure getTrNrDataSegy; end interface
	interface operator (.trid.); module procedure getTrIdDataSegy; end interface
	interface operator (.recele.); module procedure getRecEleDataSegy; end interface
	interface operator (.surfele.); module procedure getSurfEleDataSegy; end interface
	interface operator (.srcdepth.); module procedure getSrcDepthDataSegy; end interface
	interface operator (.elepow.); module procedure getElePowDataSegy; end interface
	interface operator (.coordpow.); module procedure getCoordPowDataSegy; end interface
	interface operator (.timepow.); module procedure getTimePowDataSegy; end interface
	interface operator (.xsrccoord.); module procedure getXSrcCoordDataSegy; end interface
	interface operator (.ysrccoord.); module procedure getYSrcCoordDataSegy; end interface
	interface operator (.xgrpcoord.); module procedure getXGrpCoordDataSegy; end interface
	interface operator (.ygrpcoord.); module procedure getYGrpCoordDataSegy; end interface
	interface operator (.intcs.); module procedure getIntCoordUnitsDataSegy; end interface
	interface operator (.charcs.); module procedure getCharCoordUnitsDataSegy; end interface
	interface operator (.srccorr.); module procedure getSrcCorrDataSegy; end interface
	interface operator (.grpcorr.); module procedure getGrpCorrDataSegy; end interface
	interface operator (.recdelay.); module procedure getRecDelayDataSegy; end interface
	interface operator (.numsamp.); module procedure getNumSampDataSegy; end interface
	interface operator (.recyear.); module procedure getRecYearDataSegy; end interface
	interface operator (.sampint.); module procedure getSampIntDataSegy; end interface
	interface operator (.doy.); module procedure getDoyDataSegy; end interface
	interface operator (.hod.); module procedure getHodDataSegy; end interface
	interface operator (.moh.); module procedure getMohDataSegy; end interface
	interface operator (.som.); module procedure getSomDataSegy; end interface
	interface operator (.sam.); module procedure getSamDataSegy; end interface
	interface operator (.timecode.); module procedure getTimecodeDataSegy; end interface
	interface operator (.trweight.); module procedure getTrWeightDataSegy; end interface
	interface operator (.valunits.); module procedure getValUnitsDataSegy; end interface
	interface operator (.srctype.); module procedure getSrcTypeDataSegy; end interface
	interface operator (.tracelongint.); module procedure getTrace4IntDataSegy; end interface
	interface operator (.traceshortint.); module procedure getTrace2IntDataSegy; end interface
	interface operator (.traceieee.); module procedure getTraceIeeeDataSegy; end interface
	interface operator (.traceibm.); module procedure getTraceIbmDataSegy; end interface
!
	contains
!-------------------------------------------------------------
!> \brief Fill data_segy object with data from file <filename>
!> \param this		data_segy object
!> \param filename	filename of input-file (memory space for data will be allocated)
!
	subroutine readDataSegy(this,filename)
!
		type (data_segy), intent(inout) :: this
		character(len=*), intent(in) :: filename
		integer (kind = 1) :: status
!
		open(unit=11, file=filename, status='old', access='stream', action='read', iostat=status, convert='big_endian')
		if (status /= 0) then
			stop "mod dataSegy: Unable to open input-file"
		endif
!
		call readBinarySegy(this%binary_segy)
		call readTraceSegy(this%trace_segy,.numtr.this%binary_segy,.gennumsamp.this%binary_segy,.sampformat.this%binary_segy,&
				 & .bytespersamp.this%binary_segy,.fixedlen.this%binary_segy,.skipext.this%binary_segy)
!
		close (11, status='keep', iostat=status)
		if (status /= 0) then
			stop "mod dataSegy: Unable to close input-file"
		endif
!
	end subroutine readDataSegy
!-------------------------------------------------------------
!> \brief Fill binary_segy object in given data_segy object with given data
!> \param this			data_segy object
!> \param numtr			number of traces
!> \param fixedlen		indicates whether or not all traces have same sample interval and number of samples (1: yes, 0: no)
!> \param skipext		indicates whether or not file contains extended file headers (x: x extended file headers, 0: no)
!> \param gensampint		generic sampling interval
!> \param gennumsamp		generic number of samples per trace
!> \param revnr			segy revision number
!> \param numtren_opt		(optional) number of traces per ensemble
!> \param enfold_opt		(optional) ensemble fold
!> \param trsort_opt		(optional) trace sorting code
!> \param meassys_opt		(optional) measurement system
!> \param sigpol_opt		(optional) impulse signal polarity
!
	subroutine buildBinaryDataSegy(this,numtr,fixedlen,skipext,gensampint,gennumsamp,numtren_opt,enfold_opt,trsort_opt,&
				     & meassys_opt,sigpol_opt)
!
		type (data_segy), intent(inout) :: this
		integer (kind = 2), intent(in) :: numtr
		integer (kind = 2), intent(in) :: fixedlen
		integer (kind = 2), intent(in) :: skipext
		integer (kind = 2), intent(in) :: gensampint
		integer (kind = 2), intent(in) :: gennumsamp
		integer (kind = 2), intent(in), optional :: numtren_opt
		integer (kind = 2) :: numtren
		integer (kind = 2), intent(in), optional :: enfold_opt
		integer (kind = 2) :: enfold
		integer (kind = 2), intent(in), optional :: trsort_opt
		integer (kind = 2) :: trsort
		integer (kind = 2), intent(in), optional :: meassys_opt
		integer (kind = 2) :: meassys
		integer (kind = 2), intent(in), optional :: sigpol_opt
		integer (kind = 2) :: sigpol
		integer (kind = 2) :: revnr
!
		numtren = 1
		enfold = 1
		trsort = 1
		meassys = 1
		sigpol = 1
		revnr = 256
!
		if (present(numtren_opt)) numtren = numtren_opt
		if (present(enfold_opt)) enfold = enfold_opt
		if (present(trsort_opt)) trsort = trsort_opt
		if (present(meassys_opt)) meassys = meassys_opt
		if (present(sigpol_opt)) sigpol = sigpol_opt
!
		call buildBinarySegy(this%binary_segy,numtr,fixedlen,skipext,numtren,gensampint,gennumsamp,revnr,enfold,trsort,meassys,sigpol)
!
	end subroutine buildBinaryDataSegy
!-------------------------------------------------------------
!> \brief Convert trace of type integer (kind = 4) to real (kind = 4)
!> \param trace4int	input
!> \param traceieee	in-/output (pass array of same size as input)
!
	subroutine int4ToReal4DataSegy(trace4int,traceieee)
		integer (kind = 4), intent(in), dimension(:) :: trace4int
		real (kind = 4), intent(inout), dimension(:) :: traceieee

		integer (kind = 4) :: i

		do i=1, size(trace4int)
			traceieee(i) = float(trace4int(i))
		enddo		
	end subroutine int4ToReal4DataSegy
!-------------------------------------------------------------
!> \brief Convert trace of type integer (kind = 2) to real (kind = 4)
!> \param trace2int	input
!> \param traceieee	in-/output (pass array of same size as input)
!
	subroutine int2ToReal4DataSegy(trace2int,traceieee)
		integer (kind = 2), intent(in), dimension(:) :: trace2int
		real (kind = 4), intent(inout), dimension(:) :: traceieee

		integer (kind = 4) :: i

		do i=1, size(trace2int)
			traceieee(i) = float(trace2int(i))
		enddo			
	end subroutine int2ToReal4DataSegy
!-------------------------------------------------------------
!> \brief Fill trace_segy object in given data_segy object with given data
!> \param this		data_segy object
!> \param trnr		tracenumber
!> \param trace		array consisting of tracevalues
!> \param trid_opt	(optional) trace identification code
!> \param recele_opt	(optional) receiver group elevation
!> \param surfele_opt	(optional) surface elevation at source
!> \param srcdepth_opt	(optional) source depth below surface (a positive number)
!> \param elepow_opt	(optional) scalar to be applied to elevations to give real value (positiv->multiplier,negativ->divisor)
!> \param coordpow_opt	(optional) scalar to be applied to coordinates to give real value (positiv->multiplier,negativ->divisor)
!> \param timepow_opt	(optional) scalar to be applied to times to give real value (positiv->multiplier,negativ->divisor)
!> \param xsrccoord_opt	(optional) source coordinate: X
!> \param ysrccoord_opt	(optional) source coordinate: Y
!> \param xgrpcoord_opt	(optional) group coordinate: X
!> \param ygrpcoord_opt	(optional) group coordinate: Y
!> \param coordunits_opt(optional) coordinate units
!> \param srccorr_opt	(optional) source static correction in milliseconds
!> \param grpcorr_opt	(optional) group static correction in milliseconds
!> \param recdelay_opt	(optional) delay recording time
!> \param sampint	(optional) sampling interval (indispensable if fixedlen was set to 0 in binary header)
!> \param recyear_opt	(optional) year data recorded
!> \param doy_opt	(optional) day of year
!> \param hod_opt	(optional) hour of day
!> \param moh_opt	(optional) minute of hour
!> \param som_opt	(optional) second of minute
!> \param timecode_opt	(optional) time basis code
!> \param trweight_opt	(optional) trace weighting factor
!> \param valunits_opt	(optional) trace value measurement unit
!> \param srctype_opt	(optional) source type
!
	subroutine buildTraceDataSegy(this,trnr,trace,trid_opt,recele_opt,surfele_opt,srcdepth_opt,elepow_opt,coordpow_opt,&
				    & timepow_opt,xsrccoord_opt,ysrccoord_opt,xgrpcoord_opt,ygrpcoord_opt,coordunits_opt,&
				    & srccorr_opt,grpcorr_opt,recdelay_opt,numsamp_opt,sampint_opt,recyear_opt,doy_opt,&
				    & hod_opt,moh_opt,som_opt,timecode_opt,trweight_opt,valunits_opt,srctype_opt)
!
		type (data_segy), intent(inout) :: this
		integer (kind = 4), intent(in) :: trnr
		real (kind = 4), intent(in), dimension(:) :: trace
		integer (kind = 2), intent(in), optional :: trid_opt
		integer (kind = 2) :: trid
		integer (kind = 4), intent(in), optional :: recele_opt
		integer (kind = 4) :: recele
		integer (kind = 4), intent(in), optional :: surfele_opt
		integer (kind = 4) :: surfele
		integer (kind = 4), intent(in), optional :: srcdepth_opt
		integer (kind = 4) :: srcdepth
		integer (kind = 2), intent(in), optional :: elepow_opt
		integer (kind = 2) :: elepow
		integer (kind = 2), intent(in), optional :: coordpow_opt
		integer (kind = 2) :: coordpow
		integer (kind = 2), intent(in), optional :: timepow_opt
		integer (kind = 2) :: timepow
		integer (kind = 4), intent(in), optional :: xsrccoord_opt
		integer (kind = 4) :: xsrccoord
		integer (kind = 4), intent(in), optional :: ysrccoord_opt
		integer (kind = 4) :: ysrccoord
		integer (kind = 4), intent(in), optional :: xgrpcoord_opt
		integer (kind = 4) :: xgrpcoord
		integer (kind = 4), intent(in), optional :: ygrpcoord_opt
		integer (kind = 4) :: ygrpcoord
		integer (kind = 2), intent(in), optional :: coordunits_opt
		integer (kind = 2) :: coordunits
		integer (kind = 2), intent(in), optional :: srccorr_opt
		integer (kind = 2) :: srccorr
		integer (kind = 2), intent(in), optional :: grpcorr_opt
		integer (kind = 2) :: grpcorr
		integer (kind = 2), intent(in), optional :: recdelay_opt
		integer (kind = 2) :: recdelay
		integer (kind = 2), intent(in), optional :: recyear_opt
		integer (kind = 2) :: recyear
		integer (kind = 2), intent(in), optional :: numsamp_opt
		integer (kind = 2) :: numsamp
		integer (kind = 2), intent(in), optional :: sampint_opt
		integer (kind = 2) :: sampint
		integer (kind = 2), intent(in), optional :: doy_opt
		integer (kind = 2) :: doy
		integer (kind = 2), intent(in), optional :: hod_opt
		integer (kind = 2) :: hod
		integer (kind = 2), intent(in), optional :: moh_opt
		integer (kind = 2) :: moh
		integer (kind = 2), intent(in), optional :: som_opt
		integer (kind = 2) :: som
		integer (kind = 2), intent(in), optional :: timecode_opt
		integer (kind = 2) :: timecode
		integer (kind = 2), intent(in), optional :: trweight_opt
		integer (kind = 2) :: trweight
		integer (kind = 2), intent(in), optional :: valunits_opt
		integer (kind = 2) :: valunits
		integer (kind = 2), intent(in), optional :: srctype_opt
		integer (kind = 2) :: srctype
!
		trid = 0
		recele = 0
		surfele = 0
		srcdepth = 0
		elepow = 1
		coordpow = 1
		timepow = 1
		xsrccoord = 0
		ysrccoord = 0
		xgrpcoord = 0
		ygrpcoord = 0
		coordunits = 1
		srccorr = 0
		grpcorr= 0
		recdelay = 0
		numsamp = size(trace)
		sampint = .gensampint.this%binary_segy
		recyear = 0
		doy = 0
		hod = 0
		moh = 0
		som = 0
		timecode = 1
		trweight = 0
		valunits = 0
		srctype = 0
!
		if (present(trid_opt)) trid = trid_opt
		if (present(recele_opt)) recele = recele_opt
		if (present(surfele_opt)) surfele = surfele_opt
		if (present(srcdepth_opt)) srcdepth = srcdepth_opt
		if (present(elepow_opt)) elepow = elepow_opt
		if (present(coordpow_opt)) coordpow = coordpow_opt
		if (present(timepow_opt)) timepow = timepow_opt
		if (present(xsrccoord_opt)) xsrccoord = xsrccoord_opt
		if (present(ysrccoord_opt)) ysrccoord = ysrccoord_opt
		if (present(xgrpcoord_opt)) xgrpcoord = xgrpcoord_opt
		if (present(ygrpcoord_opt)) ygrpcoord = ygrpcoord_opt
		if (present(coordunits_opt)) coordunits = coordunits_opt
		if (present(srccorr_opt)) srccorr = srccorr_opt
		if (present(grpcorr_opt)) grpcorr = grpcorr_opt
		if (present(recdelay_opt)) recdelay = recdelay_opt
		if (present(numsamp_opt)) then
			numsamp = numsamp_opt
		elseif(.not. present(numsamp_opt) .and. .fixedlen.this == 0) then
			stop "mod dataSegy: if you set Fixed Length Trace Flag to 0, you have to state Number of Samples in each Traceheader"
		endif
		if (present(sampint_opt)) then
			sampint = sampint_opt
		elseif(.not. present(sampint_opt) .and. .fixedlen.this == 0) then
			stop "mod dataSegy: if you set Fixed Length Trace Flag to 0, you have to state Sampling Interval in each Traceheader"
		endif
		if (present(recyear_opt)) recyear = recyear_opt
		if (present(doy_opt)) doy = doy_opt
		if (present(hod_opt)) hod = hod_opt
		if (present(moh_opt)) moh = moh_opt
		if (present(som_opt)) som = som_opt
		if (present(timecode_opt)) timecode = timecode_opt
		if (present(trweight_opt)) trweight = trweight_opt
		if (present(valunits_opt)) valunits = valunits_opt
		if (present(srctype_opt)) srctype = srctype_opt
!
		call allocTraceSegy(this%trace_segy,.numtr.this%binary_segy)
!
		call buildTraceSegy(this%trace_segy(trnr),trnr,trace,trid,recele,surfele,srcdepth,elepow,coordpow,timepow,xsrccoord,&
				  & ysrccoord,xgrpcoord,ygrpcoord,coordunits,srccorr,grpcorr,recdelay,numsamp,sampint,recyear,doy,hod,moh,&
				  & som,timecode,trweight,valunits,srctype,.sampformat.this%binary_segy)
!
	end subroutine buildTraceDataSegy
!-------------------------------------------------------------
!> \brief Write content of data_segy object to file <filename>
!> \param this		data_segy object
!> \param filename	filename of output-file
!
	subroutine writeDataSegy(this,lu,filename)
!
		type (data_segy), intent(in) :: this
		integer, intent(in) :: lu
		character(len=*), intent(in) :: filename
		integer (kind = 1) :: status
!
		open(unit=lu, file=filename, status="replace", access="stream", action="write", convert="big_endian", iostat=status)
		if (status /= 0) then
			stop "mod dataSegy: Unable to open input-file"
		endif
!
		call writeBinarySegy(this%binary_segy,lu)
!
		call writeTraceSegy(this%trace_segy,lu,.numtr.this%binary_segy,.gennumsamp.this%binary_segy,&
			& .sampformat.this%binary_segy,.bytespersamp.this%binary_segy,.skipext.this%binary_segy,&
			& .fixedlen.this%binary_segy)
!
		close (lu, status="keep", iostat=status)
		if (status /= 0) then
			stop "mod dataSegy: Unable to close input-file"
		endif
!
	end subroutine writeDataSegy
!-------------------------------------------------------------
!> \brief Dealloc memory associated with data_segy object
!> \param this		data_segy object
!
	subroutine deallocDataSegy(this)
!
		type (data_segy) :: this
!
		call deallocTraceSegy(this%trace_segy,.numtr.this%binary_segy,.sampformat.this%binary_segy)
!
	end subroutine deallocDataSegy
!-------------------------------------------------------------
!> \brief Get number of traces
!> \param this		data_segy object
!
	integer (kind = 2) function getNumTrDataSegy(this)
		type (data_segy), intent(in) :: this
		getNumTrDataSegy = getNumTrBinarySegy(this%binary_segy)
	end function getNumTrDataSegy
!-------------------------------------------------------------
!> \brief Get number of traces per ensemble
!> \param this		data_segy object
!
	integer (kind = 2) function getNumTrEnDataSegy(this)
		type (data_segy), intent(in) :: this
		getNumTrEnDataSegy = getNumTrenBinarySegy(this%binary_segy)
	end function getNumTrEnDataSegy
!-------------------------------------------------------------
!> \brief Get generic sampling interval (if fixedlen in binary header is set to 1)
!> \param this		data_segy object
!
	integer (kind = 2) function getGenSampIntDataSegy(this)
		type (data_segy), intent(in) :: this
		getGenSampIntDataSegy = getGenSampIntBinarySegy(this%binary_segy)
	end function getGenSampIntDataSegy
!-------------------------------------------------------------
!> \brief Get number of samples per trace (if fixedlen in binary header is set to 1)
!> \param this		data_segy object
!
	integer (kind = 2) function getGenNumSampDataSegy(this)
		type (data_segy), intent(in) :: this
		getGenNumSampDataSegy = getGenNumSampBinarySegy(this%binary_segy)
	end function getGenNumSampDataSegy
!-------------------------------------------------------------
!> \brief Get sampling format code
!> \param this		data_segy object
!
	integer (kind = 2) function getSampFormatDataSegy(this)
		type (data_segy), intent(in) :: this
		getSampFormatDataSegy = getSampFormatBinarySegy(this%binary_segy)
	end function getSampFormatDataSegy
!-------------------------------------------------------------
!> \brief Get ensemble fold
!> \param this		data_segy object
!
	integer (kind = 2) function getEnFoldDataSegy(this)
		type (data_segy), intent(in) :: this
		getEnFoldDataSegy = getEnFoldBinarySegy(this%binary_segy)
	end function getEnFoldDataSegy
!-------------------------------------------------------------
!> \brief Get trace sorting code
!> \param this		data_segy object
!
	integer (kind = 2) function getTrSortDataSegy(this)
		type (data_segy), intent(in) :: this
		getTrSortDataSegy = getTrSortBinarySegy(this%binary_segy)
	end function getTrSortDataSegy
!-------------------------------------------------------------
!> \brief Get measurement system
!> \param this		data_segy object
!
	integer (kind = 2) function getMeasSysDataSegy(this)
		type (data_segy), intent(in) :: this
		getMeasSysDataSegy = getMeasSysBinarySegy(this%binary_segy)
	end function getMeasSysDataSegy
!-------------------------------------------------------------
!> \brief Get impulse signal polarity
!> \param this		data_segy object
!
	integer (kind = 2) function getSigPolDataSegy(this)
		type (data_segy), intent(in) :: this
		getSigPolDataSegy = getSigPolBinarySegy(this%binary_segy)
	end function getSigPolDataSegy
!-------------------------------------------------------------
!> \brief Get segy format revision number
!> \param this		data_segy object
!
	real (kind = 4) function getRevNrDataSegy(this)
		type (data_segy), intent(in) :: this
		integer(kind=2) :: q
		q = getRevNrBinarySegy(this%binary_segy)
		getRevNrDataSegy = q/256
	end function getRevNrDataSegy
!-------------------------------------------------------------
!> \brief Get fixed lenght indicator
!> \param this		data_segy object
!
	integer (kind = 2) function getFixedLenDataSegy(this)
		type (data_segy), intent(in) :: this
		getFixedLenDataSegy = getFixedLenBinarySegy(this%binary_segy)
	end function getFixedLenDataSegy
!-------------------------------------------------------------
!> \brief Get number of extended file headers
!> \param this		data_segy object
!
	integer (kind = 2) function getNumExtDataSegy(this)
		type (data_segy), intent(in) :: this
		getNumExtDataSegy = getNumExtBinarySegy(this%binary_segy)
	end function getNumExtDataSegy
!-------------------------------------------------------------
!> \brief Get trace number
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getTrNrDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getTrNrDataSegy = getTrNrTraceSegy(this%trace_segy(n))
	end function getTrNrDataSegy
!-------------------------------------------------------------
!> \brief Get trace identification code
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getTrIdDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getTrIdDataSegy = getTrIdTraceSegy(this%trace_segy(n))
	end function getTrIdDataSegy
!-------------------------------------------------------------
!> \brief Get receiver group elevation
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getRecEleDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getRecEleDataSegy = getRecEleTraceSegy(this%trace_segy(n))
	end function getRecEleDataSegy
!-------------------------------------------------------------
!> \brief Get surface elevation at source
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSurfEleDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSurfEleDataSegy = getSurfEleTraceSegy(this%trace_segy(n))
	end function getSurfEleDataSegy
!-------------------------------------------------------------
!> \brief Get source depth below surface (a positive number)
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSrcDepthDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSrcDepthDataSegy = getSrcDepthTraceSegy(this%trace_segy(n))
	end function getSrcDepthDataSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to elevations to give real value (positiv->multiplier,negativ->divisor)
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getElePowDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getElePowDataSegy = getElePowTraceSegy(this%trace_segy(n))
	end function getElePowDataSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to coordinates to give real value (positiv->multiplier,negativ->divisor)
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getCoordPowDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getCoordPowDataSegy = getCoordPowTraceSegy(this%trace_segy(n))
	end function getCoordPowDataSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to times to give real value (positiv->multiplier,negativ->divisor)
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getTimePowDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getTimePowDataSegy = getTimePowTraceSegy(this%trace_segy(n))
	end function getTimePowDataSegy
!-------------------------------------------------------------
!> \brief Get source coordinate: X
!> \param this		data_segy object
!> \param n		trace number
!
	real (kind = 4) function getXSrcCoordDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getXSrcCoordDataSegy = getXSrcCoordTraceSegy(this%trace_segy(n))
	end function getXSrcCoordDataSegy
!-------------------------------------------------------------
!> \brief Get source coordinate: Y
!> \param this		data_segy object
!> \param n		trace number
!
	real (kind = 4) function getYSrcCoordDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getYSrcCoordDataSegy = getYSrcCoordTraceSegy(this%trace_segy(n))
	end function getYSrcCoordDataSegy
!-------------------------------------------------------------
!> \brief Get group coordinate: X
!> \param this		data_segy object
!> \param n		trace number
!
	real (kind = 4) function getXGrpCoordDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getXGrpCoordDataSegy = getXGrpCoordTraceSegy(this%trace_segy(n))
	end function getXGrpCoordDataSegy
!-------------------------------------------------------------
!> \brief Get group coordinate: Y
!> \param this		data_segy object
!> \param n		trace number
!
	real (kind = 4) function getYGrpCoordDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getYGrpCoordDataSegy = getYGrpCoordTraceSegy(this%trace_segy(n))
	end function getYGrpCoordDataSegy
!-------------------------------------------------------------
!> \brief Get coordinate units (stored in integer)
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getIntCoordUnitsDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getIntCoordUnitsDataSegy = getCoordUnitsTraceSegy(this%trace_segy(n))
	end function getIntCoordUnitsDataSegy
!-------------------------------------------------------------
!> \brief Get coordinate units (stored in character)
!> \param this		data_segy object
!> \param n		trace number
!
	character (len=1) function getCharCoordUnitsDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		if (getCoordUnitsTraceSegy(this%trace_segy(n)) == 1) then
			getCharCoordUnitsDataSegy = 'C'
		elseif (getCoordUnitsTraceSegy(this%trace_segy(n)) == 2) then
			getCharCoordUnitsDataSegy = 'A'
		elseif (getCoordUnitsTraceSegy(this%trace_segy(n)) == 3) then
			getCharCoordUnitsDataSegy = 'S'
		elseif (getCoordUnitsTraceSegy(this%trace_segy(n)) == 4) then
			getCharCoordUnitsDataSegy = 'D'
		else
			getCharCoordUnitsDataSegy = 'X'
		endif
	end function getCharCoordUnitsDataSegy
!-------------------------------------------------------------
!> \brief Get source static correction in milliseconds
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSrcCorrDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSrcCorrDataSegy = getSrcCorrTraceSegy(this%trace_segy(n))
	end function getSrcCorrDataSegy
!-------------------------------------------------------------
!> \brief Get group static correction in milliseconds
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getGrpCorrDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getGrpCorrDataSegy = getGrpCorrTraceSegy(this%trace_segy(n))
	end function getGrpCorrDataSegy
!-------------------------------------------------------------
!> \brief Get delay recording time
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getRecDelayDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getRecDelayDataSegy = getRecDelayTraceSegy(this%trace_segy(n))
	end function getRecDelayDataSegy
!-------------------------------------------------------------
!> \brief Get sampling interval
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSampIntDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		if (getFixedLenBinarySegy(this%binary_segy) == 1) then
			getSampIntDataSegy = getGenSampIntBinarySegy(this%binary_segy)
		else
			getSampIntDataSegy = getSampIntTraceSegy(this%trace_segy(n))
		endif
	end function getSampIntDataSegy
!-------------------------------------------------------------
!> \brief Get number of samples per trace
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getNumSampDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		if (getFixedLenBinarySegy(this%binary_segy) == 1) then
			getNumSampDataSegy = getGenNumSampBinarySegy(this%binary_segy)
		else
			getNumSampDataSegy = getNumSampTraceSegy(this%trace_segy(n))
		endif
	end function getNumSampDataSegy
!-------------------------------------------------------------
!> \brief Get year data recorded
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getRecYearDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getRecYearDataSegy = getRecYearTraceSegy(this%trace_segy(n))
	end function getRecYearDataSegy
!-------------------------------------------------------------
!> \brief Get day of year
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getDoyDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getDoyDataSegy = getDoyTraceSegy(this%trace_segy(n))
	end function getDoyDataSegy
!-------------------------------------------------------------
!> \brief Get hour of day
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getHodDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getHodDataSegy = getHodTraceSegy(this%trace_segy(n))
	end function getHodDataSegy
!-------------------------------------------------------------
!> \brief Get minute of hour
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getMohDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getMohDataSegy = getMohTraceSegy(this%trace_segy(n))
	end function getMohDataSegy
!-------------------------------------------------------------
!> \brief Get seconds of minute
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSomDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSomDataSegy = getSomTraceSegy(this%trace_segy(n))
	end function getSomDataSegy
!-------------------------------------------------------------
!> \brief Get seconds after midnight
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSamDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSamDataSegy = (this.som.n) + 60*((this.moh.n) + 60*(this.hod.n))
	end function getSamDataSegy
!-------------------------------------------------------------
!> \brief Get time basis code
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getTimecodeDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getTimecodeDataSegy = getTimecodeTraceSegy(this%trace_segy(n))
	end function getTimecodeDataSegy
!-------------------------------------------------------------
!> \brief Get trace weighting factor
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getTrWeightDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getTrWeightDataSegy = getTrWeightTraceSegy(this%trace_segy(n))
	end function getTrWeightDataSegy
!-------------------------------------------------------------
!> \brief Get trace value measurement unit
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getValUnitsDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getValUnitsDataSegy = getValUnitsTraceSegy(this%trace_segy(n))
	end function getValUnitsDataSegy
!-------------------------------------------------------------
!> \brief Get source type
!> \param this		data_segy object
!> \param n		trace number
!
	integer (kind = 4) function getSrcTypeDataSegy(this,n)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		getSrcTypeDataSegy = getSrcTypeTraceSegy(this%trace_segy(n))
	end function getSrcTypeDataSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 2; you will receive a pointer)
!> \param this		data_segy object
!> \param n		trace number
!
	function getTrace4IntDataSegy(this,n) result(p)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		integer (kind = 4), dimension(:), pointer :: p
		p => getTrace4IntTraceSegy(this%trace_segy,n)
	end function getTrace4IntDataSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 3; you will receive a pointer)
!> \param this		data_segy object
!> \param n		trace number
!
	function getTrace2IntDataSegy(this,n) result(p)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		integer (kind = 2), dimension(:), pointer :: p
		p => getTrace2IntTraceSegy(this%trace_segy,n)
	end function getTrace2IntDataSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 5; you will receive a pointer)
!> \param this		data_segy object
!> \param n		trace number
!
	function getTraceIeeeDataSegy(this,n) result(p)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		real (kind = 4), dimension(:), pointer :: p
		p => getTraceIeeeTraceSegy(this%trace_segy,n)
	end function getTraceIeeeDataSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 1; you will receive a pointer)
!> \param this		data_segy object
!> \param n		trace number
!
	function getTraceIbmDataSegy(this,n) result(p)
		type(data_segy), intent(in) :: this
		integer (kind = 4), intent(in) :: n
		real (kind = 4), dimension(:), pointer :: p
		p => getTraceIbmTraceSegy(this%trace_segy,n)
	end function getTraceIbmDataSegy
!-------------------------------------------------------------
end module dataSegy
