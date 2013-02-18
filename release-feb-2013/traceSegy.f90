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
!-------------------------------------------------------------
!> \brief Sub-Module with segyformat utilities
!
module traceSegy
!
	implicit none
!
	interface operator (.trnr.); module procedure getTrNrTraceSegy; end interface
	interface operator (.trid.); module procedure getTrIdTraceSegy; end interface
	interface operator (.recele.); module procedure getRecEleTraceSegy; end interface
	interface operator (.surfele.); module procedure getSurfEleTraceSegy; end interface
	interface operator (.srcdepth.); module procedure getSrcDepthTraceSegy; end interface
	interface operator (.elepow.); module procedure getElePowTraceSegy; end interface
	interface operator (.coordpow.); module procedure getCoordPowTraceSegy; end interface
	interface operator (.timepow.); module procedure getTimePowTraceSegy; end interface
	interface operator (.xsrccoord.); module procedure getXSrcCoordTraceSegy; end interface
	interface operator (.ysrccoord.); module procedure getYSrcCoordTraceSegy; end interface
	interface operator (.xgrpcoord.); module procedure getXGrpCoordTraceSegy; end interface
	interface operator (.ygrpcoord.); module procedure getYGrpCoordTraceSegy; end interface
	interface operator (.coordunits.); module procedure getCoordUnitsTraceSegy; end interface
	interface operator (.srccorr.); module procedure getSrcCorrTraceSegy; end interface
	interface operator (.grpcorr.); module procedure getGrpCorrTraceSegy; end interface
	interface operator (.recdelay.); module procedure getRecDelayTraceSegy; end interface
	interface operator (.numsamp.); module procedure getNumSampTraceSegy; end interface
	interface operator (.recyear.); module procedure getRecYearTraceSegy; end interface
	interface operator (.doy.); module procedure getDoyTraceSegy; end interface
	interface operator (.hod.); module procedure getHodTraceSegy; end interface
	interface operator (.moh.); module procedure getMohTraceSegy; end interface
	interface operator (.som.); module procedure getSomTraceSegy; end interface
	interface operator (.timecode.); module procedure getTimecodeTraceSegy; end interface
	interface operator (.trweight.); module procedure getTrWeightTraceSegy; end interface
	interface operator (.valunits.); module procedure getValUnitsTraceSegy; end interface
	interface operator (.srctype.); module procedure getSrcTypeTraceSegy; end interface
	interface operator (.tracelongint.); module procedure getTrace4IntTraceSegy; end interface
	interface operator (.traceshortint.); module procedure getTrace2IntTraceSegy; end interface
	interface operator (.traceieee.); module procedure getTraceIeeeTraceSegy; end interface
	interface operator (.traceibm.); module procedure getTraceIbmTraceSegy; end interface
!
!-------------------------------------------------------------
!> \brief traceheader object
!
	type traceheader
		private
!> number of trace
		integer (kind = 4) :: trnr
!> trace identification code
		integer (kind = 2) :: trid
!> receiver group elevation
		integer (kind = 4) :: recele
!> surface elevation at source
		integer (kind = 4) :: surfele
!> source depth below surface
		integer (kind = 4) :: srcdepth
		integer (kind = 2) :: elepow
		integer (kind = 2) :: coordpow
		integer (kind = 2) :: timepow
!> source coordinate x
		integer (kind = 4) :: xsrccoord
!> source coordinate y
		integer (kind = 4) :: ysrccoord
!> group coordinate x
		integer (kind = 4) :: xgrpcoord
!> group coordinate y
		integer (kind = 4) :: ygrpcoord
!> coordinate units
		integer (kind = 2) :: coordunits
!> source static correction in milliseconds
		integer (kind = 2) :: srccorr
!> group static correction in milliseconds
		integer (kind = 2) :: grpcorr
!> delay recording time
		integer (kind = 2) :: recdelay
!> number of samples in this trace
		integer (kind = 2) :: numsamp
!> sample interval in microseconds
		integer (kind = 2) :: sampint
!> year data recorded
		integer (kind = 2) :: recyear
!> day of year
		integer (kind = 2) :: doy
!> hour of minute
		integer (kind = 2) :: hod
!> minute of hour
		integer (kind = 2) :: moh
!> second of minute
		integer (kind = 2) :: som
!> time basis code
		integer (kind = 2) :: timecode
!> trace weighting factor
		integer (kind = 2) :: trweight
!> trace value measurement unit
		integer (kind = 2) :: valunits
!> source type/orientation
		integer (kind = 2) :: srctype
	end type traceheader
!
!-------------------------------------------------------------
!> \brief trace_segy object
!
	type trace_segy
		private
		type (traceheader) :: traceheader
		integer (kind = 4), allocatable, dimension(:) :: trace4int
		integer (kind = 2), allocatable, dimension(:) :: trace2int
		real (kind = 4), allocatable, dimension(:) :: traceieee
		real (kind = 4), allocatable, dimension(:) :: traceibm
	end type trace_segy
!
	contains
!-------------------------------------------------------------
!> \brief Fill trace_segy object with data from opened file
!> \param this		trace_segy object
!> \param numtr		number of traces
!> \param numsamp	number of samples
!> \param sampformat	sampling format code
!> \param bytespersamp	bytes per sample
!> \param fixedlen	indicates whether or not all traces have same sample interval and number of samples (1: yes, 0: no)
!> \param skipext	indicates whether or not file contains extended file headers (x: x extended file headers, 0: no)
!
	subroutine readTraceSegy(this,numtr,numsamp,sampformat,bytespersamp,fixedlen,skipext)
!
		type (trace_segy), allocatable, dimension(:) :: this
		integer (kind = 2) :: numtr
		integer (kind = 2) :: numsamp
		integer (kind = 2) :: sampformat
		integer (kind = 2) :: bytespersamp
		integer (kind = 2) :: fixedlen
		integer (kind = 2) :: skipext
		integer (kind = 1) :: status
		integer (kind = 2) :: i,j
		integer (kind = 8) :: alltrlen,tmp
!
		allocate (this(numtr), stat=status)
		if (status /= 0) then
			stop "mod traceSegy: unable to allocate memory"
		endif
!
		alltrlen = 0
		do i=1, numtr
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+1) this(i)%traceheader%trnr
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+29) this(i)%traceheader%trid
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+41) this(i)%traceheader%recele
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+45) this(i)%traceheader%surfele
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+49) this(i)%traceheader%srcdepth
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+69) this(i)%traceheader%elepow
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+71) this(i)%traceheader%coordpow
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+73) this(i)%traceheader%xsrccoord
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+77) this(i)%traceheader%ysrccoord
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+81) this(i)%traceheader%xgrpcoord
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+85) this(i)%traceheader%ygrpcoord
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+89) this(i)%traceheader%coordunits
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+99) this(i)%traceheader%srccorr
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+101) this(i)%traceheader%grpcorr
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+109) this(i)%traceheader%recdelay
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+115) this(i)%traceheader%numsamp
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+117) this(i)%traceheader%sampint
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+157) this(i)%traceheader%recyear
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+159) this(i)%traceheader%doy
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+161) this(i)%traceheader%hod
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+163) this(i)%traceheader%moh
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+165) this(i)%traceheader%som
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+167) this(i)%traceheader%timecode
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+169) this(i)%traceheader%trweight
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+203) this(i)%traceheader%valunits
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+215) this(i)%traceheader%timepow
			read(11, pos=3600+skipext+(i-1)*240+alltrlen+217) this(i)%traceheader%srctype
!
			if (sampformat == 1) then
				if (fixedlen == 1) then
					allocate (this(i)%traceibm(numsamp), stat=status)
				else
					allocate (this(i)%traceibm(this(i)%traceheader%numsamp), stat=status)
				endif
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			elseif (sampformat == 2) then
				if (fixedlen == 1) then
					allocate (this(i)%trace4int(numsamp), stat=status)
				else
					allocate (this(i)%trace4int(this(i)%traceheader%numsamp), stat=status)
				endif
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			elseif (sampformat == 3) then
				if (fixedlen == 1) then
					allocate (this(i)%trace2int(numsamp), stat=status)
				else
					allocate (this(i)%trace2int(this(i)%traceheader%numsamp), stat=status)
				endif
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			elseif (sampformat == 5) then
				if (fixedlen == 1) then
					allocate (this(i)%traceieee(numsamp), stat=status)
				else
					allocate (this(i)%traceieee(this(i)%traceheader%numsamp), stat=status)
				endif
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			endif
!
			if (sampformat == 1) then
				if (fixedlen == 1) then
					do j=1, numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceibm(j)
						this(i)%traceibm(j) = ibm2ieeeBigEndianTraceSegy(this(i)%traceibm(j))
					end do
				else
					do j=1, this(i)%traceheader%numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceibm(j)
						this(i)%traceibm(j) = ibm2ieeeBigEndianTraceSegy(this(i)%traceibm(j))
					end do
				endif
			elseif (sampformat == 2) then
				if (fixedlen == 1) then
					do j=1, numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%trace4int(j)
					end do
				else
					do j=1, this(i)%traceheader%numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%trace4int(j)
					end do
				endif
			elseif (sampformat == 3) then
				if (fixedlen == 1) then
					do j=1, numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%trace2int(j)
					end do
				else
					do j=1, this(i)%traceheader%numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%trace2int(j)
					end do
				endif
			elseif (sampformat == 5) then
				if (fixedlen == 1) then
					do j=1, numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceieee(j)
					end do
				else
					do j=1, this(i)%traceheader%numsamp
						read(11, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceieee(j)
					end do
				endif
			endif
			if (fixedlen == 1) then
				tmp = numsamp
				tmp = tmp*bytespersamp
				alltrlen = alltrlen + tmp
			else
				tmp = this(i)%traceheader%numsamp
				tmp = tmp*bytespersamp
				alltrlen = alltrlen + tmp
			endif
		end do
!
	end subroutine readTraceSegy
!-------------------------------------------------------------
!> \brief Allocate memory for the necessary amount of trace_segy objects
!> \param this		trace_segy object
!> \param numtr		amount of traces
!
	subroutine allocTraceSegy(this,numtr)
!
	type (trace_segy), allocatable, dimension(:), intent(inout) :: this
	integer (kind = 2), intent(in) :: numtr
!
	integer (kind = 1) :: status
!
	if(.not. allocated(this)) then
		allocate (this(numtr), stat=status)
		if (status /= 0) then
			stop "mod traceSegy: unable to allocate memory"
		endif
	endif
!
	end subroutine allocTraceSegy
!-------------------------------------------------------------
!> \brief Fill trace_segy object with given data
!> \param this		data_segy object
!> \param trnr		tracenumber
!> \param trace		array consisting of tracevalues
!> \param trid		trace identification code
!> \param recele	receiver group elevation
!> \param surfele	surface elevation at source
!> \param srcdepth	source depth below surface (a positive number)
!> \param elepow	scalar to be applied to elevations to give real value (positiv->multiplier,negativ->divisor)
!> \param coordpow	scalar to be applied to coordinates to give real value (positiv->multiplier,negativ->divisor)
!> \param timepow	scalar to be applied to times to give real value (positiv->multiplier,negativ->divisor)
!> \param xsrccoord	source coordinate: X
!> \param ysrccoord	source coordinate: Y
!> \param xgrpcoord	group coordinate: X
!> \param ygrpcoord	group coordinate: Y
!> \param coordunits	coordinate units
!> \param srccorr	source static correction in milliseconds
!> \param grpcorr	group static correction in milliseconds
!> \param recdelay	delay recording time
!> \param sampint	sampling interval (indispensable if fixedlen was set to 0 in binary header)
!> \param recyear	year data recorded
!> \param doy		day of year
!> \param hod		hour of day
!> \param moh		minute of hour
!> \param som		second of minute
!> \param timecode	time basis code
!> \param trweight	trace weighting factor
!> \param valunits	trace value measurement unit
!> \param srctype	source type
!> \param sampformat	sampling format code
!
	subroutine buildTraceSegy(this,trnr,trace,trid,recele,surfele,srcdepth,elepow,coordpow,timepow,xsrccoord,ysrccoord,xgrpcoord,&
				& ygrpcoord,coordunits,srccorr,grpcorr,recdelay,numsamp,sampint,recyear,doy,hod,moh,som,timecode,&
				& trweight,valunits,srctype,sampformat)
!---------------------------------------------------------
		type (trace_segy), intent(inout) :: this
		integer (kind = 4), intent(in) :: trnr
		integer (kind = 2), intent(in) :: trid
		integer (kind = 4), intent(in) :: recele
		integer (kind = 4), intent(in) :: surfele
		integer (kind = 4), intent(in) :: srcdepth
		integer (kind = 2), intent(in) :: elepow
		integer (kind = 2), intent(in) :: coordpow
		integer (kind = 2), intent(in) :: timepow
		integer (kind = 4), intent(in) :: xsrccoord
		integer (kind = 4), intent(in) :: ysrccoord
		integer (kind = 4), intent(in) :: xgrpcoord
		integer (kind = 4), intent(in) :: ygrpcoord
		integer (kind = 2), intent(in) :: coordunits
		integer (kind = 2), intent(in) :: srccorr
		integer (kind = 2), intent(in) :: grpcorr
		integer (kind = 2), intent(in) :: recdelay
		integer (kind = 2), intent(in) :: numsamp
		integer (kind = 2), intent(in) :: sampint
		integer (kind = 2), intent(in) :: recyear
		integer (kind = 2), intent(in) :: doy
		integer (kind = 2), intent(in) :: hod
		integer (kind = 2), intent(in) :: moh
		integer (kind = 2), intent(in) :: som
		integer (kind = 2), intent(in) :: timecode
		integer (kind = 2), intent(in) :: trweight
		integer (kind = 2), intent(in) :: valunits
		integer (kind = 2), intent(in) :: srctype
		integer (kind = 2), intent(in) :: sampformat
		real (kind = 4), intent(in), dimension(:) :: trace
		integer (kind = 1) :: status
		integer (kind = 2) :: i
!
		this%traceheader%trnr = trnr
		this%traceheader%trid = trid
		this%traceheader%recele = recele
		this%traceheader%surfele = surfele
		this%traceheader%srcdepth = srcdepth
		this%traceheader%elepow = elepow
		this%traceheader%coordpow = coordpow
		this%traceheader%timepow = timepow
		this%traceheader%xsrccoord = xsrccoord
		this%traceheader%ysrccoord = ysrccoord
		this%traceheader%xgrpcoord = xgrpcoord
		this%traceheader%ygrpcoord = ygrpcoord
		this%traceheader%coordunits = coordunits
		this%traceheader%srccorr = srccorr
		this%traceheader%grpcorr = grpcorr
		this%traceheader%recdelay = recdelay
		this%traceheader%numsamp = numsamp
		this%traceheader%sampint = sampint
		this%traceheader%recyear = recyear
		this%traceheader%doy = doy
		this%traceheader%hod = hod
		this%traceheader%moh = moh
		this%traceheader%som = som
		this%traceheader%timecode = timecode
		this%traceheader%trweight = trweight
		this%traceheader%valunits = valunits
		this%traceheader%srctype = srctype
!
		if (sampformat == 5) then
			if(.not. allocated(this%traceieee)) then
				allocate (this%traceieee(this%traceheader%numsamp), stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			else
				deallocate(this%traceieee, stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to deallocate memory"
				endif
				allocate (this%traceieee(this%traceheader%numsamp), stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to allocate memory"
				endif
			endif
		else
			stop "mod traceSegy: can not process this sampling code"
		endif
!
		if (sampformat == 5) then
			do i=1, this%traceheader%numsamp
				this%traceieee(i) = trace(i)
			enddo
		else
			stop "mod traceSegy: can not process this sampling code"
		endif
!
	end subroutine buildTraceSegy
!-------------------------------------------------------------
!> \brief Write content of trace_segy object to opened file
!> \param this		trace_segy object
!> \param numtr		number of traces
!> \param numsamp	number of samples
!> \param sampformat	sampling format code
!> \param bytespersamp	bytes per sample
!> \param skipext	indicates whether or not file contains extended file headers (x: x extended file headers, 0: no)
!> \param fixedlen	indicates whether or not all traces have same sample interval and number of samples (1: yes, 0: no)
!
	subroutine writeTraceSegy(this,lu,numtr,numsamp,sampformat,bytespersamp,skipext,fixedlen)
!---------------------------------------------------------
		type (trace_segy), dimension(:), intent(in) :: this
		integer :: lu
		integer (kind = 2), intent(in) :: numtr
		integer (kind = 2), intent(in) :: numsamp
		integer (kind = 2), intent(in) :: sampformat
		integer (kind = 2), intent(in) :: bytespersamp
		integer (kind = 2), intent(in) :: skipext
		integer (kind = 2), intent(in) :: fixedlen
		integer (kind = 2) :: i,j
		integer (kind = 8) :: alltrlen,tmp
!
		alltrlen = 0
		do i=1, numtr
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+1) this(i)%traceheader%trnr
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+29) this(i)%traceheader%trid
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+41) this(i)%traceheader%recele
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+45) this(i)%traceheader%surfele
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+49) this(i)%traceheader%srcdepth
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+69) this(i)%traceheader%elepow
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+71) this(i)%traceheader%coordpow
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+73) this(i)%traceheader%xsrccoord
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+77) this(i)%traceheader%ysrccoord
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+81) this(i)%traceheader%xgrpcoord
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+85) this(i)%traceheader%ygrpcoord
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+89) this(i)%traceheader%coordunits
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+99) this(i)%traceheader%srccorr
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+101) this(i)%traceheader%grpcorr
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+109) this(i)%traceheader%recdelay
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+115) this(i)%traceheader%numsamp
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+117) this(i)%traceheader%sampint
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+157) this(i)%traceheader%recyear
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+159) this(i)%traceheader%doy
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+161) this(i)%traceheader%hod
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+163) this(i)%traceheader%moh
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+165) this(i)%traceheader%som
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+167) this(i)%traceheader%timecode
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+169) this(i)%traceheader%trweight
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+203) this(i)%traceheader%valunits
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+215) this(i)%traceheader%timepow
			write(lu, pos=3600+skipext+(i-1)*240+alltrlen+217) this(i)%traceheader%srctype

			if (sampformat == 5) then
				if (fixedlen == 1) then
					do j=1, numsamp
					   write(lu,pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceieee(j)
					end do
				else
					do j=1, this(i)%traceheader%numsamp
					   write(lu, pos=3600+skipext+(i)*240+alltrlen+bytespersamp*(j-1)+1) this(i)%traceieee(j)
					end do
				endif
			endif
			if (fixedlen == 1) then
				tmp = numsamp
				tmp = tmp*bytespersamp
				alltrlen = alltrlen + tmp
			else
				tmp = this(i)%traceheader%numsamp
				tmp = tmp*bytespersamp
				alltrlen = alltrlen + tmp
			endif
		enddo
!
	end subroutine writeTraceSegy
!-------------------------------------------------------------
!> \brief Dealloc memory associated with trace_segy object
!> \param this		trace_segy object
!> \param numtr		number of traces
!> \param sampformat	sampling format code
!
	subroutine deallocTraceSegy(this,numtr,sampformat)
!---------------------------------------------------------
		type (trace_segy), allocatable, dimension(:) :: this
		integer (kind = 2), intent(in) :: numtr
		integer (kind = 2), intent(in) :: sampformat
!
		integer (kind = 1) :: status
		integer (kind = 2) :: i
!
		if (sampformat == 1) then
			do i=1, numtr - 1
				deallocate (this(i)%traceibm, stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to deallocate memory"
				endif
			enddo
		endif
		if (sampformat == 2) then
			do i=1, numtr - 1
				deallocate (this(i)%trace4int, stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to deallocate memory"
				endif
			enddo
		endif
		if (sampformat == 3) then
			do i=1, numtr - 1
				deallocate (this(i)%trace2int, stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to deallocate memory"
				endif
			enddo
		endif
		if (sampformat == 5) then
			do i=1, numtr - 1
				deallocate (this(i)%traceieee, stat=status)
				if (status /= 0) then
					stop "mod traceSegy: unable to deallocate memory"
				endif
			enddo
		endif
!
		deallocate (this, stat=status)
		if (status /= 0) then
			stop "mod traceSegy: unable to deallocate memory"
		endif
!
	end subroutine deallocTraceSegy
!-------------------------------------------------------------
!> \brief Convert values which are read incorrectly due to the different bitordering of ibm floating point numbers
!> \param val		value to be converted
!
real (kind=4) function ibm2ieeeBigEndianTraceSegy(val)
	real (kind=4) :: val,Mantisse
	integer (kind=4) :: val_int,Fraction
	integer (kind=4) :: sign
	integer (kind=2) :: A,B,exponent,e
	integer (kind=1) :: i

	A = 16
	B = 64
	
	val_int = transfer(val,val_int)

!	bit 31 ist sign
	if (btest(val_int,31)) then
		sign = -1
	else
		sign = +1
	endif
!	bit 30-24 ist exponent
	exponent = ibits(val_int,24,7)
	e = exponent - B
!	bit 23-00 sind Fraction (als linksbündigen Binärbruch zu lesen, d.h. z.B. 101 = 1*1/2 + 0*1/4 + 1*1/8 = 5/8)
	Fraction = ibits(val_int,0,23)
	Mantisse = 0.;i=0
	do while(i.le.23)
		if (btest(val_int,i)) then
			Mantisse = Mantisse + 1./(2.**(24-i))
		endif
		i=i+1
	enddo

	ibm2ieeeBigEndianTraceSegy = sign*Mantisse*(float(A)**e)
end function ibm2ieeeBigEndianTraceSegy
!-------------------------------------------------------------
!> \brief Get trace number
!> \param this		trace_segy object
!
	integer (kind = 4) function getTrNrTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getTrNrTraceSegy = this%traceheader%trnr
	end function getTrNrTraceSegy
!-------------------------------------------------------------
!> \brief Get trace identification code
!> \param this		trace_segy object
!
	integer (kind = 4) function getTrIdTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getTrIdTraceSegy = this%traceheader%trid
	end function getTrIdTraceSegy
!-------------------------------------------------------------
!> \brief Get receiver group elevation
!> \param this		trace_segy object
!
	integer (kind = 4) function getRecEleTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.elepow.this > 0) then
			getRecEleTraceSegy = (this%traceheader%recele)*(.elepow.this)
		elseif (.elepow.this < 0) then
			getRecEleTraceSegy = (this%traceheader%recele)/(.elepow.this)*(-1)
		else
			getRecEleTraceSegy = this%traceheader%recele
		endif	
	end function getRecEleTraceSegy
!-------------------------------------------------------------
!> \brief Get surface elevation at source
!> \param this		trace_segy object
!
	integer (kind = 4) function getSurfEleTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.elepow.this > 0) then
			getSurfEleTraceSegy = (this%traceheader%surfele)*(.elepow.this)
		elseif (.elepow.this < 0) then
			getSurfEleTraceSegy = (this%traceheader%surfele)/(.elepow.this)*(-1)
		else
			getSurfEleTraceSegy = this%traceheader%surfele
		endif
	end function getSurfEleTraceSegy
!-------------------------------------------------------------
!> \brief Get source depth below surface (a positive number)
!> \param this		trace_segy object
!
	integer (kind = 4) function getSrcDepthTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.elepow.this > 0) then
			getSrcDepthTraceSegy = (this%traceheader%srcdepth)*(.elepow.this)
		elseif (.elepow.this < 0) then
			getSrcDepthTraceSegy = (this%traceheader%srcdepth)/(.elepow.this)*(-1)
		else
			getSrcDepthTraceSegy = this%traceheader%srcdepth
		endif
	end function getSrcDepthTraceSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to elevations to give real value (positiv->multiplier,negativ->divisor)
!> \param this		trace_segy object
!
	integer (kind = 4) function getElePowTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getElePowTraceSegy = this%traceheader%elepow
	end function getElePowTraceSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to coordinates to give real value (positiv->multiplier,negativ->divisor)
!> \param this		trace_segy object
!
	integer (kind = 4) function getCoordPowTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getCoordPowTraceSegy = this%traceheader%coordpow
	end function getCoordPowTraceSegy
!-------------------------------------------------------------
!> \brief Get scalar to be applied to times to give real value (positiv->multiplier,negativ->divisor)
!> \param this		trace_segy object
!
	integer (kind = 4) function getTimePowTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getTimePowTraceSegy = this%traceheader%timepow
	end function getTimePowTraceSegy
!-------------------------------------------------------------
!> \brief Get source coordinate: X
!> \param this		trace_segy object
!
	real (kind = 4) function getXSrcCoordTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.coordpow.this > 0) then
			getXSrcCoordTraceSegy = real((this%traceheader%xsrccoord),4)*real(.coordpow.this,4)
		elseif (.coordpow.this < 0) then
			getXSrcCoordTraceSegy = real((this%traceheader%xsrccoord),4)/real(-1*(.coordpow.this),4)
		else
			getXSrcCoordTraceSegy = real(this%traceheader%xsrccoord,4)
		endif
	end function getXSrcCoordTraceSegy
!-------------------------------------------------------------
!> \brief Get source coordinate: Y
!> \param this		trace_segy object
!
	real (kind = 4) function getYSrcCoordTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.coordpow.this > 0) then
			getYSrcCoordTraceSegy = real((this%traceheader%ysrccoord),4)*real(.coordpow.this,4)
		elseif (.coordpow.this < 0) then
			getYSrcCoordTraceSegy = real((this%traceheader%ysrccoord),4)/real(-1*(.coordpow.this),4)
		else
			getYSrcCoordTraceSegy = real(this%traceheader%ysrccoord,4)
		endif
	end function getYSrcCoordTraceSegy
!-------------------------------------------------------------
!> \brief Get group coordinate: X
!> \param this		trace_segy object
!
	real (kind = 4) function getXGrpCoordTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.coordpow.this > 0) then
			getXGrpCoordTraceSegy = real((this%traceheader%xgrpcoord),4)*real(.coordpow.this,4)
		elseif (.coordpow.this < 0) then
			getXGrpCoordTraceSegy = real((this%traceheader%xgrpcoord),4)/real(-1*(.coordpow.this),4)
		else
			getXGrpCoordTraceSegy = real(this%traceheader%xgrpcoord,4)
		endif
	end function getXGrpCoordTraceSegy
!-------------------------------------------------------------
!> \brief Get group coordinate: Y
!> \param this		trace_segy object
!
	real (kind = 4) function getYGrpCoordTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.coordpow.this > 0) then
			getYGrpCoordTraceSegy = real((this%traceheader%ygrpcoord),4)*real(.coordpow.this,4)
		elseif (.coordpow.this < 0) then
			getYGrpCoordTraceSegy = real((this%traceheader%ygrpcoord),4)/real(-1*(.coordpow.this),4)
		else
			getYGrpCoordTraceSegy = real(this%traceheader%ygrpcoord,4)
		endif
	end function getYGrpCoordTraceSegy
!-------------------------------------------------------------
!> \brief Get coordinate units
!> \param this		trace_segy object
!
	integer (kind = 4) function getCoordUnitsTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getCoordUnitsTraceSegy = this%traceheader%coordunits
	end function getCoordUnitsTraceSegy
!-------------------------------------------------------------
!> \brief Get source static correction in milliseconds
!> \param this		trace_segy object
!
	integer (kind = 4) function getSrcCorrTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.timepow.this > 0) then
			getSrcCorrTraceSegy = (this%traceheader%srccorr)*(.coordpow.this)
		elseif (.timepow.this < 0) then
			getSrcCorrTraceSegy = (this%traceheader%srccorr)/(.coordpow.this)*(-1)
		else
			getSrcCorrTraceSegy = this%traceheader%srccorr
		endif
	end function getSrcCorrTraceSegy
!-------------------------------------------------------------
!> \brief Get group static correction in milliseconds
!> \param this		trace_segy object
!
	integer (kind = 4) function getGrpCorrTraceSegy(this)
		type(trace_segy), intent(in) :: this
		if (.timepow.this > 0) then
			getGrpCorrTraceSegy = (this%traceheader%grpcorr)*(.coordpow.this)
		elseif (.timepow.this < 0) then
			getGrpCorrTraceSegy = (this%traceheader%grpcorr)/(.coordpow.this)*(-1)
		else
			getGrpCorrTraceSegy = this%traceheader%grpcorr
		endif
	end function getGrpCorrTraceSegy
!-------------------------------------------------------------
!> \brief Get delay recording time
!> \param this		trace_segy object
!
	integer (kind = 4) function getRecDelayTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getRecDelayTraceSegy = this%traceheader%recdelay
	end function getRecDelayTraceSegy
!-------------------------------------------------------------
!> \brief Get sampling interval
!> \param this		trace_segy object
!
	integer (kind = 4) function getSampIntTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getSampIntTraceSegy = this%traceheader%sampint
	end function getSampIntTraceSegy
!-------------------------------------------------------------
!> \brief Get number of samples per trace
!> \param this		trace_segy object
!
	integer (kind = 4) function getNumSampTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getNumSampTraceSegy = this%traceheader%numsamp
	end function getNumSampTraceSegy
!-------------------------------------------------------------
!> \brief Get year data recorded
!> \param this		trace_segy object
!
	integer (kind = 4) function getRecYearTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getRecYearTraceSegy = this%traceheader%recyear
	end function getRecYearTraceSegy
!-------------------------------------------------------------
!> \brief Get day of year
!> \param this		trace_segy object
!
	integer (kind = 4) function getDoyTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getDoyTraceSegy = this%traceheader%doy
	end function getDoyTraceSegy
!-------------------------------------------------------------
!> \brief Get hour of day
!> \param this		trace_segy object
!
	integer (kind = 4) function getHodTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getHodTraceSegy = this%traceheader%hod
	end function getHodTraceSegy
!-------------------------------------------------------------
!> \brief Get minute of hour
!> \param this		trace_segy object
!
	integer (kind = 4) function getMohTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getMohTraceSegy = this%traceheader%moh
	end function getMohTraceSegy
!-------------------------------------------------------------
!> \brief Get seconds of minute
!> \param this		trace_segy object
!
	integer (kind = 4) function getSomTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getSomTraceSegy = this%traceheader%som
	end function getSomTraceSegy
!-------------------------------------------------------------
!> \brief Get time basis code
!> \param this		trace_segy object
!
	integer (kind = 4) function getTimecodeTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getTimecodeTraceSegy = this%traceheader%timecode
	end function getTimecodeTraceSegy
!-------------------------------------------------------------
!> \brief Get trace weighting factor
!> \param this		trace_segy object
!
	integer (kind = 4) function getTrWeightTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getTrWeightTraceSegy = this%traceheader%trweight
	end function getTrWeightTraceSegy
!-------------------------------------------------------------
!> \brief Get trace value measurement unit
!> \param this		trace_segy object
!
	integer (kind = 4) function getValUnitsTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getValUnitsTraceSegy = this%traceheader%valunits
	end function getValUnitsTraceSegy
!-------------------------------------------------------------
!> \brief Get source type
!> \param this		trace_segy object
!
	integer (kind = 4) function getSrcTypeTraceSegy(this)
		type(trace_segy), intent(in) :: this
		getSrcTypeTraceSegy = this%traceheader%srctype
	end function getSrcTypeTraceSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 2; you will receive a pointer)
!> \param this		trace_segy object
!> \param n		trace number
!
	function getTrace4IntTraceSegy(this,n) result(p)
		type(trace_segy), dimension(:), target, intent(in) :: this
		integer (kind = 4), intent(in) :: n
		integer (kind = 4), dimension(:), pointer :: p
		p => this(n)%trace4int
	end function getTrace4IntTraceSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 3; you will receive a pointer)
!> \param this		trace_segy object
!> \param n		trace number
!

	function getTrace2IntTraceSegy(this,n) result(p)
		type(trace_segy), dimension(:), target, intent(in) :: this
		integer (kind = 4), intent(in) :: n
		integer (kind = 2), dimension(:), pointer :: p
		p => this(n)%trace2int
	end function getTrace2IntTraceSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 5; you will receive a pointer)
!> \param this		trace_segy object
!> \param n		trace number
!
	function getTraceIeeeTraceSegy(this,n) result(p)
		type(trace_segy), dimension(:), target, intent(in) :: this
		integer (kind = 4), intent(in) :: n
		real (kind = 4), dimension(:), pointer :: p
		p => this(n)%traceieee
	end function getTraceIeeeTraceSegy
!-------------------------------------------------------------
!> \brief Get trace values (use if sampling format code is equal to 1; you will receive a pointer)
!> \param this		trace_segy object
!> \param n		trace number
!
	function getTraceIbmTraceSegy(this,n) result(p)
		type(trace_segy), dimension(:), target, intent(in) :: this
		integer (kind = 4), intent(in) :: n
		real (kind = 4), dimension(:), pointer :: p
		p => this(n)%traceibm
	end function getTraceIbmTraceSegy
!-------------------------------------------------------------
end module traceSegy
