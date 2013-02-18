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
module binarySegy
!
	implicit none
!
	interface operator (.numtr.); module procedure getNumTrBinarySegy; end interface
	interface operator (.numtren.); module procedure getNumTrEnBinarySegy; end interface
	interface operator (.gensampint.); module procedure getGenSampIntBinarySegy; end interface
	interface operator (.gennumsamp.); module procedure getGenNumSampBinarySegy; end interface
	interface operator (.sampformat.); module procedure getSampFormatBinarySegy; end interface
	interface operator (.enfold.); module procedure getEnFoldBinarySegy; end interface
	interface operator (.trsort.); module procedure getTrSortBinarySegy; end interface
	interface operator (.meassys.); module procedure getMeasSysBinarySegy; end interface
	interface operator (.sigpol.); module procedure getSigPolBinarySegy; end interface
	interface operator (.revnr.); module procedure getRevNrBinarySegy; end interface
	interface operator (.fixedlen.); module procedure getFixedLenBinarySegy; end interface
	interface operator (.numext.); module procedure getNumExtBinarySegy; end interface
	interface operator (.skipext.); module procedure getSkipExtBinarySegy; end interface
	interface operator (.bytespersamp.); module procedure getBytesPerSampBinarySegy; end interface
!
!-------------------------------------------------------------
!> \brief binary_segy object
!
	type binary_segy
		private
!> number of traces in file
		integer (kind = 2) :: numtr
!> number of auxiliary traces per ensemble
		integer (kind = 2) :: numtren
!> generic sample interval in microseconds (applies to all traces if fixed length trace flag is set)
		integer (kind = 2) :: gensampint
!> generic number of samples per trace (applies to all traces if fixed length trace flag is set)
		integer (kind = 2) :: gennumsamp
!> sampling format code
		integer (kind = 2) :: sampformat
!> ensemble fold - the expected number of data traces per trace ensemble (e.g. CMP fold)
		integer (kind = 2) :: enfold
!> trace sorting code
		integer (kind = 2) :: trsort
!> measurement system
		integer (kind = 2) :: meassys
!> impulse signal polarity
		integer (kind = 2) :: sigpol
!> SEG Y format revision number
		integer (kind = 2) :: revnr
!> fixed length trace flag
		integer (kind = 2) :: fixedlen
!> number of extended textual file header records
		integer (kind = 2) :: skipext
	end type binary_segy
!
	contains
!-------------------------------------------------------------
!> \brief Fill binary_segy object data from opened file
!> \param this		binary_segy object
!
	subroutine readBinarySegy(this)
!
		type (binary_segy), intent(inout) :: this
		integer (kind = 1) :: status
		integer (kind = 2) :: i
		integer (kind = 2) :: numsamp,dummy
		integer (kind = 8) :: alltrlen,tmp
!
		read(11, pos=3213) this%numtren
		read(11, pos=3217) this%gensampint
		read(11, pos=3221) this%gennumsamp
		read(11, pos=3225) this%sampformat
		read(11, pos=3227) this%enfold
		read(11, pos=3229) this%trsort
		read(11, pos=3255) this%meassys
		read(11, pos=3257) this%sigpol
		read(11, pos=3501) this%revnr
		read(11, pos=3503) this%fixedlen
		read(11, pos=3505) this%skipext
!
		this%numtr = 0
		alltrlen = 0
		i = 0
		if (.fixedlen.this == 1) then
			numsamp = .gennumsamp.this
			do
				tmp = .bytespersamp.this
				tmp = numsamp*tmp
				alltrlen = alltrlen+tmp
				read(11, pos=3600+(.skipext.this)+(i+1)*240+alltrlen+1,iostat=status) dummy
				i=i+1;this%numtr=this%numtr+1
				if (status .ne. 0) exit
			enddo
		else
			do
				read(11, pos=3600+(.skipext.this)+(i)*240+alltrlen+115) numsamp
				tmp = .bytespersamp.this
				tmp = numsamp*tmp
				alltrlen = alltrlen+tmp
				read(11, pos=3600+(.skipext.this)+(i+1)*240+alltrlen+1,iostat=status) numsamp
				i=i+1;this%numtr=this%numtr+1
				if (status .ne. 0) exit
			enddo
		endif
!
	end subroutine readBinarySegy
!-------------------------------------------------------------
!> \brief Fill binary_segy object with given data
!> \param this		binary_segy object
!> \param numtr		number of traces
!> \param numtren	number of traces per ensemble
!> \param gensampint	sampling interval
!> \param gennumsamp	number of samples per trace
!> \param revnr		segy revision number
!> \param fixedlen	indicates whether or not all traces have same sample interval and number of samples (1: yes, 0: no)
!> \param skipext	indicates whether or not file contains extended file headers (x: x extended file headers, 0: no)
!> \param enfold	ensemble fold
!> \param trsort	trace sorting code
!> \param meassys	measurement system
!> \param sigpol	impulse signal polarity
!
	subroutine buildBinarySegy(this,numtr,fixedlen,skipext,numtren,gensampint,gennumsamp,revnr,enfold,trsort,meassys,sigpol)
!
		type (binary_segy), intent(inout) :: this
		integer (kind = 2), intent(in) :: numtr
		integer (kind = 2), intent(in) :: numtren
		integer (kind = 2), intent(in) :: gensampint
		integer (kind = 2), intent(in) :: gennumsamp
		integer (kind = 2), intent(in) :: revnr
		integer (kind = 2), intent(in) :: fixedlen
		integer (kind = 2), intent(in) :: skipext
		integer (kind = 2), intent(in) :: enfold
		integer (kind = 2), intent(in) :: trsort
		integer (kind = 2), intent(in) :: meassys
		integer (kind = 2), intent(in) :: sigpol
!
		this%numtr = numtr
		this%numtren = numtren
		this%gensampint = gensampint
		this%gennumsamp = gennumsamp
		this%revnr = revnr
		this%sampformat = 5
		this%enfold = enfold
		this%trsort = trsort
		this%meassys = meassys
		this%sigpol = sigpol
		this%fixedlen = fixedlen
		this%skipext = skipext
!
	end subroutine buildBinarySegy
!-------------------------------------------------------------
!> \brief Write content of binary_segy object to opened file
!> \param this		binary_segy object
!
	subroutine writeBinarySegy(this,lu)
!
		type (binary_segy), intent(in) :: this
		integer :: lu
!
		write(lu, pos=3213) this%numtren
		write(lu, pos=3217) this%gensampint
		write(lu, pos=3221) this%gennumsamp
		write(lu, pos=3225) this%sampformat
		write(lu, pos=3227) this%enfold
		write(lu, pos=3229) this%trsort
		write(lu, pos=3255) this%meassys
		write(lu, pos=3501) this%revnr
		write(lu, pos=3503) this%fixedlen
		write(lu, pos=3505) this%skipext
!
	end subroutine writeBinarySegy
!-------------------------------------------------------------
!> \brief Get number of traces
!> \param this		binary_segy object
!
	integer (kind = 2) function getNumTrBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getNumTrBinarySegy = this%numtr
	end function getNumTrBinarySegy
!-------------------------------------------------------------
!> \brief Get number of traces per ensemble
!> \param this		binary_segy object
!
	integer (kind = 2) function getNumTrEnBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getNumTrEnBinarySegy = this%numtren
	end function getNumTrEnBinarySegy
!-------------------------------------------------------------
!> \brief Get sampling interval
!> \param this		binary_segy object
!
	integer (kind = 2) function getGenSampIntBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getGenSampIntBinarySegy = this%gensampint
	end function getGenSampIntBinarySegy
!-------------------------------------------------------------
!> \brief Get number of samples per trace
!> \param this		binary_segy object
!
	integer (kind = 2) function getGenNumSampBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getGenNumSampBinarySegy = this%gennumsamp
	end function getGenNumSampBinarySegy
!-------------------------------------------------------------
!> \brief Get sampling format code
!> \param this		binary_segy object
!
	integer (kind = 2) function getSampFormatBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getSampFormatBinarySegy = this%sampformat
	end function getSampFormatBinarySegy
!-------------------------------------------------------------
!> \brief Get ensemble fold
!> \param this		binary_segy object
!
	integer (kind = 2) function getEnFoldBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getEnFoldBinarySegy = this%enfold
	end function getEnFoldBinarySegy
!-------------------------------------------------------------
!> \brief Get trace sorting code
!> \param this		binary_segy object
!
	integer (kind = 2) function getTrSortBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getTrSortBinarySegy = this%trsort
	end function getTrSortBinarySegy
!-------------------------------------------------------------
!> \brief Get measurement system
!> \param this		binary_segy object
!
	integer (kind = 2) function getMeasSysBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getMeasSysBinarySegy = this%meassys
	end function getMeasSysBinarySegy
!-------------------------------------------------------------
!> \brief Get impulse signal polarity
!> \param this		binary_segy object
!
	integer (kind = 2) function getSigPolBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getSigPolBinarySegy = this%sigpol
	end function getSigPolBinarySegy
!-------------------------------------------------------------
!> \brief Get segy format revision number
!> \param this		binary_segy object
!
	integer (kind = 2) function getRevNrBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getRevNrBinarySegy = this%revnr
	end function getRevNrBinarySegy
!-------------------------------------------------------------
!> \brief Get fixed lenght indicator
!> \param this		binary_segy object
!
	integer (kind = 2) function getFixedLenBinarySegy(this)
		type (binary_segy), intent(in) :: this
		getFixedLenBinarySegy = this%fixedlen
	end function getFixedLenBinarySegy
!-------------------------------------------------------------
!> \brief Get number of extended file headers
!> \param this		binary_segy object
!
	integer (kind = 2) function getNumExtBinarySegy(this)
		type (binary_segy), intent(in) :: this
		if (this%skipext .ge. 0) then
			getNumExtBinarySegy = this%skipext
		else
			stop "mod binarySegy: unknown number of extended file headers"
		endif
	end function getNumExtBinarySegy
!-------------------------------------------------------------
!> \brief Get byte shift due to existence of extended file headers
!> \param this		binary_segy object
!
	integer (kind = 2) function getSkipExtBinarySegy(this)
		type (binary_segy), intent(in) :: this
		if (this%skipext .ge. 0) then 
			getSkipExtBinarySegy = this%skipext*3200
		else
			stop "mod binarySegy: unknown number of extended file headers"
		endif
	end function getSkipExtBinarySegy
!-------------------------------------------------------------
!> \brief Get number of bytes per sample
!> \param this		binary_segy object
!
	integer (kind = 2) function getBytesPerSampBinarySegy(this)
		type (binary_segy), intent(in) :: this
		if (this%sampformat == 2 .or. this%sampformat == 5 .or. this%sampformat == 1) then
			getBytesPerSampBinarySegy = 4
		elseif (this%sampformat == 3) then
			getBytesPerSampBinarySegy = 2
		else
			stop "mod binarySegy: can not process this sampling code"
		endif
	end function getBytesPerSampBinarySegy
!-------------------------------------------------------------
end module binarySegy
