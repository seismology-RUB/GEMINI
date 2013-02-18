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
!----------------------------------------------------------------------
!  Module with routine to handle frequency wavenumber spectra
!----------------------------------------------------------------------
 module greenFrequencyWavenumber
	use mathConstants
	use streamAccess
	use wavenumberIntegralMapping
	use errorMessage
	use fourierTransform
	use fileUnitHandler
	use geminiEarthModel
	implicit none
	interface readBodyGreenFrequencyWavenumber
		module procedure readBodyZsGreenFrequencyWavenumber
		module procedure readBodyJrGreenFrequencyWavenumber
	end interface
	interface readParallelBodyGreenFrequencyWavenumber
		module procedure readParallelBodyZsGreenFrequencyWavenumber
		module procedure readParallelBodyJrGreenFrequencyWavenumber
	end interface
	interface operator (.styp.); module procedure getSourceTypeGreenFrequencyWavenumber; end interface
	interface operator (.nfa.); module procedure getNf1GreenFrequencyWavenumber; end interface
	interface operator (.nfb.); module procedure getNf2GreenFrequencyWavenumber; end interface
	interface operator (.nf.); module procedure getNfGreenFrequencyWavenumber; end interface
	interface operator (.nkfmax.); module procedure getNkfmaxGreenFrequencyWavenumber; end interface
	interface operator (.nwnmax.); module procedure getNwnmaxGreenFrequencyWavenumber; end interface
	interface operator (.nwn.); module procedure getNwnSelectedGreenFrequencyWavenumber; end interface
	interface operator (.nnod.); module procedure getNnodGreenFrequencyWavenumber; end interface
	interface operator (.df.); module procedure getDfGreenFrequencyWavenumber; end interface
	interface operator (.dwn.); module procedure getDwnGreenFrequencyWavenumber; end interface
	interface operator (.sigma.); module procedure getSigmaGreenFrequencyWavenumber; end interface
	interface operator (.rearth.); module procedure getRerdeGreenFrequencyWavenumber; end interface
	interface operator (.ze.); module procedure getZeGreenFrequencyWavenumber; end interface
	interface operator (.gem.); module procedure getGeminiEarthMomdelGreenFrequencyWavenumber; end interface
	interface dealloc; module procedure deallocGreenFrequencyWavenumber; end interface
	type green_frequency_wavenumber
		private
		type (gemini_earth_model) :: gem                 !< model parameters for Green function
		integer :: istyp                                 !< source type (0 = force, 1 = moment)
		integer :: nf1                                   !< index of first frequency
		integer :: nf2                                   !< index of last frequency
		integer :: nkfmax                                !< size of frequency-wavenumber spectra
		integer :: nwnmax                                !< max number of wavenumbers
		integer :: nsrcsph                               !< number of source jumps spheroidal
		integer :: nsrctor                               !< number of source jumps toroidal
		integer, dimension(6) :: dsvmask_sph             !< mask of available DSV-components within (U,R,V,S,UP,VP)
		integer, dimension(3) :: dsvmask_tor             !< mask of available DSV-components within (W,T,WP)
		integer :: derivflag                             !< flag whether gf-derivatives have been calculated
		integer :: dsvstep                               !< stepping in DSVs, either (URVS,WT) or (UV,W)
		integer, dimension(9,4) :: ksp                   !< mapping fom isp,jsp to 1D spectra count
		integer, dimension(:), pointer :: nwn            !< number of wavenumbers vs frequency
		integer :: numtasks                              !< number of processes by which fk-spectrum was calculated
		real :: df                                       !< frequency spacing
		real :: dwn                                      !< wavenumber spacing
		real :: sigma                                    !< imaginary part of omega
		real :: rearth                                   !< earth radius
		real :: ze                                       !< receiver depth OR source depth
		double precision :: rcur                         !< radius for which spectra are currently read in
		complex, dimension(:,:), pointer :: greenkfsp    !< array for fk-spectra for current radius (rcur)
		type (file_stream_access) :: fda                 !< DA-file objectes
		type (group_stream_access) :: root               !< DA root group objects
		type (file_stream_access), dimension(:), pointer :: fda_par => null()    !< DA file objects for parallel files
		type (group_stream_access), dimension(:), pointer :: root_par => null()  !< DA root group objects
	end type
!
 contains
!-----------------------------------------------------------------------
!  read frequency wavenumber spectra written by sequential gfdsvrkf
!
	function readHeaderGreenFrequencyWavenumber(this,fuh,basename,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	type (file_unit_handler) :: fuh
	integer :: nsp
	logical, optional :: vsflag
	character (len=*) :: basename
	type (error_message) :: errmsg
!
	call new(errmsg,'readHeaderGreenFrequencyWavenumber')
	errmsg = readHeaderBasicGreenFrequencyWavenumber(this,fuh,trim(basename),vsflag)
	if (.level.errmsg == 2) then
		call addTraceErrorMessage(errmsg,'readHeaderGreenFrequencyWavenumber')
		return
	endif
	this%numtasks = 1
!
!  allocate space for greenkfsp to be reused for all fk-spectra of different radii
!
	nsp = sum(this%dsvmask_sph)*this%nsrcsph+sum(this%dsvmask_tor)*this%nsrctor
	print *,'Total number of fk-spectra in file: ',nsp
	allocate(this%greenkfsp(this%nkfmax,nsp))
	end function readHeaderGreenFrequencyWavenumber
!-----------------------------------------------------------------------
!  read frequency wavenumber spectra written by parallel gfdsvrkf
!
	function readHeaderParallelGreenFrequencyWavenumber(this,fuh,basename,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	type (file_unit_handler) :: fuh
	integer :: j,ierr,if,nsp
	character (len=*) :: basename
	logical, optional :: vsflag
	type (error_message) :: errmsg
	integer, dimension(2) :: path
	type (group_stream_access), pointer :: group => null()
	type (data_stream_access), pointer :: dset => null()
	integer, dimension(:), pointer :: id
	character (len=132) :: dsvfile
	character (len=3) :: crank
	character (len=132) :: myname = 'readHeaderParallelGreenFrequencyWavenumber'
!
	call new(errmsg,myname)
	errmsg = readHeaderBasicGreenFrequencyWavenumber(this,fuh,trim(basename)//'.000',vsflag)
	if (.level.errmsg == 2) then
		call addTraceErrorMessage(errmsg,'readHeaderParallelGreenFrequencyWavenumber')
		return
	endif
!
!  read numtasks from integer header
!
	path = (/ 5, 0 /)
	call traversePathStreamAccess(this%root,0,path,group,dset)
	call readDatasetVectorStreamAccess(dset,this%fda,id)
	this%numtasks = id(10)
	deallocate(id)
	print *,'numtasks = ',this%numtasks
!
!  open all other files and read group tree
!
	allocate(this%fda_par(this%numtasks),this%root_par(this%numtasks))
	this%fda_par(1) = this%fda
	this%root_par(1) = this%root
	do j = 2,this%numtasks
		write(crank,'(i3.3)') j-1
		dsvfile = trim(basename)//'.'//crank
		ierr = openFileStreamAccess(this%fda_par(j),get(fuh),trim(dsvfile))
		if (ierr /= 0) then
			call new(errmsg,2,trim(dsvfile)//' can not be opened',myname)
			call undo(fuh)
			return
		endif
		call readGroupStreamAccess(this%root_par(j),this%fda_par(j))
	enddo
!
!  collect nwn(if) - array, read in here because not set in readBody-routine
!
	this%nwn = 0
	path = (/ 6,0 /)
	do j = 1,this%numtasks
		call traversePathStreamAccess(this%root_par(j),0,path,group,dset)
		call readDatasetVectorStreamAccess(dset,this%fda_par(j),id)
		do if = this%nf1+j-1,this%nf2,this%numtasks
			this%nwn(if) = id(if)
		enddo
		deallocate(id)
	enddo
!
!  collect nkfmax
!
	path = (/ 5, 0 /)
	this%nkfmax = 0
	do j = 1,this%numtasks
		call traversePathStreamAccess(this%root_par(j),0,path,group,dset)
		call readDatasetVectorStreamAccess(dset,this%fda_par(j),id)
		this%nkfmax = this%nkfmax+id(4)
		deallocate(id)
	enddo
!
!  allocate space for greenkfsp to be reused for all fk-spectra of different radii
!
	nsp = sum(this%dsvmask_sph)*this%nsrcsph+sum(this%dsvmask_tor)*this%nsrctor
	print *,'Total number of fk-spectra in file: ',nsp
	allocate(this%greenkfsp(this%nkfmax,nsp))
	end function readHeaderParallelGreenFrequencyWavenumber
!-----------------------------------------------------------------------
!  open frequency-wavenumber spectrum file and read out header
!  data, fill frequency-wavenumber objects
!
	function readHeaderBasicGreenFrequencyWavenumber(this,fuh,basename,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	type (file_unit_handler) :: fuh
	character (len=*) :: basename
	logical, optional :: vsflag
	type (error_message) :: errmsg
	character (len=132) :: myname = 'readHeaderBasicGreenFrequencyWavenumber'
	integer :: ierr,lu
	integer, dimension(2) :: path
	type (group_stream_access), pointer :: group => null()
	type (data_stream_access), pointer :: dset => null()
	real, dimension(:), pointer :: d
	integer, dimension(:), pointer :: id
	double precision, dimension(:,:), pointer :: rp
	double complex, dimension(:,:), pointer :: zp
	type (flexible), dimension(:), pointer :: ft
!
	call new(errmsg,'myname')
	lu = get(fuh)
	ierr = openFileStreamAccess(this%fda,lu,trim(basename))
	if (ierr /= 0) then
		call new(errmsg,2,trim(basename)//' can not be opened',myname)
		call add(fuh,lu)
		return
	endif
	print *,'FrequencyWavenumber file ',trim(basename),' opened'
!
!  read in group tree
!
	if (present(vsflag)) then; verboseStreamAccess = vsflag; else; verboseStreamAccess = .false.; endif
	call readGroupStreamAccess(this%root,this%fda)
	print *,'group tree read in'
!
!  read out header information
!
	call traversePathStreamAccess(this%root,0,(/ 1, 0 /),group,dset)
	call readDatasetVectorStreamAccess(dset,this%fda,d)
	this%ze = d(1); this%sigma = d(2); this%df = d(3); this%dwn = d(4)
	print *,'ze sigma tlen xlen: ',this%ze,this%sigma,1./this%df,2.*mc_pi/this%dwn
	deallocate(d)
!
!  read out gemini earth model
!
	call traversePathStreamAccess(this%root,0,(/ 2, 0 /),group,dset)
	call readDatasetVectorStreamAccess(dset,this%fda,ft)
	call traversePathStreamAccess(this%root,0,(/ 3, 0 /),group,dset)
	call readDatasetDouble2DArrayStreamAccess(dset,this%fda,rp)
	call traversePathStreamAccess(this%root,0,(/ 4, 0 /),group,dset)
	call readDatasetDoubleComplex2DArrayStreamAccess(dset,this%fda,zp)
	call createFromPackedGeminiEarthModel(this%gem,ft,rp,zp)
	deallocate(ft,rp,zp)
	this%rearth = .rearth.(this%gem)
!
!  read out integer header, except numtasks
!
	path = (/ 5, 0 /)
	call traversePathStreamAccess(this%root,0,path,group,dset)
	call readDatasetVectorStreamAccess(dset,this%fda,id)
	this%istyp = id(1); this%nf1 = id(2); this%nf2 = id(3)
	this%nkfmax = id(4); this%nwnmax = id(5)
	this%nsrcsph = id(6); this%nsrctor = id(7); this%dsvstep = id(8)
	this%derivflag = id(9)
	print *,'int header sph: istyp,nf1,nf2,nkfmax,nwnmax,nsrcsph,nsrctor,dsvstep,derivflag:'
	print *, id(1:9)
	deallocate(id)
	print *,'Frequency information: fmin = ',(this%nf1-1)*this%df,', fmax = ',(this%nf2-1)*this%df
!
!  fill dsvmasks, saying which DSV components are contained in GFK-file
!
	this%dsvmask_sph = (/ 1,0,1,0,0,0 /)
	this%dsvmask_tor = (/ 1,0,0 /)
	if (this%derivflag == 1) then
		this%dsvmask_sph(5:6) = 1
		this%dsvmask_tor(3) = 1
	endif
	if (this%dsvstep == 1) then
		this%dsvmask_sph(2:4:2) = 1
		this%dsvmask_tor(2) = 1
	endif
!
!  allocate space for nwn(if)
!
	allocate(this%nwn(this%nf2))
!
	end function readHeaderBasicGreenFrequencyWavenumber
!------------------------------------------------------------------------
!  deallocate
!
	subroutine deallocGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber) :: this
	integer :: j
	call dealloc(this%gem)
	if (associated(this%nwn)) deallocate(this%nwn)
	if (associated(this%greenkfsp)) deallocate(this%greenkfsp)
	call clearGroupTree(this%root)
	call dealloc(this%fda)
	do j = 2,this%numtasks
		call dealloc(this%fda_par(j))
		call clearGroupTree(this%root_par(j))
	enddo
	end subroutine deallocGreenFrequencyWavenumber
!-----------------------------------------------------------------
!  read out all frequency wavenumber spectra written by sequential GEMINI run
!  for fixed radial node index jr. 
!  spheroidal: (DSV = 1,2,3,4,5,6 and JUMP = 1,2,3,4)
!  toroidal:   (DSV = 1,2,3   and JUMP = 1,2)
!  for given radial node index
!
	function readBodyJrGreenFrequencyWavenumber(this,jr,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	integer :: jr
	character (len=132) :: myname = 'readBodyJrGreenFrequencyWavenumber'
	integer :: ic,icp,iwe,if,nwn,isp,jsp,nsrc,cnt,offset,ndsvsph
	logical, optional :: vsflag
	character (len=3) :: motion
	integer, dimension(3) :: path
	complex, dimension(:), pointer :: d
	type (group_stream_access), pointer :: group
	type (data_stream_access), pointer :: dset
	type (error_message) :: errmsg
!
	call new(errmsg,'myname')
	if (jr <= 0 .or. jr > .nnod.(this%gem)) then
		call new(errmsg,2,'Source depth node index out of range',myname)
		return
	endif
	this%rcur = (this%gem).rnod.jr
!
	if (present(vsflag)) then; verboseStreamAccess = vsflag; else; verboseStreamAccess = .false.; endif
!
!  zero this%ksp
!
	this%ksp = 0
!
	cnt = 0                                                     ! 1D count of existing spectra
	ndsvsph = sum(this%dsvmask_sph)
	do isp = 1,9                                                ! loop over first fk-spectra index (U,R,V,S,W,T,UP,VP,WP)
		motion = getMotionFrequencyWavenumberMapping(isp)
		ic = getDsvFrequencyWavenumberMapping(isp)
		if (motion == 'sph') then
			if (this%dsvmask_sph(ic) == 0) cycle
			icp = sum(this%dsvmask_sph(1:ic))               ! position in spheroidal part
			nsrc = this%nsrcsph
			offset = 0
		else if (motion == 'tor') then
			if (this%dsvmask_tor(ic) == 0) cycle
			icp = sum(this%dsvmask_tor(1:ic))              ! position in toroidal part
			nsrc = this%nsrctor
			offset = ndsvsph*this%nsrcsph
		else
			print *,trim(myname),': invalid wave motion type'; stop
		endif
		do jsp = 1,nsrc
			cnt = cnt+1
			this%ksp(isp,jsp) = cnt                  ! mapping from isp,jsp to 1D count
			iwe = 0
			do if = this%nf1,this%nf2
				path = (/ if-this%nf1+1, jr, offset+jsp+(icp-1)*nsrc /)
				call traversePathStreamAccess(this%root,2,path,group,dset)
				call readDatasetVectorStreamAccess(dset,this%fda,d)
				nwn = size(d)
				this%greenkfsp(iwe+1:iwe+nwn,cnt) = d(1:nwn)
				deallocate(d)
				iwe = iwe+nwn
				call setNwnSelectedGreenFrequencyWavenumber(this,if,nwn)
			enddo
		enddo
	enddo
	end function readBodyJrGreenFrequencyWavenumber
!-----------------------------------------------------------------
!  read out all frequency wavenumber spectra written by parallel GEMINI run
!  for fixed radial node index jr. 
!  spheroidal: (DSV = 1,2,3,4,5,6 and JUMP = 1,2,3,4)
!  toroidal:   (DSV = 1,2,3   and JUMP = 1,2)
!  for given radial node index
!
	function readParallelBodyJrGreenFrequencyWavenumber(this,jr,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	integer :: jr
	character (len=132) :: myname = 'readParallelBodyJrGreenFrequencyWavenumber'
	integer :: ic,icp,iwe,if,nwn,isp,jsp,nsrc,cnt,j,ndsvsph,offset
	logical, optional :: vsflag
	character (len=3) :: motion
	integer, dimension(3) :: path
	complex, dimension(:), pointer :: d
	type (group_stream_access), pointer :: group
	type (data_stream_access), pointer :: dset
	type (error_message) :: errmsg
!
	call new(errmsg,'myname')
	if (jr <= 0 .or. jr > .nnod.(this%gem)) then
		call new(errmsg,2,'Source depth node index out of range',myname)
		return
	endif
	this%rcur = (this%gem).rnod.jr
!
	if (present(vsflag)) then; verboseStreamAccess = vsflag; else; verboseStreamAccess = .false.; endif
!
!  zero this%ksp
!
	this%ksp = 0
!
	cnt = 0                                                     ! 1D count of existing spectra
	ndsvsph = sum(this%dsvmask_sph)
	do isp = 1,9                                                ! loop over first fk-spectra index
		motion = getMotionFrequencyWavenumberMapping(isp)
		ic = getDsvFrequencyWavenumberMapping(isp)
		if (motion == 'sph') then
			if (this%dsvmask_sph(ic) == 0) cycle
			icp = sum(this%dsvmask_sph(1:ic))        ! position in file
			nsrc = this%nsrcsph
			offset = 0
		else if (motion == 'tor') then
			if (this%dsvmask_tor(ic) == 0) cycle
			icp = sum(this%dsvmask_tor(1:ic))        ! position in file
			nsrc = this%nsrctor
			offset = ndsvsph*this%nsrcsph
		else
			print *,'<readBodyGreenFrequencyWavenumber>: invalid wave motion type'; stop
		endif
		do jsp = 1,nsrc
			cnt = cnt+1
			this%ksp(isp,jsp) = cnt                  ! mapping from isp,jsp to 1D count
			do j = 1,this%numtasks
				do if = this%nf1+j-1,this%nf2,this%numtasks
					path = (/ (if-this%nf1)/this%numtasks+1, jr, offset+jsp+(icp-1)*nsrc /)
					call traversePathStreamAccess(this%root_par(j),2,path,group,dset)
					call readDatasetVectorStreamAccess(dset,this%fda_par(j),d)
					nwn = size(d)
					iwe = sum(this%nwn(this%nf1:if-1))
					this%greenkfsp(iwe+1:iwe+nwn,cnt) = d(1:nwn)
					deallocate(d)
				enddo
			enddo
		enddo
	enddo
	end function readParallelBodyJrGreenFrequencyWavenumber
!-------------------------------------------------------------------
!  read out all different frequency wavenumber spectra, sequential fk-spectra
!  for fixed radial node index jr. 
!  spheroidal: (DSV = 1,2,3,4,5,6 and JUMP = 1,2,3,4)
!  toroidal:   (DSV = 1,2,3   and JUMP = 1,2)
!  zs:		source depth in meters
!
! for given source depth in
!
	function readBodyZsGreenFrequencyWavenumber(this,zs,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	real :: zs
	logical, optional :: vsflag
	type (error_message) :: errmsg
	integer :: jr
	character (len=132) :: myname = 'readBodyZsGreenFrequencyWavenumber'
!
	call new(errmsg,'myname')
	jr = radialNodeIndexGreenFrequencyWavenumber(this,zs)
	if (jr <= 0 .or. jr  > .nnod.(this%gem)) then
		call new(errmsg,2,'Source depth outside range of depth nodes',myname)
		return
	endif
	print *,'Depth: ',zs,' Node index: ',jr,' True radius: ',(this%gem).rnod.jr
	errmsg = readBodyJrGreenFrequencyWavenumber(this,jr,vsflag)
	if (.level.errmsg /= 0) then
		call addTraceErrorMessage(errmsg,myname)
		return
	endif
	end function readBodyZsGreenFrequencyWavenumber
!-------------------------------------------------------------------
!  read out all different frequency wavenumber spectra, parallel fk-spectra
!  for fixed radial node index jr. 
!  spheroidal: (DSV = 1,2,3,4,5,6 and JUMP = 1,2,3,4)
!  toroidal:   (DSV = 1,2,3   and JUMP = 1,2)
!  zs:		source depth in meters
!
! for given source depth in
!
	function readParallelBodyZsGreenFrequencyWavenumber(this,zs,vsflag) result(errmsg)
	type (green_frequency_wavenumber) :: this
	real :: zs
	logical, optional :: vsflag
	type (error_message) :: errmsg
	integer :: jr
	character (len=132) :: myname = 'readParallelBodyZsGreenFrequencyWavenumber'
!
	call new(errmsg,'myname')
	jr = radialNodeIndexGreenFrequencyWavenumber(this,zs)
	if (jr <= 0 .or. jr  > .nnod.(this%gem)) then
		call new(errmsg,2,'Source depth outside range of depth nodes',myname)
		return
	endif
	print *,'Depth: ',zs,' Node index: ',jr,' True radius: ',(this%gem).rnod.jr
	errmsg = readParallelBodyJrGreenFrequencyWavenumber(this,jr,vsflag)
	if (.level.errmsg /= 0) then
		call addTraceErrorMessage(errmsg,myname)
		return
	endif
	end function readParallelBodyZsGreenFrequencyWavenumber
!--------------------------------------------------------------------
!  set nwn for frequency if
!
	subroutine setNwnSelectedGreenFrequencyWavenumber(this,if,nwn)
	type (green_frequency_wavenumber) :: this
	integer :: if,nwn
!
	this%nwn(if) = nwn
	end subroutine setNwnSelectedGreenFrequencyWavenumber
!-------------------------------------------------------------
!  Compute a cos**2 taper for given frequency
!  internal use
!
	subroutine taperGreenFrequencyWavenumber(this,if,tapfrac,greentap)
	type (green_frequency_wavenumber) :: this
	real, dimension(:) :: greentap
	real :: tapfrac,c1,wn,c2
	integer :: if,iwn1,iwn
!
	iwn1=min(nint(this%nwn(if)*(1.-tapfrac)),this%nwn(if))
	iwn1=max(1,iwn1)
	c1=(this%nwn(if)-1)*this%dwn*tapfrac*2./mc_pi
	c2=(1.-tapfrac)/tapfrac*0.5*mc_pi
	do iwn=1,iwn1-1
		greentap(iwn)=1.
	enddo
	wn=(iwn1-1)*this%dwn
	do iwn=iwn1,this%nwn(if)
		greentap(iwn)=cos( wn/c1-c2 )**2
		wn=wn+this%dwn
	enddo
	end subroutine taperGreenFrequencyWavenumber
!-----------------------------------------------------------
!  Compute wavenumber integrals for single force source
!
!  kint= 1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
!  kint= 2:       R^2/(2*pi)*int ( dk k (-2/kR) U^2 J_1 )  
!  kint= 3:       R^2/(2*pi)*int ( dk k (-kR) V^1 J_1 )  
!  kint= 4:       R^2/(2*pi)*int ( dk k (-2 V^2)(J_0-J_1/(kx)) )  
!  kint= 5:       R^2/(2*pi)*int ( dk k (-2 W^1) J_1/(kx) )  
!  kint= 6:       R^2/(2*pi)*int ( dk k (+2 W^1)(J_0-J_1/(kx)) )  
!  kint= 7:       R^2/(2*pi)*int ( dk k (+2 V^2) J_1/(kx) )
!
!  kint= 8:       R^2/(2*pi)*int ( dk k U^1_r J_0 )  
!  kint= 9:       R^2/(2*pi)*int ( dk k (-2/kR) U^2_r J_1 )  
!  kint=10:       R^2/(2*pi)*int ( dk k (-kR) V^1_r J_1 )  
!  kint=11:       R^2/(2*pi)*int ( dk k (-2 V^2_r)(J_0-J_1/(kx)) )  
!  kint=12:       R^2/(2*pi)*int ( dk k (-2 W^1_r) J_1/(kx) )  
!  kint=13:       R^2/(2*pi)*int ( dk k (+2 W^1_r)(J_0-J_1/(kx)) )  
!  kint=14:       R^2/(2*pi)*int ( dk k (+2 V^2_r) J_1/(kx) )
!
!  kint=15:       R^2/(2*pi)*int ( dk k U^1 (-kR) J_1 )  
!  kint=16:       R^2/(2*pi)*int ( dk k -2 U^2 (J_0-J_1/(kx)) )  
!  kint=17:       R^2/(2*pi)*int ( dk k -(kR)^2 V^1 (J_0-J_1/(kx)) )  
!  kint=18:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR)(-J_1+J_2/(kx)) )  
!  kint=19:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR) J_2/(kx) )  
!  kint=20:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR)(-J_1+J_2/(kx)) )  
!  kint=21:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR) J_2(kx)/(kx) )
!
!  kint=22:       R^2/(2*pi)*int ( dk k R^1 J_0 )  
!  kint=23:       R^2/(2*pi)*int ( dk k (-2/kR) R^2 J_1 )  
!
!  for one receiver at x and ALL frequencies. This saves time
!  because the Bessel functions do not depend on frequency
!  and the wavenumber sampling is independent of frequency.
!  A complex array zdis(f) is returned containing the
!  value of the integral.
!  tapfra! is the fraction of wavenumbers multiplied by a cosine taper.
!
!  Note: Unit of zdis(if) is returned in nanometers/N
!
	subroutine wnintForceGreenFrequencyWavenumber(this,kint,xx,tapfrac,besselj,zdis)
	type (green_frequency_wavenumber) :: this
	complex, dimension(:) :: zdis
	real, dimension(:,0:) :: besselj
	real, dimension(:), allocatable :: wnbesj,greentap
	real :: x,xx,wn3,wn,wn3besy,xdk,wnbesy,tapfrac,by0,by1,by2
	integer :: kint,j3,iwn,if,j33,j,nwn,nwe,jsp,is,js
	complex :: zsum
	type (error_message) :: errmsg
	character (len=132) :: myname = 'wnintForceGreenFrequencyWavenumber'
	real :: bessy0,bessy1
	external :: bessy0,bessy1
!
	call getIspJspWavenumberIntegralMapping(kint,this%istyp,is,js)
	jsp = this%ksp(is,js)
	if (jsp == 0) then
		call new(errmsg,2,'No fk-spectra calculated for desired wavenumber integral',myname)
		call print(errmsg)
		stop
	endif
!
	if(xx < 3./(this%dwn*1.e9)) then
		x = 3./(this%dwn*1.e9)
		print *,'<wnintForceGreenFrequencyWavenumber>: WARNING: distance set to ',x,' instead less'
	else
		x = xx
	endif
	wn3=3./x
!
!  find index of first wn-point greater than 3
!
	j3=int(wn3/this%dwn)+2
	wn3=(j3-1)*this%dwn
!
!  Bessel terms over whole interval
!  division by zero is unproblematic because we avoid wn=0
!  Bessel terms at k=0 vanish anyway
!  and x shouldn't be zero either
!
	allocate(wnbesj(this%nwnmax))
	wn=this%dwn
	do iwn=2,this%nwnmax
		if(kint.eq.1.or.kint.eq.8.or.kint.eq.22)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.2.or.kint.eq.9.or.kint.eq.23)  wnbesj(iwn)=-2*besselj(iwn,1)/this%rearth
		if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wnbesj(iwn)=-wn*wn*besselj(iwn,1)*this%rearth
		if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wnbesj(iwn)=-2.*wn*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.5.or.kint.eq.12) wnbesj(iwn)=-2.*wn*besselj(iwn,1)/(wn*x)
		if(kint.eq.6.or.kint.eq.13) wnbesj(iwn)=+2.*wn*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.7.or.kint.eq.14) wnbesj(iwn)=+2.*wn*besselj(iwn,1)/(wn*x)
		if(kint.eq.17) wnbesj(iwn)=-wn*(wn*this%rearth)**2*( besselj(iwn,0)-besselj(iwn,1)/(wn*x) )
		if(kint.eq.18) wnbesj(iwn)=-2*wn*wn*this%rearth*( -besselj(iwn,1)+besselj(iwn,2)/(wn*x) )
		if(kint.eq.19) wnbesj(iwn)=+2*wn*wn*this%rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.20) wnbesj(iwn)=+2*wn*wn*this%rearth*( -besselj(iwn,1)+besselj(iwn,2)/(wn*x) )
		if(kint.eq.21) wnbesj(iwn)=-2*wn*wn*this%rearth*besselj(iwn,2)/(wn*x)
		wn=wn+this%dwn
	enddo
!
!  Neumann terms at wn3
!
	by0=bessy0(wn3*x)
	by1=bessy1(wn3*x)
	by2=-by0+2./(wn3*x)*by1
	wn3besy = 0.0
	if(kint.eq.1.or.kint.eq.8.or.kint.eq.22) wn3besy=wn3*by0
	if(kint.eq.2.or.kint.eq.9.or.kint.eq.23) wn3besy=-2*by1/this%rearth
	if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wn3besy=-wn3*wn3*by1*this%rearth
	if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wn3besy=-2.*wn3*(by0-by1/(wn3*x))
	if(kint.eq.5.or.kint.eq.12) wn3besy=-2.*wn3*by1/(wn3*x)
	if(kint.eq.6.or.kint.eq.13) wn3besy=+2.*wn3*(by0-by1/(wn3*x))
	if(kint.eq.7.or.kint.eq.14) wn3besy=+2.*wn3*by1/(wn3*x)
	if(kint.eq.17) wn3besy=-wn3*(wn3*this%rearth)**2*(by0-by1/(wn3*x))
	if(kint.eq.18) wn3besy=-2*wn3*wn3*this%rearth*(-by1+by2/(wn3*x))
	if(kint.eq.19) wn3besy=+2*wn3*wn3*this%rearth*by2/(wn3*x)
	if(kint.eq.20) wn3besy=+2*wn3*wn3*this%rearth*(-by1+by2/(wn3*x))
	if(kint.eq.21) wn3besy=-2*wn3*wn3*this%rearth*by2/(wn3*x)
!
!  perform integration for each frequency
!  Use trapezoidal rule if wn*x < 3 and Filon else
!
	allocate(greentap(this%nwnmax))
	nwe=0
	do if=this%nf1,this%nf2
		nwn=this%nwn(if)
		call taperGreenFrequencyWavenumber(this,if,tapfrac,greentap)
		j33=min(j3,nwn)
		zdis(if)=0.
		do j=2,j33-1
			zdis(if)=zdis(if)+this%greenkfsp(nwe+j,jsp)*wnbesj(j)*greentap(j)
		enddo
		zdis(if)=zdis(if)+0.5*this%greenkfsp(nwe+j33,jsp)*wnbesj(j33)*greentap(j33)
		zdis(if)=zdis(if)*this%dwn
		if(j33.eq.nwn) goto 11
!
!  use the Filon integration for kx > 3
!  Neumann terms at the end of integration interval (k_N)
!
		wn=(this%nwn(if)-1)*this%dwn
		by0=bessy0(wn*x)
		by1=bessy1(wn*x)
		by2=-by0+2./(wn*x)*by1
		wnbesy = 0.0
		if(kint.eq.1.or.kint.eq.8.or.kint.eq.22) wnbesy=wn*by0
		if(kint.eq.2.or.kint.eq.9.or.kint.eq.23) wnbesy=-2*by1/this%rearth
		if(kint.eq.3.or.kint.eq.10.or.kint.eq.15) wnbesy=-wn*wn*by1*this%rearth
		if(kint.eq.4.or.kint.eq.11.or.kint.eq.16) wnbesy=-2.*wn*(by0-by1/(wn*x))
		if(kint.eq.5.or.kint.eq.12) wnbesy=-2.*wn*by1/(wn*x)
		if(kint.eq.6.or.kint.eq.13) wnbesy=+2.*wn*(by0-by1/(wn*x))
		if(kint.eq.7.or.kint.eq.14) wnbesy=+2.*wn*by1/(wn*x)
		if(kint.eq.17) wnbesy=-wn*(wn*this%rearth)**2*(by0-by1/(wn*x))
		if(kint.eq.18) wnbesy=-2*wn*wn*this%rearth*(-by1+by2/(wn*x))
		if(kint.eq.19) wnbesy=+2*wn*wn*this%rearth*by2/(wn*x)
		if(kint.eq.20) wnbesy=+2*wn*wn*this%rearth*(-by1+by2/(wn*x))
		if(kint.eq.21) wnbesy=-2*wn*wn*this%rearth*by2/(wn*x)
!
		xdk=x*this%dwn
		zsum=this%greenkfsp(nwe+j3,jsp)*wnbesj(j3)*greentap(j3)
		do j=j3+1,nwn-1
			zsum=zsum+2.*this%greenkfsp(nwe+j,jsp)*wnbesj(j)*greentap(j)
		enddo
		zsum=zsum+this%greenkfsp(nwe+nwn,jsp)*wnbesj(nwn)*greentap(nwn)
		zdis(if)=zdis(if)+zsum/x*(1-cos(xdk))/xdk
		zdis(if)=zdis(if)+(this%greenkfsp(nwe+nwn,jsp)*wnbesy*greentap(nwn) &
		               &   -this%greenkfsp(nwe+j3,jsp)*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
!
!  multiply by R^2/(2*pi)
!
 11		zdis(if)=zdis(if)*this%rearth**2/(2.*mc_pi)
!
!  convert to nanometers per N or milli-Pa per N
!
		zdis(if)=zdis(if)*1.e-3
!
!  update address of last wavenumber spectrum value
!	
		nwe=nwe+nwn
	enddo
	deallocate(greentap,wnbesj)
	end subroutine wnintForceGreenFrequencyWavenumber
!-----------------------------------------------------------
!  Compute wavenumber integrals for moment tensor source
!
!  kint=1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
!  kint=2:       R^2/(2*pi)*int ( dk k U^2 J_0 )
!  kint=3:       R^2/(2*pi)*int ( dk k U^3 2/(kR) J_1 )
!  kint=4:       R^2/(2*pi)*int ( dk k U^4 2J_2 )
!  kint=5:       R^2/(2*pi)*int ( dk k V^1 -kR J_1 )
!  kint=6:       R^2/(2*pi)*int ( dk k V^2 -kR J_1 )
!  kint=7:       R^2/(2*pi)*int ( dk k V^3 2 (J_0-J_1/(kx)) )
!  kint=8:       R^2/(2*pi)*int ( dk k V^4 2kR (J_1-2J_2/(kx)) )
!  kint=9:       R^2/(2*pi)*int ( dk k W^1 2J_1/(kx) )
!  kint=10:      R^2/(2*pi)*int ( dk k W^2 4kRJ_2/(kx) )
!  kint=11:      R^2/(2*pi)*int ( dk k W^1 2(J_0-J_1/(kx)) )
!  kint=12:      R^2/(2*pi)*int ( dk k W^2 -2kR(J_1-2J_2/(kx)) )
!  kint=13:      R^2/(2*pi)*int ( dk k V^3 2J_1/(kx) )
!  kint=14:      R^2/(2*pi)*int ( dk k V^4 -4kRJ_2/(kx) )
!
!  kint=15:      R^2/(2*pi)*int ( dk k U^1_r J_0 )
!  kint=16:      R^2/(2*pi)*int ( dk k U^2_r J_0 )
!  kint=17:      R^2/(2*pi)*int ( dk k U^3_r 2/(kR) J_1 )
!  kint=18:      R^2/(2*pi)*int ( dk k U^4_r 2J_2 )
!  kint=19:      R^2/(2*pi)*int ( dk k V^1_r -kR J_1 )
!  kint=20:      R^2/(2*pi)*int ( dk k V^2_r -kR J_1 )
!  kint=21:      R^2/(2*pi)*int ( dk k V^3_r 2(J_0-J_1/(kx)) )
!  kint=22:      R^2/(2*pi)*int ( dk k V^4_r 2kR(J_1-2J_2/(kx)) )
!  kint=23:      R^2/(2*pi)*int ( dk k W^1_r 2J_1/(kx) )
!  kint=24:      R^2/(2*pi)*int ( dk k W^2_r 4kRJ_2/(kx) )
!  kint=25:      R^2/(2*pi)*int ( dk k W^1_r 2(J_0-J_1/(kx)) )
!  kint=26:      R^2/(2*pi)*int ( dk k W^2_r -2kR(J_1-2J_2/(kx)) )
!  kint=27:      R^2/(2*pi)*int ( dk k V^3_r 2J_1/(kx) )
!  kint=28:      R^2/(2*pi)*int ( dk k V^4_r -4kRJ_2/(kx) )
!
!  kint=29:      R^2/(2*pi)*int ( dk k U^1 -kR J_1 )  
!  kint=30:      R^2/(2*pi)*int ( dk k U^2 -kR J_1 )
!  kint=31:      R^2/(2*pi)*int ( dk k U^3 2 (J_0-J_1/(kx)) )
!  kint=32:      R^2/(2*pi)*int ( dk k U^4 2kR (J_1-2J_2/(kx)) )
!  kint=33:      R^2/(2*pi)*int ( dk k V^1 -(kR)^2 (J_0-J_1/(kx)) )
!  kint=34:      R^2/(2*pi)*int ( dk k V^2 -(kR)^2 (J_0-J_1/(kx)) )
!  kint=35:      R^2/(2*pi)*int ( dk k V^3 -2kR (J_1-J_2/(kx)) )
!  kint=36:      R^2/(2*pi)*int ( dk k V^4 2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
!  kint=37:      R^2/(2*pi)*int ( dk k W^1 -2kR J_2/(kx) )
!  kint=38:      R^2/(2*pi)*int ( dk k W^2 4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
!  kint=39:      R^2/(2*pi)*int ( dk k W^1 -2kR (J_1-J_2/(kx)) )
!  kint=40:      R^2/(2*pi)*int ( dk k W^2 -2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
!  kint=41:      R^2/(2*pi)*int ( dk k V^3 -2kR J_2/(kx) )
!  kint=42:      R^2/(2*pi)*int ( dk k V^4 -4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
!
!  kint=43:       R^2/(2*pi)*int ( dk k R^1 J_0 )  
!  kint=44:       R^2/(2*pi)*int ( dk k R^2 J_0 )
!  kint=45:       R^2/(2*pi)*int ( dk k R^3 2/(kR) J_1 )
!  kint=46:       R^2/(2*pi)*int ( dk k R^4 2J_2 )
!
!  for one receiver at x and ALL frequencies. This saves time
!  because the Bessel functions do not depend on frequency
!  and the wavenumber sampling is independent of frequency.
!  A complex array zdis(f) is returned containing the
!  value of the integral.
!  tapfra! is the fraction of wavenumbers multiplied by a cosine taper.
!
!  Note: Unit of zdis(if) is returned in nanometers/(Nm)
!
	subroutine wnintMomentGreenFrequencyWavenumber(this,kint,xx,tapfrac,besselj,zdis)
	type (green_frequency_wavenumber) :: this
	complex, dimension(:) :: zdis
	real, dimension(:,0:) :: besselj
	real, dimension(:), allocatable :: greentap,wnbesj
	real :: xx,x,wn3,wn,wn3besy,xdk,wnbesy,tapfrac,by0,by1,by2
	integer :: kint,j3,iwn,if,j33,j,nwn,nwe,jsp,is,js
	complex :: zsum
	type (error_message) :: errmsg
	character (len=132) :: myname = 'wnintMomentGreenFrequencyWavenumber'
	real :: bessy0,bessy1
	external :: bessy0,bessy1
!
	call getIspJspWavenumberIntegralMapping(kint,this%istyp,is,js)
	jsp = this%ksp(is,js)
	if (jsp == 0) then
		call new(errmsg,2,'No fk-spectra calculated for desired wavenumber integral',myname)
		call print(errmsg)
		stop
	endif
!
!  if x is too small j3 may become larger than representable by integer*4
!  therefore choose x such that j3 = 1e9
!
	if(xx < 3./(this%dwn*1.e9)) then
		x = 3./(this%dwn*1.e9)
		print *,'<wnintMomentGreenFrequencyWavenumber>: WARNING: distance set to ',x,' instead less'
	else
		x = xx
	endif
	wn3=3./x
!
!  find index of first wn-point greater than 3
!
	j3=int(wn3/this%dwn)+2
	wn3=(j3-1)*this%dwn
!
!  Bessel terms over whole interval
!  division by zero is unproblemati! because we avoid wn=0
!  Bessel terms at k=0 vanish anyway
!  and x shouldn't be zero either
!
	allocate(wnbesj(this%nwnmax))
	wn=this%dwn
	do iwn=2,this%nwnmax
		if(kint.eq.1.or.kint.eq.15.or.kint.eq.43)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.2.or.kint.eq.16.or.kint.eq.44)  wnbesj(iwn)=wn*besselj(iwn,0)
		if(kint.eq.3.or.kint.eq.17.or.kint.eq.45)  wnbesj(iwn)=2./this%rearth*besselj(iwn,1)
		if(kint.eq.4.or.kint.eq.18.or.kint.eq.46)  wnbesj(iwn)=wn*2.*besselj(iwn,2)
		if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wnbesj(iwn)=wn*(-wn*this%rearth)*besselj(iwn,1)
		if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wnbesj(iwn)=wn*(-wn*this%rearth)*besselj(iwn,1)
		if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wnbesj(iwn)=wn*2.*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wnbesj(iwn)=wn*2.*wn*this%rearth*(besselj(iwn,1)-2.*besselj(iwn,2)/(wn*x))
		if(kint.eq.9.or.kint.eq.23)  wnbesj(iwn)=wn*2.*besselj(iwn,1)/(wn*x)
		if(kint.eq.10.or.kint.eq.24) wnbesj(iwn)=wn*4.*wn*this%rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.11.or.kint.eq.25) wnbesj(iwn)=wn*2.*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.12.or.kint.eq.26) wnbesj(iwn)=-wn*2.*wn*this%rearth*(besselj(iwn,1)-2.*besselj(iwn,2)/(wn*x))
		if(kint.eq.13.or.kint.eq.27) wnbesj(iwn)=wn*2.*besselj(iwn,1)/(wn*x)
		if(kint.eq.14.or.kint.eq.28) wnbesj(iwn)=-wn*4.*wn*this%rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.33.or.kint.eq.34) wnbesj(iwn)=-wn*(wn*this%rearth)**2*(besselj(iwn,0)-besselj(iwn,1)/(wn*x))
		if(kint.eq.35.or.kint.eq.39) wnbesj(iwn)=-wn*2*wn*this%rearth*(besselj(iwn,1)-besselj(iwn,2)/(wn*x))
		if(kint.eq.36) wnbesj(iwn)=+wn*2*(wn*this%rearth)**2*((6./(wn*x)**2-1.)*besselj(iwn,2)-besselj(iwn,1)/(wn*x))
		if(kint.eq.37.or.kint.eq.41) wnbesj(iwn)=-wn*2*wn*this%rearth*besselj(iwn,2)/(wn*x)
		if(kint.eq.38) wnbesj(iwn)=+wn*4.*(wn*this%rearth)**2*(-3.*besselj(iwn,2)/(wn*x)**2+besselj(iwn,1)/(wn*x))
		if(kint.eq.40) wnbesj(iwn)=-wn*2*(wn*this%rearth)**2*((6./(wn*x)**2-1.)*besselj(iwn,2)-besselj(iwn,1)/(wn*x))
		if(kint.eq.42) wnbesj(iwn)=-wn*4.*(wn*this%rearth)**2*(-3.*besselj(iwn,2)/(wn*x)**2+besselj(iwn,1)/(wn*x))
		wn=wn+this%dwn
	enddo
!
!  Neumann terms at wn3
!
	by0=bessy0(wn3*x)
	by1=bessy1(wn3*x)
	by2=-by0+2./(wn3*x)*by1
	wn3besy = 0.0
	if(kint.eq.1.or.kint.eq.15.or.kint.eq.43)  wn3besy=wn3*by0
	if(kint.eq.2.or.kint.eq.16.or.kint.eq.44)  wn3besy=wn3*by0
	if(kint.eq.3.or.kint.eq.17.or.kint.eq.45)  wn3besy=2./this%rearth*by1
	if(kint.eq.4.or.kint.eq.18.or.kint.eq.46)  wn3besy=wn3*2.*by2
	if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wn3besy=wn3*(-wn3*this%rearth)*by1
	if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wn3besy=wn3*(-wn3*this%rearth)*by1
	if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wn3besy=wn3*2.*(by0-by1/(wn3*x))
	if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wn3besy=wn3*2.*wn3*this%rearth*(by1-2.*by2/(wn3*x))
	if(kint.eq.9.or.kint.eq.23)  wn3besy=wn3*2.*by1/(wn3*x)
	if(kint.eq.10.or.kint.eq.24) wn3besy=wn3*4.*wn3*this%rearth*by2/(wn3*x)
	if(kint.eq.11.or.kint.eq.25) wn3besy=wn3*2.*(by0-by1/(wn3*x))
	if(kint.eq.12.or.kint.eq.26) wn3besy=-wn3*2.*wn3*this%rearth*(by1-2.*by2/(wn3*x))
	if(kint.eq.13.or.kint.eq.27) wn3besy=wn3*2.*by1/(wn3*x)
	if(kint.eq.14.or.kint.eq.28) wn3besy=-wn3*4.*wn3*this%rearth*by2/(wn3*x)
	if(kint.eq.33.or.kint.eq.34) wn3besy=-wn3*(wn3*this%rearth)**2*(by0-by1/(wn3*x))
	if(kint.eq.35.or.kint.eq.39) wn3besy=-wn3*2*wn3*this%rearth*(by1-by2/(wn3*x))
	if(kint.eq.36) wn3besy=+wn3*2*(wn3*this%rearth)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
	if(kint.eq.37.or.kint.eq.41) wn3besy=-wn3*2*wn3*this%rearth*by2/(wn3*x)
	if(kint.eq.38) wn3besy=+wn3*4.*(wn3*this%rearth)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
	if(kint.eq.40) wn3besy=-wn3*2*(wn3*this%rearth)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
	if(kint.eq.42) wn3besy=-wn3*4.*(wn3*this%rearth)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
!
!  perform integration for each frequency
!  Use trapezoidal rule if wn*x < 3 and Filon else
!
	allocate(greentap(this%nwnmax))
	nwe=0
	do if=this%nf1,this%nf2
		nwn=this%nwn(if)
		call taperGreenFrequencyWavenumber(this,if,tapfrac,greentap)
		j33=min(j3,nwn)
		zdis(if)=0.
		do j=2,j33-1
			zdis(if)=zdis(if)+this%greenkfsp(nwe+j,jsp)*wnbesj(j)*greentap(j)
		enddo
		zdis(if)=zdis(if)+0.5*this%greenkfsp(nwe+j33,jsp)*wnbesj(j33)*greentap(j33)
		zdis(if)=zdis(if)*this%dwn
		if(j33.eq.nwn) goto 11
!
!  use the Filon integration for kx > 3
!  Neumann terms at the end of integration interval (k_N)
!
		wn=(this%nwn(if)-1)*this%dwn
		by0=bessy0(wn*x)
		by1=bessy1(wn*x)
		by2=-by0+2./(wn*x)*by1
		wnbesy = 0.0
		if(kint.eq.1.or.kint.eq.15.or.kint.eq.43)  wnbesy=wn*by0
		if(kint.eq.2.or.kint.eq.16.or.kint.eq.44)  wnbesy=wn*by0
		if(kint.eq.3.or.kint.eq.17.or.kint.eq.45)  wnbesy=2./this%rearth*by1
		if(kint.eq.4.or.kint.eq.18.or.kint.eq.46)  wnbesy=wn*2.*by2
		if(kint.eq.5.or.kint.eq.19.or.kint.eq.29)  wnbesy=wn*(-wn*this%rearth)*by1
		if(kint.eq.6.or.kint.eq.20.or.kint.eq.30)  wnbesy=wn*(-wn*this%rearth)*by1
		if(kint.eq.7.or.kint.eq.21.or.kint.eq.31)  wnbesy=wn*2.*(by0-by1/(wn*x))
		if(kint.eq.8.or.kint.eq.22.or.kint.eq.32)  wnbesy=wn*2.*wn*this%rearth*(by1-2.*by2/(wn*x))
		if(kint.eq.9.or.kint.eq.23)  wnbesy=wn*2.*by1/(wn*x)
		if(kint.eq.10.or.kint.eq.24) wnbesy=wn*4.*wn*this%rearth*by2/(wn*x)
		if(kint.eq.11.or.kint.eq.25) wnbesy=wn*2.*(by0-by1/(wn*x))
		if(kint.eq.12.or.kint.eq.26) wnbesy=-wn*2.*wn*this%rearth*(by1-2.*by2/(wn*x))
		if(kint.eq.13.or.kint.eq.27) wnbesy=wn*2.*by1/(wn*x)
		if(kint.eq.14.or.kint.eq.28) wnbesy=-wn*4.*wn*this%rearth*by2/(wn*x)
		if(kint.eq.33.or.kint.eq.34) wnbesy=-wn*(wn*this%rearth)**2*(by0-by1/(wn*x))
		if(kint.eq.35.or.kint.eq.39) wnbesy=-wn*2*wn*this%rearth*(by1-by2/(wn*x))
		if(kint.eq.36) wnbesy=+wn*2*(wn*this%rearth)**2*((6./(wn*x)**2-1.)*by2-by1/(wn*x))
		if(kint.eq.37.or.kint.eq.41) wnbesy=-wn*2*wn*this%rearth*by2/(wn*x)
		if(kint.eq.38) wnbesy=+wn*4.*(wn*this%rearth)**2*(-3.*by2/(wn*x)**2+by1/(wn*x))
		if(kint.eq.40) wnbesy=-wn*2*(wn*this%rearth)**2*((6./(wn*x)**2-1.)*by2-by1/(wn*x))
		if(kint.eq.42) wnbesy=-wn*4.*(wn*this%rearth)**2*(-3.*by2/(wn*x)**2+by1/(wn*x))
		xdk=x*this%dwn
		zsum=this%greenkfsp(nwe+j3,jsp)*wnbesj(j3)*greentap(j3)
		do j=j3+1,nwn-1
			zsum=zsum+2.*this%greenkfsp(nwe+j,jsp)*wnbesj(j)*greentap(j)
		enddo
		zsum=zsum+this%greenkfsp(nwe+nwn,jsp)*wnbesj(nwn)*greentap(nwn)
		zdis(if)=zdis(if)+zsum/x*(1-cos(xdk))/xdk
		zdis(if)=zdis(if)+(this%greenkfsp(nwe+nwn,jsp)*wnbesy*greentap(nwn) &
		                 & -this%greenkfsp(nwe+j3,jsp)*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
!
!  multiply by R^2/(2*pi)
!
 11		zdis(if)=zdis(if)*this%rearth**2/(2.*mc_pi)
!
!  convert to nanometer per Nm or milli-Pa per Nm
!
		zdis(if)=zdis(if)*1.e-6
!
!  update address of last wavenumber-spectrum value
!
		nwe=nwe+nwn
	enddo
	deallocate(greentap,wnbesj)
	end subroutine wnintMomentGreenFrequencyWavenumber
!-----------------------------------------------------------
!  Compute wavenumber integrals for single force source in 2D
!  Integration from zero to kmax
!  There are two types:
!  int_0^kmax S(k) cos(kx) dk  with S(-k) = +S(k) and
!  int_0^kmax A(k) sin(kx) dk  with A(-k) = -A(k)
!  The first is equal to   1/2   *int_-kmax^kmax S(k) exp(ikx) dk
!  The second is equal to  1/(2i)*int_-kmax^kmax A(k) exp(ikx) dk
!
!  kint= 1:       R^2/(pi)*int ( dk U^1 cos(kx) )
!  kint= 2:       R^2/(pi)*int ( dk (-2/kR) U^2 sin(kx) )
!  kint= 3:       R^2/(pi)*int ( dk (-kR) V^1 sin(kx) )
!  kint= 4:       R^2/(pi)*int ( dk (-2 V^2) cos(kx) )
!  kint= 6:       R^2/(pi)*int ( dk (+2 W^1) cos(kx) )
!
! d/dr:
!
!  kint= 8:       R^2/(pi)*int ( dk U^1_r cos(kx) )
!  kint= 9:       R^2/(pi)*int ( dk (-2/kR) U^2_r sin(kx) )
!  kint=10:       R^2/(pi)*int ( dk (-kR) V^1_r sin(kx) )
!  kint=11:       R^2/(pi)*int ( dk (-2 V^2_r) cos(kx) )
!  kint=13:       R^2/(pi)*int ( dk (+2 W^1_r) cos(kx) )
!
! d/dx
!
!  kint=15:       R^2/(pi)*int ( dk (-k U^1) sin(kx) )
!  kint=16:       R^2/(pi)*int ( dk (-2/R) U^2 cos(kx) )
!  kint=17:       R^2/(pi)*int ( dk (-kR) kV^1 cos(kx) )
!  kint=18:       R^2/(pi)*int ( dk (+2 kV^2) sin(kx) )
!  kint=20:       R^2/(pi)*int ( dk (-2 kW^1) sin(kx) )
!
! Pressure
!
!  kint=22:       R^2/(pi)*int ( dk R^1 cos(kx) )  
!  kint=23:       R^2/(pi)*int ( dk (-2/kR) R^2 sin(kx) )  
!
!  for all receivers at x_j and ALL frequencies. This saves time
!  because the wavenumber sampling is independent of frequency.
!  A complex array zdis(f,x_j) is returned containing the
!  value of the integral.
!  tapfrac is the fraction of wavenumbers multiplied by a cosine taper.
!
!  Note: Unit of zdis(if) is returned in nanometers/N/dim(k)
!
!	subroutine wnint2DForceGreenFrequencyWavenumber(this,kint,x,tapfrac,zdis)
!	type (green_frequency_wavenumber) :: this
!	real, dimension(:) :: x
!	complex, dimension(:,:) :: zdis
!	real, dimension(:), allocatable :: greentap,wn
!	double precision, dimension(:), allocatable :: xi
!	double complex, dimension(:), allocatable :: fk,zy
!	real :: tapfrac
!	integer :: kint,if,j,nwn,nwe,jsp,is,js,ier,nr
!	type (error_message) :: errmsg
!	character (len=132) :: myname = 'wnint2DForceGreenFrequencyWavenumber'
!
!	call getIspJspWavenumberIntegralMapping(kint,this%istyp,is,js)
!	jsp = this%ksp(is,js)
!	if (jsp == 0) then
!		call new(errmsg,2,'No fk-spectra calculated for desired wavenumber integral',myname)
!		call print(errmsg)
!		stop
!	endif
!
!	nr = size(x)
!
!  perform integration for each frequency, use NUFFT
!
!	allocate(greentap(this%nwnmax),wn(this%nwnmax),xi(nr),fk(-this%nwnmax:this%nwnmax),zy(nr))
!	xi = this%dwn*x                                                                         ! normalized x-values
!	nwe=0
!	do if = this%nf1,this%nf2
!		nwn = this%nwn(if)
!		call taperGreenFrequencyWavenumber(this,if,tapfrac,greentap)
!		wn = (/ ((j-1)*this%dwn,j=1,nwn) /)
!		fk = 0.d0
!
!  fill spectrum for NUFFT, fk(k=0) = fk(0) = 0
!
!		if (kint == 1 .or. kint == 22) fk(1:nwn-1) = this%greenkfsp(nwe+2:nwe+nwn,jsp)*greentap(2:nwn)
!		if (kint == 2 .or. kint == 23) fk(1:nwn-1) = -2./this%rearth/wn(2:nwn)*this%greenkfsp(nwe+2:nwe+nwn,jsp)*greentap(2:nwn)
!		if (kint == 3) fk(1:nwn-1) = -wn(2:nwn)*this%rearth*this%greenkfsp(nwe+2:nwe+nwn,jsp)*greentap(2:nwn)
!		if (kint == 4) fk(1:nwn-1) = -2.*this%greenkfsp(nwe+2:nwe+nwn,jsp)*greentap(2:nwn)
!		if (kint == 6) fk(1:nwn-1) = +2.*this%greenkfsp(nwe+2:nwe+nwn,jsp)*greentap(2:nwn)
!
!  negative wavenumbers
!
!		if (kint == 1 .or. kint == 22) fk(-1:-nwn+1:-1) = +fk(1:nwn-1)
!		if (kint == 2 .or. kint == 23) fk(-1:-nwn+1:-1) = -fk(1:nwn-1)
!		if (kint == 3) fk(-1:-nwn+1:-1) = -fk(1:nwn-1)
!		if (kint == 4) fk(-1:-nwn+1:-1) = +fk(1:nwn-1)
!		if (kint == 6) fk(-1:-nwn+1:-1) = +fk(1:nwn-1)
!		fk(-nwn) = 0.
!
!  Integration
!
!		call nonUniformOutputSamplesFourierTransform(nr,xi,zy,1,1.d-6,2*nwn,fk,ier)
!
!		if (kint == 1 .or. kint == 4 .or. kint == 6 .or. kint == 22) zdis(if,1:nr) = zy(1:nr)*this%dwn*this%rearth**2/(2.*mc_pi)
!		if (kint == 2 .or. kint == 3 .or. kint == 23) zdis(if,1:nr) = -mc_ci*zy(1:nr)*this%dwn*this%rearth**2/(2.*mc_pi)
!
!  convert to nanometers per N or milli-Pa per N
!
!		zdis(if,1:nr)=zdis(if,1:nr)*1.e-3
!		if (kint == 2 .and. if == 20) then
!			print *,kint,if,zdis(if,1)
!		endif
!
!  update address of last wavenumber spectrum value
!	
!		nwe=nwe+nwn
!	enddo
!	deallocate(greentap,wn,xi,fk,zy)
!	end subroutine wnint2DForceGreenFrequencyWavenumber
!--------------------------------------------------------------------------
!  Precompute values of Bessel functions at wavenumber points
!  and receivers
!
	subroutine besselGreenFrequencyWavenumber(this,x,besselj)
	type (green_frequency_wavenumber) :: this
	real, dimension(:,0:) :: besselj
	real :: x,wn
	integer iwn
	real :: bessj0,bessj1
	external :: bessj0,bessj1
!
	wn=this%dwn
	do iwn=2,this%nwnmax
		besselj(iwn,0)=bessj0(wn*x)
		besselj(iwn,1)=bessj1(wn*x)
		besselj(iwn,2)=-besselj(iwn,0)+2./(wn*x)*besselj(iwn,1)
		wn=wn+this%dwn
	enddo
	end subroutine besselGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get wavenumber intergral index array
!
	subroutine getWavenumberIntegralIndexArrayGreenFrequencyWavenumber(this,nwint,kint)
	type (green_frequency_wavenumber) :: this
	integer :: nwint,i
	integer, dimension(:), pointer :: kint
!
!  Moment tensor source
!
	if (this%istyp == 1) then
		if (this%dsvstep == 1 .and. this%derivflag == 1) then
			nwint = 46
			allocate(kint(nwint))
			kint = (/ (i,i=1,46) /)
		else if (this%dsvstep == 1 .and. this%derivflag == 0) then
			nwint = 18
			allocate(kint(nwint))
			kint = (/ (i,i=1,14),(i,i=43,46) /)
		else if (this%dsvstep == 2 .and. this%derivflag == 1) then
			nwint = 42
			allocate(kint(nwint))
			kint = (/ (i,i=1,42) /)
		else if (this%dsvstep == 2 .and. this%derivflag == 0) then
			nwint = 14
			allocate(kint(nwint))
			kint = (/ (i,i=1,14) /)
		endif
!
!  Force source
!
	else if (this%istyp == 0) then
		if (this%dsvstep == 1 .and. this%derivflag == 1) then
			nwint = 23
			allocate(kint(nwint))
			kint = (/ (i,i=1,23) /)
		else if (this%dsvstep == 1 .and. this%derivflag == 1) then
			nwint = 9
			allocate(kint(nwint))
			kint = (/ (i,i=1,7),(i,i=22,23) /)
		else if (this%dsvstep == 2 .and. this%derivflag == 1) then
			nwint = 21
			allocate(kint(nwint))
			kint = (/ (i,i=1,21) /)
		else if (this%dsvstep == 2 .and. this%derivflag == 1) then
			nwint = 7
			allocate(kint(nwint))
			kint = (/ (i,i=1,7) /)
		endif
	endif
	end subroutine getWavenumberIntegralIndexArrayGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get a pointer to gfk-spectrum for given second index
!
	function getSelectedSpectrumGreenFrequencyWavenumber(this,jsp) result(fk)
 	type (green_frequency_wavenumber) :: this
	complex, dimension(:), pointer :: fk
	integer :: jsp
	fk => this%greenkfsp(:,jsp)
	end function getSelectedSpectrumGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get linear spectra count from DSV-index and jump index
!
	integer function getLinearCountGreenFrequencyWavenumber(this,isp,jsp)
	type (green_frequency_wavenumber), intent(in) :: this
	integer :: isp,jsp
	getLinearCountGreenFrequencyWavenumber = this%ksp(isp,jsp)
	end function getLinearCountGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get source type
!
	integer function getSourceTypeGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getSourceTypeGreenFrequencyWavenumber = this%istyp
	end function getSourceTypeGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nf1
!
	integer function getNf1GreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNf1GreenFrequencyWavenumber = this%nf1
	end function getNf1GreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nf2
!
	integer function getNf2GreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNf2GreenFrequencyWavenumber = this%nf2
	end function getNf2GreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nkfmax
!
	integer function getNkfmaxGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNkfmaxGreenFrequencyWavenumber = this%nkfmax
	end function getNkfmaxGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nwnmax
!
	integer function getNwnmaxGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNwnmaxGreenFrequencyWavenumber = this%nwnmax
	end function getNwnmaxGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nwn for given frequency index
!
	integer function getNwnSelectedGreenFrequencyWavenumber(this,if)
	type (green_frequency_wavenumber), intent(in) :: this
	integer, intent(in) :: if
	getNwnSelectedGreenFrequencyWavenumber = this%nwn(if)
	end function getNwnSelectedGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get nf
!
	integer function getNfGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNfGreenFrequencyWavenumber = this%nf2-this%nf1+1
	end function getNfGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get number of radial source nodes
!
	integer function getNnodGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getNnodGreenFrequencyWavenumber = .nnod.(this%gem)
	end function getNnodGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get frequency spacing
!
	real function getDfGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getDfGreenFrequencyWavenumber = this%df
	end function getDfGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get wavenumber spacing
!
	real function getDwnGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getDwnGreenFrequencyWavenumber = this%dwn
	end function getDwnGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get sigma
!
	real function getSigmaGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getSigmaGreenFrequencyWavenumber = this%sigma
	end function getSigmaGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get first complex frequency
!
	complex function getZomaGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getZomaGreenFrequencyWavenumber = cmplx(2.*mc_pi*(this%nf1-1)*this%df,-this%sigma)
	end function getZomaGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get rearth
!
	real function getRerdeGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getRerdeGreenFrequencyWavenumber = .rearth.(this%gem)
	end function getRerdeGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get ze
!
	real function getZeGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getZeGreenFrequencyWavenumber = this%ze
	end function getZeGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get current radius for which spectra have been read in
!
	double precision function getRcurGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber), intent(in) :: this
	getRcurGreenFrequencyWavenumber = this%rcur
	end function getRcurGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  derivatives available ?
!
	logical function isDerivativesGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber) :: this
	isDerivativesGreenFrequencyWavenumber = (this%derivflag == 1)
	end function isDerivativesGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  DSV stepping
!
	integer function dsvSteppingGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber) :: this
	dsvSteppingGreenFrequencyWavenumber = this%dsvstep
	end function dsvSteppingGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  calculate index of radial node closest to given depth
!
	integer function radialNodeIndexGreenFrequencyWavenumber(this,zs) result(jr)
	type (green_frequency_wavenumber) :: this
	real :: zs
	jr = radialNodeIndexGeminiEarthModel(this%gem,dble(zs))
	end function radialNodeIndexGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get total number of non-zero spectra
!
	integer function getNspGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber) :: this
	getNspGreenFrequencyWavenumber = sum(this%dsvmask_sph)*this%nsrcsph+sum(this%dsvmask_tor)*this%nsrctor
	end function getNspGreenFrequencyWavenumber
!-------------------------------------------------------------------------
!  get linear spectral index form DSV-index (isp) amd jump index (jsp)
!
	integer function getLinearSpectralCountGreenFrequencyWavenumber(this,isp,jsp)
	type (green_frequency_wavenumber) :: this
	integer :: isp,jsp
	getLinearSpectralCountGreenFrequencyWavenumber = this%ksp(isp,jsp)
	end function getLinearSpectralCountGreenFrequencyWavenumber
!--------------------------------------------------------------------------
!> \brief Get geminiEarthModel
!
	function getGeminiEarthMomdelGreenFrequencyWavenumber(this) result(gem)
	type (green_frequency_wavenumber), intent(in) :: this
	type (gemini_earth_model) :: gem
	gem = this%gem
	end function getGeminiEarthMomdelGreenFrequencyWavenumber
!--------------------------------------------------------------------------
!  print table of contents of FK spectra file
!
	subroutine tableContentsGreenFrequencyWavenumber(this)
	type (green_frequency_wavenumber) :: this
	integer :: isp,jsp
	character (len=2) :: cdsv
!
	write(6,'(a)') 'Table of contents:'
	write(6,'(a6,a8,a8,a8)') 'Name','DSV','Jump','Count'
	do isp = 1,9
		cdsv = getNameDsvFrequencyWavenumberMapping(isp)
		do jsp = 1,4
			if (this%ksp(isp,jsp) > 0) write(6,'(a6,i8,i8,i8)') cdsv,isp,jsp,this%ksp(isp,jsp)
		enddo
	enddo
	end subroutine tableContentsGreenFrequencyWavenumber
!
 end module
