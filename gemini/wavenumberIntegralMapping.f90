!----------------------------------------------------------------------------------------
!  Module with mapping from wavenumber integral indexing to motion type, DSV component
!  an jump vector index
!-----------------------------------------------------------------------------------------
 module wavenumberIntegralMapping
	implicit none
!
!  current max number of wavenumber integrals for force excitation (twice as much for moment tensor)
!
	integer, parameter :: maxWavenumberIntegrals = 23
!
!  Frequency-wavenumber spectra ordered by two indices: first gives DSV and second jump,
!  DSV in the following order: U,R,V,S,W,T,(UP,VP,WP) 
!  Mapping from first fk-spectra index to DSV component and motion
!
	integer, dimension(9) :: fksp_dsv =  (/ 1,2,3,4,1,2,5,6,3 /)
	character (len=3), dimension(9) :: fksp_motion = &
		& (/'sph','sph','sph','sph','tor','tor','sph','sph','tor' /)
	character (len=2), dimension(9) :: fksp_cdsv = (/' U',' R',' V',' S',' W',' T','UP','VP','WP' /) 
!
!                                 FORCE EXCITATION
!
!  Mapping from the 23 wavenumber integrals to the two frequency-wavenumber indices:
!
	integer, dimension(23,2) :: wnint_fksp_f = reshape( &
		& (/ 1,1,3,3,5,5,3, 7,7,8,8,9,9,8, 1,1,3,3,5,5,3, 2,2, &
	       &    1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2 /) &
		& ,(/ 23,2 /))
!                                 MOMENT EXCITATION
!
!  Mapping from the 46 wavenumber integrals to the two frequency-wavenumber indices:
!
	integer, dimension(46,2) :: wnint_fksp_m = reshape( &
		 & (/ 1,1,1,1,3,3,3,3,5,5,5,5,3,3, 7,7,7,7,8,8,8,8,9,9,9,9,8,8, 1,1,1,1,3,3,3,3,5,5,5,5,3,3, 2,2,2,2,  &
		 &    1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4 /) &
		 & ,(/ 46,2 /))
!
	private :: fksp_dsv,fksp_motion,wnint_fksp_f,wnint_fksp_m,maxWavenumberIntegrals
!
 contains
!-----------------------------------------------------------------------------------------
!  given first frequency-wavenumber spectrum index, return DSV index
!
	integer function getDsvFrequencyWavenumberMapping(isp) result(idsv)
	integer :: isp
	idsv = fksp_dsv(isp)
	end function getDsvFrequencyWavenumberMapping
!-----------------------------------------------------------------------------------------
!  given first frequency-wavenumber spectrum index, return DSV name
!
	character (len=2) function getNameDsvFrequencyWavenumberMapping(isp) result(cdsv)
	integer :: isp
	cdsv = fksp_cdsv(isp)
	end function getNameDsvFrequencyWavenumberMapping
!----------------------------------------------------------------------------------
!  given first frequency wavenumber spectrum index, return motion type
!
	function getMotionFrequencyWavenumberMapping(isp) result(motion)
	character (len=3) :: motion
	integer :: isp
	motion = fksp_motion(isp)
	end function getMotionFrequencyWavenumberMapping
!-----------------------------------------------------------------------------------------------
!  given wavenumber integral index and source type, return frequency-wavenumber spectrum indices 
!
	subroutine getIspJspWavenumberIntegralMapping(kint,istyp,isp,jsp)
	integer :: kint,istyp,isp,jsp
	if (istyp == 0) then
		isp = wnint_fksp_f(kint,1)
		jsp = wnint_fksp_f(kint,2)
	else if(istyp == 1) then
		isp = wnint_fksp_m(kint,1)
		jsp = wnint_fksp_m(kint,2)
	else
		print *,'getKintFrequencyWavenumberMapping: invalid source type'
		stop
	endif
	end subroutine getIspJspWavenumberIntegralMapping
!----------------------------------------------------------------------------------------
!  get max number of wavenumber integrals
!
	integer function getMaxWavenumberIntegrals()
	getMaxWavenumberIntegrals = maxWavenumberIntegrals
	end function getMaxWavenumberIntegrals
!
 end module wavenumberIntegralMapping
