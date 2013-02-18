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
!  Module with routines to calculate displacement spectra from
!  wavenumber integrals
!----------------------------------------------------------------------
 module greenDisplacementSpectrum
	implicit none
!
 contains
!---------------------------------------------------------------
!  displacement spectrum for single force and one component
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,N,E,H,R,T)
!  azi:  propagation direction
!  phi:  receiver azimuth
!
	subroutine forceGreenDisplacementSpectrum(nf1,nf2,csys,comp,azi,phi,force,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: csys,comp,cdum
	real :: azi,phi,caz,saz,sphi,cphi
	real, dimension(:) :: force            !  force components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  7 wavenumber integrals vs frequency
	complex :: zspl,zspt,zspn,zspe
!
	cdum = csys
	caz = cos(azi); saz = sin(azi)
	sphi = sin(phi); cphi = cos(phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,1)*force(1)+zdis(if,2)*(force(2)*cphi+force(3)*sphi)
		enddo
	else if(comp == 'H') then
		do if=nf1,nf2
			zsp(if) = zdis(if,22)*force(1)+zdis(if,23)*(force(2)*cphi+force(3)*sphi)
		enddo
	else if(comp == 'N' .or. comp == 'E' .or. comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl = zdis(if,3)*force(1) &
			    & +(zdis(if,4)+zdis(if,5))*(force(2)*cphi+force(3)*sphi)
!
			zspt = (zdis(if,6)+zdis(if,7))*(force(2)*sphi-force(3)*cphi)
!
! if x is pointing south and y pointing east and z pointing up
! there is no difference for spherical or cartesian coordinate system
!
			zspn=-zspl*caz+zspt*saz
			zspe=+zspl*saz+zspt*caz
!
			if(comp == 'N') zsp(if)=zspn
			if(comp == 'E') zsp(if)=zspe
			if(comp == 'R') zsp(if)=zspl
			if(comp == 'T') zsp(if)=zspt
		enddo
   	endif
	end subroutine forceGreenDisplacementSpectrum
!---------------------------------------------------------------
!  radial derivative of displacement spectrum for single force and one component
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  phi:  receiver azimuth
!
	subroutine forceRadDerivGreenDisplacementSpectrum(nf1,nf2,comp,phi,force,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: phi,sphi,cphi
	real, dimension(:) :: force            !  force components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  7 wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,8)*force(1)+zdis(if,9)*(force(2)*cphi+force(3)*sphi)
		enddo
	else if(comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl = zdis(if,10)*force(1) &
			    & +(zdis(if,11)+zdis(if,12))*(force(2)*cphi+force(3)*sphi)
			zspt = (zdis(if,13)+zdis(if,14))*(force(2)*sphi-force(3)*cphi)
			if(comp == 'R') zsp(if)=zspl
			if(comp == 'T') zsp(if)=zspt
		enddo
   	endif
	end subroutine forceRadDerivGreenDisplacementSpectrum
!---------------------------------------------------------------
!  1/r*theta derivative of displacement spectrum for single force and one component
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  r:    receiver radius
!  phi:  receiver azimuth
!
	subroutine forceThetaDerivGreenDisplacementSpectrum(nf1,nf2,comp,r,phi,force,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: r,phi,sphi,cphi
	real, dimension(:) :: force            !  force components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  7 wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,15)*force(1)+zdis(if,16)*(force(2)*cphi+force(3)*sphi)
			zsp(if) = zsp(if)/r
		enddo
	else if(comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl = zdis(if,17)*force(1) &
			    & +(zdis(if,18)+zdis(if,19))*(force(2)*cphi+force(3)*sphi)
			zspt = (zdis(if,20)+zdis(if,21))*(force(2)*sphi-force(3)*cphi)
			if(comp == 'R') zsp(if)=zspl/r
			if(comp == 'T') zsp(if)=zspt/r
		enddo
   	endif
	end subroutine forceThetaDerivGreenDisplacementSpectrum
!---------------------------------------------------------------
!  1/(r*sin(theta))*phi-derivative of displacement spectrum for single force and one component
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  rst:  r*sin(theta) in csys=S or dis in csys = C
!  phi:  receiver azimuth
!
	subroutine forcePhiDerivGreenDisplacementSpectrum(nf1,nf2,comp,rst,phi,force,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: rst,phi,sphi,cphi
	real, dimension(:) :: force            !  force components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  7 wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,2)*(-force(2)*sphi+force(3)*cphi)
			zsp(if) = zsp(if)/rst
		enddo
	else if(comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl = (zdis(if,4)+zdis(if,5))*(-force(2)*sphi+force(3)*cphi)
			zspt = (zdis(if,6)+zdis(if,7))*(+force(2)*cphi+force(3)*sphi)
			if(comp == 'R') zsp(if)=zspl/rst
			if(comp == 'T') zsp(if)=zspt/rst
		enddo
   	endif
	end subroutine forcePhiDerivGreenDisplacementSpectrum
!--------------------------------------------------------------------
!  displacement spectrum for moment tensor
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,N,E,H,R,T)
!  azi:  propagation direction
!  phi:  receiver azimuth
!
	subroutine momentGreenDisplacementSpectrum(nf1,nf2,csys,comp,azi,phi,mt,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: csys,comp,cdum
	real :: azi,phi,caz,saz,sphi,cphi,c2phi,s2phi
	real, dimension(:) :: mt               !  moment tensor components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  wavenumber integrals vs frequency
	complex :: zspl,zspt,zspn,zspe
!
	cdum = csys
	caz = cos(azi); saz = sin(azi)
	sphi = sin(phi); cphi = cos(phi)
	s2phi = sin(2.*phi); c2phi = cos(2.*phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,1)* mt(1) &
				& +zdis(if,2)*(mt(2)+mt(3)) &
				& +zdis(if,3)*(mt(4)*cphi+mt(5)*sphi) &
				& +zdis(if,4)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
		enddo
	else if(comp == 'H') then
		do if=nf1,nf2
			zsp(if) = zdis(if,43)* mt(1) &
				& +zdis(if,44)*(mt(2)+mt(3)) &
				& +zdis(if,45)*(mt(4)*cphi+mt(5)*sphi) &
				& +zdis(if,46)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
		enddo
	else if(comp == 'N' .or. comp == 'E' .or. comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl =    zdis(if,5)*mt(1) &
				& +zdis(if,6)*(mt(2)+mt(3)) &
				& +(zdis(if,7)+zdis(if,9))*(mt(4)*cphi+mt(5)*sphi) &
				& +(zdis(if,8)+zdis(if,10))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
!
			zspt =    (zdis(if,11)+zdis(if,13))*(mt(5)*cphi-mt(4)*sphi) &
				& +(zdis(if,12)+zdis(if,14))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
!
!
! if x is pointing south and y pointing east and z pointing up
! there is no difference for spherical or cartesian coordinate system
!
			zspn=-zspl*caz+zspt*saz
			zspe=+zspl*saz+zspt*caz
!
			if(comp == 'N') zsp(if)=zspn
			if(comp == 'E') zsp(if)=zspe
			if (comp == 'R') zsp(if)=zspl
			if (comp == 'T') zsp(if)=zspt
		enddo
	endif
	end subroutine momentGreenDisplacementSpectrum
!--------------------------------------------------------------------
!  radial derivatives of displacement spectrum for moment tensor
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  phi:  receiver azimuth from south
!
	subroutine momentRadDerivGreenDisplacementSpectrum(nf1,nf2,comp,phi,mt,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: phi,sphi,cphi,c2phi,s2phi
	real, dimension(:) :: mt               !  moment tensor components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
	s2phi = sin(2.*phi); c2phi = cos(2.*phi)
!
	if (comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,15)* mt(1) &
				& +zdis(if,16)*(mt(2)+mt(3)) &
				& +zdis(if,17)*(mt(4)*cphi+mt(5)*sphi) &
				& +zdis(if,18)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
		enddo
	else if (comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl =    zdis(if,19)*mt(1) &
				& +zdis(if,20)*(mt(2)+mt(3)) &
				& +(zdis(if,21)+zdis(if,23))*(mt(4)*cphi+mt(5)*sphi) &
				& +(zdis(if,22)+zdis(if,24))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
!
			zspt =    (zdis(if,25)+zdis(if,27))*(mt(5)*cphi-mt(4)*sphi) &
				& +(zdis(if,26)+zdis(if,28))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
			if (comp == 'R') zsp(if)=zspl
			if (comp == 'T') zsp(if)=zspt
		enddo
	endif
	end subroutine momentRadDerivGreenDisplacementSpectrum
!--------------------------------------------------------------------
!  1/r*theta derivatives of displacement spectrum for moment tensor
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  r:    radius of receiver
!  phi:  receiver azimuth from south
!
	subroutine momentThetaDerivGreenDisplacementSpectrum(nf1,nf2,comp,r,phi,mt,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: r,phi,sphi,cphi,c2phi,s2phi
	real, dimension(:) :: mt               !  moment tensor components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
	s2phi = sin(2.*phi); c2phi = cos(2.*phi)
!
	if (comp == 'Z') then
		do if=nf1,nf2
			zsp(if) = zdis(if,29)* mt(1) &
				& +zdis(if,30)*(mt(2)+mt(3)) &
				& +zdis(if,31)*(mt(4)*cphi+mt(5)*sphi) &
				& +zdis(if,32)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
			zsp(if) = zsp(if)/r
		enddo
	else if (comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl =    zdis(if,33)*mt(1) &
				& +zdis(if,34)*(mt(2)+mt(3)) &
				& +(zdis(if,35)+zdis(if,37))*(mt(4)*cphi+mt(5)*sphi) &
				& +(zdis(if,36)+zdis(if,38))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
!
			zspt =    (zdis(if,39)+zdis(if,41))*(mt(5)*cphi-mt(4)*sphi) &
				& +(zdis(if,40)+zdis(if,42))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
			if (comp == 'R') zsp(if)=zspl/r
			if (comp == 'T') zsp(if)=zspt/r
		enddo
	endif
	end subroutine momentThetaDerivGreenDisplacementSpectrum
!--------------------------------------------------------------------
!  1/(r*sin(theta))*phi derivatives of displacement spectrum for moment tensor
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,R,T)
!  rst:  either r*sin(theta) in csys=S or dis when csys = C
!  phi:  receiver azimuth from south
!
	subroutine momentPhiDerivGreenDisplacementSpectrum(nf1,nf2,comp,rst,phi,mt,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: comp
	real :: rst,phi,sphi,cphi,c2phi,s2phi
	real, dimension(:) :: mt               !  moment tensor components
	complex, dimension(:) :: zsp           !  displacement spectrum
	complex, dimension(:,:) :: zdis        !  wavenumber integrals vs frequency
	complex :: zspl,zspt
!
	sphi = sin(phi); cphi = cos(phi)
	s2phi = sin(2.*phi); c2phi = cos(2.*phi)
!
	if (comp == 'Z') then
		do if=nf1,nf2
			zsp(if) =    zdis(if,3)*(-mt(4)*sphi+mt(5)*cphi) &
				   & +zdis(if,4)*(-2.*(mt(3)-mt(2))*s2phi-4.*mt(6)*c2phi)
			zsp(if) = zsp(if)/rst
		enddo
	else if (comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl =    (zdis(if,7)+zdis(if,9))*(-mt(4)*sphi+mt(5)*cphi) &
				& +(zdis(if,8)+zdis(if,10))*(-2.*(mt(3)-mt(2))*s2phi-4.*mt(6)*c2phi)
!
			zspt =    (zdis(if,11)+zdis(if,13))*(-mt(5)*sphi-mt(4)*cphi) &
				& +(zdis(if,12)+zdis(if,14))*(+2.*(mt(3)-mt(2))*c2phi-4.*mt(6)*s2phi)
			if (comp == 'R') zsp(if)=zspl/rst
			if (comp == 'T') zsp(if)=zspt/rst
		enddo
	endif
	end subroutine momentPhiDerivGreenDisplacementSpectrum
!---------------------------------------------------------------
!  2D displacement spectrum for single force and one component
!  nf1,nf2: start and end indices for frequency
!  comp: component (Z,N,E,H,R,T)
!  azi:  propagation direction
!  phi:  receiver azimuth
!
	subroutine force2DGreenDisplacementSpectrum(nf1,nf2,comp,azi,force,zdis,zsp)
	integer :: nf1,nf2
	character(len=1) :: comp
	real :: azi,caz
	real, dimension(:) :: force                !  force components
	complex, dimension(:) :: zsp               !  displacement spectrum
	complex, dimension(:,:) :: zdis            !  7 (2 unused) wavenumber integrals vs frequency
!
	caz = cos(azi)
	if (comp == 'Z') zsp(nf1:nf2) = zdis(nf1:nf2,1)*force(1)+zdis(nf1:nf2,2)*force(2)
	if (comp == 'H') zsp(nf1:nf2) = zdis(nf1:nf2,22)*force(1)+zdis(nf1:nf2,23)*force(2)
	if (comp == 'R') zsp(nf1:nf2) = zdis(nf1:nf2,3)*force(1)+zdis(nf1:nf2,4)*force(2)
	if (comp == 'T') zsp(nf1:nf2) = zdis(nf1:nf2,6)*(-force(3))
	if (comp == 'N') zsp(nf1:nf2) = (zdis(nf1:nf2,3)*force(1)+zdis(nf1:nf2,4)*force(2))*caz
	if (comp == 'E') zsp(nf1:nf2) = zdis(nf1:nf2,6)*(-force(3))*caz
	end subroutine force2DGreenDisplacementSpectrum
!--------------------------------------------------------------------
!  6 basis moment tensor displacement spectra
!  Basis moment tensor components are:
!
!	1., 1., 1., 0., 0., 0.
!     -1., 1., 0., 0., 0., 0.
!     -1., 0., 1., 0., 0., 0.
!      0., 0., 0., 0., 0., 1.
!      0., 0., 0., 0., 1., 0.
!      0., 0., 0., 1., 0., 0.
!
	subroutine basemomGreenDisplacementSpectrum(nf1,nf2,csys,comp,azi,phi,zdis,zsp)
	integer :: nf1,nf2,if
	character(len=1) :: csys,comp,cdum
	real :: azi,phi,caz,saz,sphi,cphi,c2phi,s2phi
	complex, dimension(:,:) :: zsp                  !  6 displacement spectra
	complex, dimension(:,:) :: zdis                 !  14 wavenumber integrals vs frequency
	complex, dimension(6) :: zspl,zspt,zspn,zspe
!
	cdum = csys
	caz = cos(azi); saz = sin(azi)
	sphi = sin(phi); cphi = cos(phi)
	s2phi = sin(2.*phi); c2phi = cos(2.*phi)
!
	if(comp == 'Z') then
		do if=nf1,nf2
			zsp(if,1) = zdis(if,1)+2.*zdis(if,2)                     ! +1,1,1,0,0,0
			zsp(if,2) = -zdis(if,1)+zdis(if,2)-zdis(if,4)*c2phi      ! -1,1,0,0,0,0
			zsp(if,3) = -zdis(if,1)+zdis(if,2)+zdis(if,4)*c2phi      ! -1,0,1,0,0,0
			zsp(if,4) = -2.*zdis(if,4)*s2phi                         !  0,0,0,0,0,1
			zsp(if,5) = zdis(if,3)*sphi                              !  0,0,0,0,1,0
			zsp(if,6) = zdis(if,3)*cphi                              !  0,0,0,1,0,0
		enddo
	else if(comp == 'H') then
		do if=nf1,nf2
			zsp(if,1) = zdis(if,43)+2.*zdis(if,44)                     ! +1,1,1,0,0,0
			zsp(if,2) = -zdis(if,43)+zdis(if,44)-zdis(if,46)*c2phi     ! -1,1,0,0,0,0
			zsp(if,3) = -zdis(if,43)+zdis(if,44)+zdis(if,46)*c2phi     ! -1,0,1,0,0,0
			zsp(if,4) = -2.*zdis(if,46)*s2phi                          !  0,0,0,0,0,1
			zsp(if,5) = zdis(if,45)*sphi                               !  0,0,0,0,1,0
			zsp(if,6) = zdis(if,45)*cphi                               !  0,0,0,1,0,0
		enddo
	else if(comp == 'N' .or. comp == 'E' .or. comp == 'R' .or. comp == 'T') then
		do if=nf1,nf2
			zspl(1) = zdis(if,5)+2.*zdis(if,6)
			zspl(2) = -zdis(if,5)+zdis(if,6)-(zdis(if,8)+zdis(if,10))*c2phi
			zspl(3) = -zdis(if,5)+zdis(if,6)+(zdis(if,8)+zdis(if,10))*c2phi
			zspl(4) = -2.*(zdis(if,8)+zdis(if,10))*s2phi
			zspl(5) = (zdis(if,7)+zdis(if,9))*sphi
			zspl(6) = (zdis(if,7)+zdis(if,9))*cphi
			zspt(1) = 0.
			zspt(2) = -(zdis(if,12)+zdis(if,14))*s2phi
			zspt(3) = +(zdis(if,12)+zdis(if,14))*s2phi
			zspt(4) = 2.*(zdis(if,12)+zdis(if,14))*c2phi
			zspt(5) = +(zdis(if,11)+zdis(if,13))*cphi
			zspt(6) = -(zdis(if,11)+zdis(if,13))*sphi
!
! if x is pointing south and y pointing east and z pointing up
! there is no difference for spherical or cartesian coordinate system
!
			zspn(:)=-zspl(:)*caz+zspt(:)*saz
			zspe(:)=+zspl(:)*saz+zspt(:)*caz
!
			if(comp == 'N') zsp(if,:)=zspn(:)
			if(comp == 'E') zsp(if,:)=zspe(:)
			if (comp == 'R') zsp(if,:)=zspl(:)
			if (comp == 'T') zsp(if,:)=zspt(:)
		enddo
	endif
	end subroutine basemomGreenDisplacementSpectrum
!---------------------------------------------------------------
!  5 basis displacement spectra for single force excitation
!  nf1,nf2: start and end indices for frequency
!
	subroutine forceBasisGreenDisplacementSpectrum(nf1,nf2,zdis,zsp)
	integer :: nf1,nf2,if
	complex, dimension(:,:) :: zsp         !  5 displacement spectra
	complex, dimension(:,:) :: zdis        !  7 wavenumber integrals vs frequency
!
	do if=nf1,nf2
		zsp(if,1) = zdis(if,1)
		zsp(if,2) = zdis(if,2)
		zsp(if,3) = zdis(if,3)
		zsp(if,4) = zdis(if,4)+zdis(if,5)
		zsp(if,5) = zdis(if,6)+zdis(if,7)
	enddo
	end subroutine forceBasisGreenDisplacementSpectrum
!---------------------------------------------------------------
!  10 basis displacement spectra for moment tensor excitation
!  nf1,nf2: start and end indices for frequency
!
	subroutine momentBasisGreenDisplacementSpectrum(nf1,nf2,zdis,zsp)
	integer :: nf1,nf2,if
	complex, dimension(:,:) :: zsp         !  10 displacement spectra
	complex, dimension(:,:) :: zdis        !  14 wavenumber integrals vs frequency
!
	do if=nf1,nf2
		zsp(if,1) = zdis(if,1)
		zsp(if,2) = zdis(if,2)
		zsp(if,3) = zdis(if,3)
		zsp(if,4) = zdis(if,4)
		zsp(if,5) = zdis(if,5)
		zsp(if,6) = zdis(if,6)
		zsp(if,7) = zdis(if,7)+zdis(if,9)
		zsp(if,8) = zdis(if,8)+zdis(if,10)
		zsp(if,9) = zdis(if,11)+zdis(if,13)
		zsp(if,10) = zdis(if,12)+zdis(if,14)
	enddo
	end subroutine momentBasisGreenDisplacementSpectrum
!-----------------------------------------------------------------------
!  convert displacment components and derivatives to global cartesion
!  coordinates
!
	subroutine convertToGlobalCartesianGreenDisplacementSpectrum(nf1,nf2,csys,tm,r,dis,rearth,zsp,zspdr,zspdt,zspdf,u,graducar)
	integer :: nf1,nf2,j,is,if,l,ir
	double precision :: r
	complex, dimension(:,:) :: zsp,zspdr,zspdt,zspdf
	complex, dimension(:,:) :: u
	complex, dimension(:,:,:) :: graducar                  ! graducar(freq,deriv,comp)
	complex, dimension(3,3) :: graduspher,work
	real, dimension(:,:) :: tm
	real :: dis,rearth,rtandel
	character (len=*) :: csys
!
	do if = nf1,nf2
!
!  transformation of displacement: U_j = M_js u_s 
!
		do j=1,3
			u(if,j) = 0.
			do is = 1,3
				u(if,j)=u(if,j)+tm(j,is)*zsp(if,is)
			enddo
		enddo
!
!  spherical gradient of displacement D_r u_s
!
		if(csys.eq.'S') then
			rtandel = r*tan(dis/rearth)
			graduspher(1,1) = zspdr(if,1)                       ! d_r u_r
			graduspher(1,2) = zspdr(if,2)                       ! d_r u_t
			graduspher(1,3) = zspdr(if,3)                       ! d_r u_f
			graduspher(2,1) = zspdt(if,1)-zsp(if,2)/r           ! 1/r d_t u_r - 1/r u_t
			graduspher(2,2) = zspdt(if,2)+zsp(if,1)/r           ! 1/r d_t u_t + 1/r u_r
			graduspher(2,3) = zspdt(if,3)                       ! 1/r d_t u_f
			graduspher(3,1) = zspdf(if,1)-zsp(if,3)/r           ! 1/(rst) d_f u_r-1/r u_f
			graduspher(3,2) = zspdf(if,2)-zsp(if,3)/rtandel     ! 1/(rst) d_f u_t-cot(theta)/r u_f
			graduspher(3,3) = zspdf(if,3)+zsp(if,1)/r+zsp(if,2)/rtandel   ! 1/rst d_f u_f+1/r u_r +cot(theta)/r u_t
		else
!
!  spherical gradient of displacement D_r u_s in shallow seismic approximation
!  spherical components are essentially treated as cylindrical ones
!  r -> inf, delta -> 0, r sin(delta) -> dis, cot(delta)/r -> 1/dis
!
			graduspher(1,1)=zspdr(if,1)                       ! d_r u_r
			graduspher(1,2)=zspdr(if,2)                       ! d_r u_t
			graduspher(1,3)=zspdr(if,3)                       ! d_r u_f
			graduspher(2,1)=zspdt(if,1)                       ! d_s u_r
			graduspher(2,2)=zspdt(if,2)                       ! d_s u_t
			graduspher(2,3)=zspdt(if,3)                       ! d_s u_f
			graduspher(3,1)=zspdf(if,1)                       ! 1/s d_f u_r
			graduspher(3,2)=zspdf(if,2)-zsp(if,3)/dis         ! 1/s d_f u_s -1/s u_f
			graduspher(3,3)=zspdf(if,3)+zsp(if,2)/dis         ! 1/s d_f u_f +1/s u_t
		endif
!
!  transformation of gradient d_j u_l = M_jr D_rs M_ls
!  first calculate D_rs M_ls
!
		do l = 1,3
			do ir = 1,3
				work(ir,l) = 0.
				do is = 1,3
					work(ir,l) = work(ir,l)+graduspher(ir,is)*tm(l,is)
				enddo
			enddo
		enddo
!
!  now calculate d_j u_l = M_jr * work(ir,l)  = graducar(if,j,l)
!
		do l = 1,3
			do j = 1,3
				graducar(if,j,l) = 0.
				do ir=1,3
					graducar(if,j,l) = graducar(if,j,l)+tm(j,ir)*work(ir,l)
				enddo
			enddo
		enddo
	enddo
	end subroutine convertToGlobalCartesianGreenDisplacementSpectrum
!
 end module greenDisplacementSpectrum
