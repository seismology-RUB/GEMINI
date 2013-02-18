c--------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
c----------------------------------------------------------------------
c  compute DSV jump solutions for sources at radial nodes and receiver at re
c
	subroutine jump_solution_tor(wn,f,sigma,ze,gfe,nsrc,zsrc,iprint,alternative,derivflag)
	include 'pi.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	double precision wn,f,sigma,ze
	integer iprint,nsrc,nnsol,nb,derivflag
	double precision om,elp1,re,rst,rearth,ratio
	double complex zom
	double complex wg1,tg1,wg2,tg2
	double complex gfe(nnd,3,2),gh(2,2,nnd),ghe(2,2,nnd),zsrc(2,2,nnd),det(2)
	integer nle,nlstrt,nr,i,j,top
	logical alternative
	common/degree/elp1
	common/omega/zom,om
c
	if(iprint.gt.0) write(6,'(a,f10.3)') ' wn/(2*pi) = ',wn/(2.*pi)
	om=2.*pi*f
	zom=dcmplx(om,-sigma)
	rearth=earthmodel_rearth
	elp1=(rearth*wn)**2
	re=rearth-ze
c
c  get layer of re
c
	call layer_getindex(re,nle,top)
c
c  zero gfe
c
	do j=1,nsrc
		do i=1,3
			do nr=1,nnod
				gfe(nr,i,j)=0.d0
			enddo
		enddo
	enddo
c
c  return if wn=0, gfe = 0
c
	if(wn.lt.1.d-8) return
c
c  return if source or receiver are in a top fluid layer
c
	call layer_isfluid(nle,i)
	if(i.eq.1) return
c
c  compute starting radius
c  for a global application, more help to find the starting radius
c  is needed (see lmspher)
c
	call startr(rnod(1),rst,nlstrt,iprint)
	if(iprint.gt.0) print *,'<greenkfsph>: rst,nlstrt :',rst,nlstrt
c
c  if starting radius is above receiver
c  stop working , because motion is observed for this k,om.
c  The starting radius may not imply a zero solution there if
c  it is in the halfspace. But startr.f ensures that the starting radius
c  is below the deepest node in that case which may also be source or receiver.
c  Thus, if rst happens to be above re it does imply a zero solution there
c  and we can safely skip the integration.

	if (rst.ge.re) return
c
c  return if starting radius is in a top fluid layer
c
	call layer_isfluid(nlstrt,i)
	if(i.eq.1) return
c
c  From now on we know that the starting radius is below receiver c
c
c  integrate toroidal equation from starting radius up to re,re and down from
c  surface to re,re. nnsol = jwd if a water layer exists
c  gh(W or T,below-above-receiver,nr)
c
	call torint(nle,re,nle,re,nlstrt,rst,gh,det,nnsol,iprint)
c
c  Calculate a solution at depth nodes for a source at re
c
	if (alternative) then
		do j=1,nsrc
c
c  nodes below re
c
			do i=1,2
				do nr=jstu,jed
					gfe(nr,i,j)=(-zsrc(1,j,1)*gh(2,2,jeu)+zsrc(2,j,1)*gh(1,2,jeu))*gh(i,1,nr)/det(2)
				enddo
			enddo
c
c  nodes above re
c
			do i=1,2
				do nr=jeu,nnsol
					gfe(nr,i,j)=(-zsrc(1,j,1)*gh(2,1,jed)+zsrc(2,j,1)*gh(1,1,jed))*gh(i,2,nr)/det(2)
				enddo
			enddo
			if (derivflag == 1) call torderivs(gfe(1,1,j))
		enddo
	else
c
c  Calculate a solution for a receiver at re for sources at depth nodes
c
c  below receiver
c
		do nr=jstu,jed
			ratio=(rnod(nr)/re)**2
c
c  construct two basis Green solutions for jumps (1,0) and (0,1) at re
c
			wg1 = -gh(2,2,jeu)/det(2)*gh(1,1,nr)       ! jump (1,0), DSV: W
			tg1 = -gh(2,2,jeu)/det(2)*gh(2,1,nr)       ! jump (1,0), DSV: T
			wg2 = +gh(1,2,jeu)/det(2)*gh(1,1,nr)       ! jump (0,1), DSV: W
			tg2 = +gh(1,2,jeu)/det(2)*gh(2,1,nr)       ! jump (0,1), DSV: T
!			print *,nr,rearth-rnod(nr),'W1(r,re) = ',ratio*abs(wg1),' T1(r,re) = ',ratio*abs(tg1),
!    1		                                  'W2(r,re) = ',ratio*abs(wg2),' T2(r,re) = ',ratio*abs(tg2)
!
c  Construct two basis solutions for a source at rnod(nr) and receiver at re
c  using the reciprocity relation
c  first index: basis solution, second index: component
! 
			ghe(1,1,nr) = -ratio*tg2         !  W1(r,rs) with jump = (1,0) at rs
			ghe(1,2,nr) = +ratio*tg1         !  T1(r,rs) with jump = (1,0) at rs
			ghe(2,1,nr) = +ratio*wg2         !  W2(r,rs) with jump = (0,1) at rs
			ghe(2,2,nr) = -ratio*wg1         !  T2(r,rs) with jump = (0,1) at rs
		enddo
c
c  same procedure above receiver
c
		do nr=jeu,nnsol
			ratio=(rnod(nr)/re)**2
			wg1 = -gh(2,1,jed)/det(2)*gh(1,2,nr)
			tg1 = -gh(2,1,jed)/det(2)*gh(2,2,nr)
			wg2 = +gh(1,1,jed)/det(2)*gh(1,2,nr)
			tg2 = +gh(1,1,jed)/det(2)*gh(2,2,nr)
			ghe(1,1,nr) = -ratio*tg2
			ghe(1,2,nr) = +ratio*tg1
			ghe(2,1,nr) = +ratio*wg2
			ghe(2,2,nr) = -ratio*wg1
		enddo
c
c  linearly combine basis solutions using prescribed jumps at source
c  nb: basis solution, jump vector = (1,0) etc
c  j:  different jump vectors needed for calculating seismograms for force or moment tensor sources
c  i:  DSV component (W,T)
c
		do i=1,2
			do j=1,nsrc
				do nr=jstu,nnsol
					gfe(nr,i,j)=0.d0
					do nb=1,2
						gfe(nr,i,j)=gfe(nr,i,j)+zsrc(nb,j,nr)*ghe(nb,i,nr)
					enddo
				enddo
			enddo
		enddo
	endif
	return
	end

