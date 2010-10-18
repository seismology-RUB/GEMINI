!-------------------------------------------------------------------
!  some utilities associated with plotting
!-------------------------------------------------------------------
 module plotUtilities
	implicit none
!
 contains
!-----------------------------------------------------------
!  Initialize a pgplot graph
!
	subroutine initPgPlotGraph(xmin,xmax,ymin,ymax,eqar,ch,xlab,ylab,title)
	real :: xmin,xmax,ymin,ymax
	real, optional :: ch
	logical :: eqar
	character (len=*), optional :: xlab,ylab,title
	character (len=132) :: xlabel,ylabel
!
	call pgsci(1);call pgscf(2)
	call pgask(.false.)
	call pgeras
	call pgsvp(0.1,0.9,0.1,0.9)
	if (.not. eqar) then
		call pgswin(xmin,xmax,ymin,ymax)
	else
		call pgwnad(xmin,xmax,ymin,ymax)
	endif
	call pgslw(3)
	if (present(ch)) then; call pgsch(ch); else; call pgsch(1.0); endif
	call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
	if (present(xlab)) then; xlabel = xlab; else; xlabel = ' '; endif
	if (present(ylab)) then; ylabel = ylab; else; ylabel = ' '; endif
	call pglab(trim(xlabel),trim(ylabel),' ')
	if (present(title)) call pgtext(0.3*xmax+0.7*xmin,1.05*ymax,trim(title))
	call pgslw(1)
	end subroutine initPgPlotGraph
!-----------------------------------------------------------
!  Draw a straight line connecting two points
!
	subroutine connectTwoPoints(x1,y1,x2,y2,ls)
	real :: x1,y1,x2,y2
	integer, optional :: ls
	if (present(ls)) call pgsls(ls)
	call pgpt1(x1,y1,-1)
	call pgdraw(x2,y2)
	call pgsls(1)
	end subroutine connectTwoPoints
!-----------------------------------------------------------
!  convert axis-parameter of pgenv to axis string of pgbox
!
	character (len=10) function convertAxisToBox(iaxis)
	integer :: iaxis
	select case (iaxis)
		case (-2); convertAxisToBox = ''
		case (-1); convertAxisToBox = 'BC'
		case (0);  convertAxisToBox = 'BCNST'
		case (1);  convertAxisToBox = 'ABCNST'
		case (2);  convertAxisToBox = 'ABCGNST'
		case default; convertAxisToBox ='BCNST' 
	end select
	end function convertAxisToBox
!--------------------------------------------------------------------
!  Draw edges of a spherical slice
!
!  np: number of points used for circular parts
!  alf1,alf2: angular limits of slice in radians counted math positive from x-axis
!  r1,r2: radii of inner and outer circle of slice
!  xcenter, ycenter: world coordinates of center of slice
!
	subroutine drawEdgesSphericalSlice(np,alf1,alf2,r1,r2,xcenter,ycenter)
	integer :: np
	real :: alf1,alf2,r1,r2,xcenter,ycenter
	real, dimension(:) :: x(2),y(2)
!
	call drawPartialCircle(np,alf1,alf2,r1,xcenter,ycenter)
	x(1) = xcenter+r1*cos(alf2); y(1) = ycenter+r1*sin(alf2)
	x(2) = xcenter+r2*cos(alf2); y(2) = ycenter+r2*sin(alf2)
	call pgline(2,x,y)
	call drawPartialCircle(np,alf1,alf2,r2,xcenter,ycenter)
	x(1) = xcenter+r1*cos(alf1); y(1) = ycenter+r1*sin(alf1)
	x(2) = xcenter+r2*cos(alf1); y(2) = ycenter+r2*sin(alf1)
	call pgline(2,x,y)
	end subroutine drawEdgesSphericalSlice
!--------------------------------------------------------------------
!  Draw part of a circle
!
!  number of points used for plotting circle
!  alf1,alf2: angular limits of circle, counted math positive from x-axis, in radians
!  r: radius of circle
!  xcenter,ycenter: world coordinates of center
!
	subroutine drawPartialCircle(np,alf1,alf2,r,xcenter,ycenter)
	integer :: np,i
	real :: alf1,alf2,alfa,dalfa,r,xcenter,ycenter
	real, dimension(:), allocatable :: x,y
!
	allocate(x(np),y(np))
	dalfa = (alf2-alf1)/(np-1)
	alfa = alf1
	do i=1,np
		x(i) = xcenter+r*cos(alfa)
		y(i) = ycenter+r*sin(alfa)
		alfa = alfa+dalfa
	enddo
	call pgline(np,x,y)
	deallocate(x,y)
	end subroutine drawPartialCircle
!--------------------------------------------------------------------
!  Draw an ellipse (also partial)
!
!  number of points used for plotting ellipse
!  alf1,alf2: angular limits of ellipse, counted math positive from x-axis, in radians
!  rx, ry: half axes of ellipse
!  xcenter,ycenter: world coordinates of center
!
	subroutine drawPartialEllipse(np,alf1,alf2,rx,ry,xcenter,ycenter)
	integer :: np,i
	real :: alf1,alf2,alfa,dalfa,rx,ry,xcenter,ycenter
	real, dimension(:), allocatable :: x,y
!
	allocate(x(np),y(np))
	dalfa = (alf2-alf1)/(np-1)
	alfa = alf1
	do i=1,np
		x(i) = xcenter+rx*cos(alfa)
		y(i) = ycenter+ry*sin(alfa)
		alfa = alfa+dalfa
	enddo
	call pgline(np,x,y)
	deallocate(x,y)
	end subroutine drawPartialEllipse
!----------------------------------------------------------------
 end module plotUtilities
