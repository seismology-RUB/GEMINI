!----------------------------------------------------------------------------------
!  General bindings for selecting new extrema of a PgPlot
!  x1,x2,y1,y2:  new extremal values for plot range
!
	subroutine zoomBindingsPgPlotSelect(pgid,x1,y1,ch,x2,y2,xmin,xmax,ymin,ymax)
	real :: x1,y1,x2,y2,xmin,xmax,ymin,ymax
	character (len=1) :: ch,ch1
	integer :: pgid
!
	call pgslct(pgid)
	call pgslw(1);call pgsci(5)
	select case (ch)
	case ('A')                                                       !  zoom to rectangle, left mouse 2 times
		call pgband(2,0,x1,y1,x,y,ch1)
		if(ch1 == ch) then
			x2 = x; y2 = y
		else
			x1=xmin; x2=xmax; y1=ymin; y2=ymax
		endif
	case ('D')                                                       !  zoom in x-direction, middle mouse 2 times
		call pgband(4,0,x1,y1,x,y,ch1)
		if(ch1 == ch) then
			x2 = x; y1 = ymin; y2 = ymax
		else
			x1=xmin; x2=xmax; y1=ymin; y2=ymax			
		endif
	case ('X')                                                       !  zoom in y-direction, right mouse 2 times
		call pgband(3,0,x1,y1,x,y,ch1)
		if(ch1 == ch) then
			x1 = xmin; x2 = xmax; y2 = y
		else
			x1=xmin; x2=xmax; y1=ymin; y2=ymax			
		endif
	case ('z')                                                       !  zoom in
		xc = 0.5*(xmin+xmax); yc = 0.5*(ymin+ymax)	
		hx = xmax-xmin; hy = ymax-ymin
		x1 = xc-0.75*hx; x2 = xc+0.75*hx
		y1 = yc-0.75*hy; y2 = yc+0.75*hy
	case ('Z')                                                       !  zoom out
		xc = 0.5*(xmin+xmax); yc = 0.5*(ymin+ymax)	
		hx = xmax-xmin; hy = ymax-ymin
		x1 = xc-0.25*hx; x2 = xc+0.25*hx
		y1 = yc-0.25*hy; y2 = yc+0.25*hy
	case ('l')                                                       !  scroll left
		hx = xmax-xmin
		x1 = xmin - 0.25*hx; y1 = ymin 
		x2 = xmax - 0.25*hx; y2 = ymax
	case ('r')                                                       !  scroll right
		hx = xmax-xmin
		x1 = xmin + 0.25*hx; y1 = ymin
		x2 = xmax + 0.25*hx; y2 = ymax
	case ('u')                                                       !  scroll up
		hy = ymax-ymin
		y1 = ymin + 0.25*hy; x1 = xmin
		y2 = ymax + 0.25*hy; x2 = xmax
	case ('d')                                                       !  scroll down
		hy = ymax-ymin
		y1 = ymin - 0.25*hy; x1 = xmin
		y2 = ymax - 0.25*hy; x2 = xmax
	case default; x1=xmin; x2=xmax; y1=ymin; y2=ymax                 !  stay as is
	end select
	end subroutine zoomBindingsPgPlotSelect
