 module pgPlotXY
 	implicit none
	private :: createPgPlotXY, displayPgPlotXY, deallocpgPlotXY
	interface new; module procedure createPgPlotXY; end interface
	interface display; module procedure displayPgPlotXY; end interface
	interface overlay; module procedure overlayPgPlotXY; end interface
	interface dealloc; module procedure deallocPgPlotXY; end interface
!
	type pgplot_xy
		private
		real, dimension(:), pointer :: x,y
		real :: xmin,xmax,ymin,ymax,ch
		integer :: ls,lw,ci
		integer :: eqar,axmode,symflag
		character (len=132) :: xlab,ylab,title,plotfile
	end type pgplot_xy
 contains
!--------------------------------------------------------------------
!  an object used to plot an xy-graph using pgline or pgpnts
!
!  x:	  	    x-values
!  y:             y-values
!  x1,xn,y1,yn:   extremal values (def = calculate from data)
!  just:          equal area or not (def = 0)
!  iaxis:         axis mode (see pgenv) (def = 0)
!  symflag:       symbol (1) or line (0) (def = 0)
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.5)
!  xlab:	    label at x-axis (def = empty)
!  ylab:	    label at y-axis (def = empty)
!  title: 	    title on top (def = empty)
!  plotfile:      output file for hardcopy (def = landscape.ps)
!
	subroutine createPgPlotXY(this,x,y,x1,xn,y1,yn,just,iaxis,symflag,lw,ci,ls,ch,xlab,ylab,title,plotfile)
	type (pgplot_xy) :: this
	real, dimension(:), target :: x,y
	real, optional :: x1,xn,y1,yn
	character (len=*), optional :: xlab,ylab,title,plotfile
	integer, optional :: just,iaxis,ls,lw,ci,symflag
	real, optional :: ch
!
	this%x => x; this%y => y
	if (present(x1)) then; this%xmin = x1; else; this%xmin = minval(x); endif
	if (present(xn)) then; this%xmax = xn; else; this%xmax = maxval(x); endif
	if (present(y1)) then; this%ymin = y1; else; this%ymin = minval(y); endif
	if (present(yn)) then; this%ymax = yn; else; this%ymax = maxval(y); endif
	if (present(just)) then; this%eqar = just; else; this%eqar = 0; endif
	if (present(iaxis)) then; this%axmode = iaxis; else; this%axmode = 0; endif
	if (present(symflag)) then; this%symflag = symflag; else; this%symflag = 0; endif
	if (present(lw)) then; this%lw = lw; else; this%lw=2; endif
	if (present(ci)) then; this%ci = ci; else; this%ci=1; endif
	if (present(ls)) then; this%ls = ls; else; this%ls=1; endif
	if (present(ch)) then; this%ch = ch; else; this%ch=1.5; endif
	if (present(xlab)) then; this%xlab = xlab; else; this%xlab=''; endif
	if (present(ylab)) then; this%ylab = ylab; else; this%ylab=''; endif
	if (present(title)) then; this%title = title; else; this%title=' '; endif
	if (present(plotfile)) then; this%plotfile = plotfile; else; this%plotfile='landscape.ps'; endif
!
	end subroutine createPgPlotXY
!--------------------------------------------------------------------
!  free memory
!
	subroutine deallocpgPlotXY(this)
	type (pgplot_xy) :: this
	nullify(this%x, this%y)
	end subroutine deallocpgPlotXY
!--------------------------------------------------------------------
!  plot graph
!
	subroutine displayPgPlotXY(this,pgid)
	type (pgplot_xy) :: this
	integer :: pgid
	real :: xmin,xmax,ymin,ymax
	real :: xr,yr,xm,ym
!
	call pgslct(pgid)
	call pgslw(1);call pgsci(1);call pgscf(2); call pgsch(1.25)
	call pgask(.false.)
	call pgenv(this%xmin,this%xmax,this%ymin,this%ymax,this%eqar,this%axmode)     ! set up viewport and axes
	call pgslw(3); call pgsch(this%ch)
	call pglab(this%xlab,this%ylab,this%title)
	call pgslw(this%lw)
	call pgsci(this%ci)
	if(this%symflag == 0) then
		call pgline(size(this%x),this%x,this%y)
	else
		call pgpt(size(this%x),this%x,this%y,this%symflag)
	endif
	end subroutine displayPgPlotXY
!----------------------------------------------------------------------
!  overlay plot
!
	subroutine overlayPgPlotXY(this,pgid)
	type (pgplot_xy) :: this
	integer :: pgid
!
	call pgslct(pgid)
	call pgslw(this%lw)
	call pgsci(this%ci)
	if(this%symflag == 0) then
		call pgline(size(this%x),this%x,this%y)
	else
		call pgpt(size(this%x),this%x,this%y,this%symflag)
	endif
	end subroutine overlayPgPlotXY
!----------------------------------------------------------------------
!  interactive choice of plot region by zooming in various ways
!
	subroutine interactivePgPlotXY(this,pgid)
	type (pgplot_xy) :: this
	real :: x1,y1,x2,y2
	character (len=1) :: ch
	integer :: pgid
	do
		call pgcurs(x1,y1,ch)
		if (ch == 'x') exit
		call zoomBindingsPgPlotSelect(pgid,x1,y1,ch,x2,y2,this%xmin,this%xmax,this%ymin,this%ymax)
		call setExtremaPgPlotXY(this,x1,x2,y1,y2)
		call display(this,pgid)
		call hardcopyBindingPgPlotXY(this,pgid,ch)
	enddo
	end subroutine interactivePgPlotXY
!-----------------------------------------------------------------------
	subroutine hardcopyBindingPgPlotXY(this,pgid,ch)
	type (pgplot_xy) :: this
	integer :: pgopen,idev,pgid
	character (len=1) :: ch
	if(ch /= 'h') return
	idev = pgopen(trim(this%plotfile)//'/CPS')
	call display(this,idev)
	call pgclos
	print *,'PgPlotXY saved to ',trim(this%plotfile)
	call pgslct(pgid)
	end subroutine hardcopyBindingPgPlotXY
!---------------------------------------------------------------------
!  set extrema
!
	subroutine setExtremaPgPlotXY(this,x1,x2,y1,y2)
	type (pgplot_xy) :: this
	real :: x1,x2,y1,y2
	this%xmin=x1; this%xmax=x2; this%ymin=y1; this%ymax=y2
	end subroutine setExtremaPgPlotXY
!
 end module pgPlotXY
