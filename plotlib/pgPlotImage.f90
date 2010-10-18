 module pgPlotImage
	use pgColorScale
	use plotUtilities
 	implicit none
	interface new; module procedure createPgPlotImage; end interface
	interface display; module procedure displayPgPlotImage;	end interface
	interface dealloc; module procedure deallocPgPlotImage;	end interface
	interface operator (.xmin.); module procedure getXminPgPlotImage; end interface
	interface operator (.xmax.); module procedure getXmaxPgPlotImage; end interface
	interface operator (.ymin.); module procedure getYminPgPlotImage; end interface
	interface operator (.ymax.); module procedure getYmaxPgPlotImage; end interface
	type pgplot_image
		real, dimension(:,:), pointer:: a
		real :: xmin,xmax,ymin,ymax
		real :: xomin,xomax,yomin,yomax
		real :: a1,a2,c1,c2,tr(6)
		integer :: nx,ny,nc,eqar
		character (len=10) :: axis
		character (len=100) :: xlab,ylab,title
		type (pg_colscale) :: cs
	end type pgplot_image
 contains
!--------------------------------------------------------------------
!  an object used to plot a shaded image using pgimag
!  optionally with contour lines
!
!  a:		2D-array
!  x1:		x-value belonging to a(1,:)
!  x2:		x-value belonging to a(n,:)
!  y1:		y-value belonging to a(:,1)
!  yn:		y-value belonging to a(:,n)
!  just:	1: axis equal; 0 : axes unequal
!  iaxis:     see pgenv; reommended is 1
!  nc:		number of contour levels
!  xlab:	label at x-axis
!  ylab:	label at y-axis
!  title: 	title on top
!
	subroutine createPgPlotImage(this,a,x1,xn,y1,yn,just,iaxis,nc,xlab,ylab,title,colscale)
	type (pgplot_image) :: this
	real, dimension(:,:), target :: a
	real, optional :: x1,xn,y1,yn
	integer, optional :: nc,just,iaxis
	character (len=*), optional :: xlab,ylab,title,colscale
	real :: xa,xb,ya,yb,dx,dy
!
	this%a => a;
	this%nx = size(a,1); this%ny = size(a,2)
	if (present(x1)) then; xa=x1; else; xa = 1.; endif
	if (present(xn)) then; xb=xn; else; xb = this%nx; endif
	if (present(y1)) then; ya=y1; else; ya = 1.; endif
	if (present(yn)) then; yb=yn; else; yb = this%ny; endif
	if (present(just)) then; this%eqar = just; else; this%eqar = 0; endif
	if (present(iaxis)) then; this%axis = convertAxisToBox(iaxis); else; this%axis = 'BCNST'; endif
	this%xmin = xa; this%xmax = xb; this%ymin=ya; this%ymax=yb
	if (xa > xb) then; this%xmin=xb; this%xmax=xa; endif
	if (ya > yb) then; this%ymin=yb; this%ymax=ya; endif
	this%a1 = minval(a); this%a2 = maxval(a)
	this%c1 = this%a1; this%c2 = this%a2
	dx = (xb-xa)/(this%nx-1); dy = (yb-ya)/(this%ny-1)
	this%tr = (/xa-dx,dx,0.,ya-dy,0.,dy/)
	this%xomax = this%xmax; this%xomin = this%xmin; this%yomin = this%ymin; this%yomax = this%ymax
	if (present(xlab)) then; this%xlab = xlab; else; this%xlab=' '; endif
	if (present(ylab)) then; this%ylab = ylab; else; this%ylab=' '; endif
	if (present(title)) then; this%title = title; else; this%title=' '; endif
	if (present(nc)) then; this%nc = nc; else; this%nc = 0; endif
	if (present(colscale)) then; call new(this%cs,trim(colscale)); else; call new(this%cs,'rwb'); endif
!
	end subroutine createPgPlotImage
!--------------------------------------------------------------------------------
!  free memory
!
	subroutine deallocpgPlotImage(this)
	type (pgplot_image) :: this
	nullify(this%a)
	call dealloc(this%cs)
	end subroutine deallocpgPlotImage
!--------------------------------------------------------------------------------
!  plot image
!
	subroutine displayPgPlotImage(this)
	type (pgplot_image) :: this
	real, dimension(:), allocatable :: cont
	integer :: i
!
	if(this%nx < 2 .or. this%ny < 2) then
		print *,'Array is essentially 1-dimensional!'
		print *,'nx = ',this%nx,' ny = ',this%ny
		print *,'I do nothing !!!'
		return
	endif
!
	call pgsci(1)
	call pgeras
	call pgsvp(0.1,0.9,0.2,0.9)                     ! leave some space at bottom for wedge
	if (this%eqar == 0) then
		call pgswin(this%xmin,this%xmax,this%ymin,this%ymax)
	else
		call pgwnad(this%xmin,this%xmax,this%ymin,this%ymax)
	endif
	call pgbox(trim(this%axis),0.0,0,trim(this%axis),0.0,0)
!
	call apply(this%cs)                                                               ! apply color scale
	call pgimag(this%a,this%nx,this%ny,1,this%nx,1,this%ny,this%a1,this%a2,this%tr)          ! plot image
	call pgwedg('BI', 5., 3., this%a1, this%a2, '')
!
	if (this%nc > 0) then                            ! plot contour lines if desired
		allocate(cont(this%nc))
		do i=1,this%nc
			cont(i) = this%c1+(this%c2-this%c1)/(this%nc-1)*(i-1)
		enddo
		call pgcont(this%a,this%nx,this%ny,1,this%nx,1,this%ny,cont,this%nc,this%tr)
		deallocate(cont)
	endif
	call pgmtxt ('B', 3., 0.5, 0.5, this%xlab)
	call pgmtxt ('L', 3., 0.5, 0.5, this%ylab)
	call pgmtxt ('T', 3., 0.5, 0.5, this%title)
!
	end subroutine displayPgPlotImage
!---------------------------------------------------------------------------
!  change lower contour bound by a fraction of the total range
!
	subroutine changeClbPgPlotImage(this,fac)
	type (pgplot_image) :: this
	real :: fac
	this%c1 = this%c1 + fac*(this%c2-this%c1)
	end subroutine changeClbPgPlotImage
!---------------------------------------------------------------------------
!  change upper contour bound by a fraction of the total range
!
	subroutine changeCubPgPlotImage(this,fac)
	type (pgplot_image) :: this
	real :: fac
	this%c2 = this%c2 + fac*(this%c2-this%c1)
	end subroutine changeCubPgPlotImage
!---------------------------------------------------------------------------
!  change lower amplitude bound by a fraction of the total range
!
	subroutine changeAlbPgPlotImage(this,fac)
	type (pgplot_image) :: this
	real :: fac
	this%a1 = this%a1 + fac*(this%a2-this%a1)
	end subroutine changeAlbPgPlotImage
!---------------------------------------------------------------------------
!  change upper amplitude bound by a fraction of the total range
!
	subroutine changeAubPgPlotImage(this,fac)
	type (pgplot_image) :: this
	real :: fac
	this%a2 = this%a2 + fac*(this%a2-this%a1)
	end subroutine changeAubPgPlotImage
!---------------------------------------------------------------------------
!  change number of contour lines by ncadd
!
	subroutine changeNContourPgPlotImage(this,ncadd)
	type (pgplot_image) :: this
	integer :: ncadd
	this%nc = max(this%nc + ncadd,0)
	end subroutine changeNContourPgPlotImage
!----------------------------------------------------------------------------
!  set new extrema
!
	subroutine setExtremaPgPlotImage(this,x1,x2,y1,y2)
	type (pgplot_image) :: this
	real :: x1,x2,y1,y2
	this%xmin = x1; this%xmax = x2; this%ymin = y1; this%ymax = y2
	end subroutine setExtremaPgPlotImage
!----------------------------------------------------------------------------
!  set original extrema
!
	subroutine setOriginalExtremaPgPlotImage(this)
	type (pgplot_image) :: this
	this%xmin = this%xomin; this%xmax = this%xomax
	this%ymin = this%yomin; this%ymax = this%yomax
	end subroutine setOriginalExtremaPgPlotImage
!---------------------------------------------------------------------------
!  get xmin
!
	real function getXminPgPlotImage(this)
	type (pgplot_image), intent(in) :: this
	getXminPgPlotImage = this%xmin
	end function getXminPgPlotImage
!---------------------------------------------------------------------------
!  get xmax
!
	real function getXmaxPgPlotImage(this)
	type (pgplot_image), intent(in) :: this
	getXmaxPgPlotImage = this%xmax
	end function getXmaxPgPlotImage
!---------------------------------------------------------------------------
!  get ymin
!
	real function getYminPgPlotImage(this)
	type (pgplot_image), intent(in) :: this
	getYminPgPlotImage = this%ymin
	end function getYminPgPlotImage
!---------------------------------------------------------------------------
!  get ymax
!
	real function getYmaxPgPlotImage(this)
	type (pgplot_image), intent(in) :: this
	getYmaxPgPlotImage = this%ymax
	end function getYmaxPgPlotImage
!
 end module pgPlotImage
