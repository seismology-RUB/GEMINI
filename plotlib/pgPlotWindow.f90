!-----------------------------------------------------------------------
!> \brief Module for a pgplot window
!
 module pgPlotWindow
	implicit none
	interface new; module procedure createPgPlotWindow; end interface
	interface dealloc; module procedure closePgPlotWindow; end interface
	interface select_win; module procedure selectPgPlotWindow; end interface
	interface select_panel; module procedure selectPanelPgPlotWindow; end interface
	type pgplot_window
		private
		integer :: pgid
		character (len=132) :: device, plotfile
		integer :: nxsub,nysub
		real :: width,aspect
		logical :: hardcopy             !< .true.: plotting to postscript
	end type
!
 contains
!-----------------------------------------------------------------
!> \brief  Create pgplot window
!
!> \param  nxsub, nysub	subdivision of window into panels (optional)
!> \param  plotfile		name of file for postscript output (optional)
!!				plots to postscript if this argument is set
!> \param  psdev		postscript device name (CPS or VCPS) (optional)
!> \param  seldev		flag if menu for selecting device is desired
!
	subroutine createPgPlotWindow(this,nxsub,nysub,width,aspect,plotfile,psdev,seldev)
	type (pgplot_window) :: this
	integer :: pgopen,selidx
	character (len=*), optional :: plotfile,psdev
	integer, optional :: nxsub,nysub
	logical, optional :: seldev
	real, optional :: width,aspect
!
	if (present(seldev)) then
		selidx = selectDevicePgPlotWindow()
		if (selidx == 1) then; this%device = '/XSERVE'; this%hardcopy = .false.; endif
		if (selidx == 2) then; this%device = 'landscape.ps/CPS'; this%hardcopy = .true.; endif
		if (selidx == 3) then; this%device = 'portrait.ps/VCPS'; this%hardcopy = .true.; endif
	else
		if (.not. present(plotfile)) then
			this%device = '/XSERVE'
			this%hardcopy = .false.
		else
			if (present(psdev)) then
				this%device = trim(plotfile)//'/'//trim(psdev)
			else
				this%device = trim(plotfile)//'/CPS'
			endif
			this%hardcopy = .true.
		endif
	endif
	this%pgid = pgopen(trim(this%device))
	if (this%pgid < 0) then
		print *,'createPgPlotWindow: Device ',this%device,' could not be opened'
		stop
	endif
	if (present(aspect)) then; this%aspect = aspect; else; this%aspect = 0.75; endif	
	if (present(width)) then; this%width = width; else; this%width = 0.; endif
	if (present(nxsub)) then; this%nxsub = nxsub; else; this%nxsub = 1; endif
	if (present(nysub)) then; this%nysub = nysub; else; this%nysub = 1; endif
	call pgpap(this%width,this%aspect)
	if (this%nxsub > 1 .or. this%nysub > 1) call pgsubp(this%nxsub,this%nysub)
	end subroutine createPgPlotWindow
!-------------------------------------------------------------------
!> \brief  close plot window
!
	subroutine closePgPlotWindow(this)
	type (pgplot_window) :: this
	call pgslct(this%pgid)
	call pgclos
	end subroutine closePgPlotWindow
!--------------------------------------------------------------------
!> \brief  dialog to select device (internal use)
!
	function selectDevicePgPlotWindow() result(seldev)
	integer :: pgid,pgopen,seldev
	character (len=1) :: ch
!
	pgid = pgopen('/XSERVE')
	call pgsch(1.0)
	call pgsvp(0.,1.,0.,1.); call pgswin(0.,1.,0.,1.)
	call pgtext(0.1,0.9,'Choose among the following output devices: ')
	call pgtext(0.15,0.8,'Screen:');                call pgtext(0.5,0.8,' -------->      1')
	call pgtext(0.15,0.7,'Postscript Landscape:');  call pgtext(0.5,0.7,' -------->      2')
	call pgtext(0.15,0.6,'Postscript Portrait:');   call pgtext(0.5,0.6,' -------->      3')
	call pgcurs(0.5,0.5,ch); read(ch,'(i1)') seldev
	call pgclos
	end function selectDevicePgPlotWindow
!--------------------------------------------------------------------
!> \brief  select window
!
	subroutine selectPgPlotWindow(this)
	type (pgplot_window), intent(in) :: this
	call pgslct(this%pgid)
	call pgpap(this%width,this%aspect)
	end subroutine selectPgPlotWindow
!--------------------------------------------------------------------
!> \brief  select subpanel
!
	subroutine selectPanelPgPlotWindow(this,ix,iy)
	type (pgplot_window), intent(in) :: this
	integer :: ix,iy
	call pgslct(this%pgid)
	call pgpanl(ix,iy)
	end subroutine selectPanelPgPlotWindow
!
 end module
