!-----------------------------------------------------------------------------
!> \brief Plot frequency wavenumber spectra
!-----------------------------------------------------------------------------
 program plotGreenFKSpectra
	use mathConstants
	use commandLine
	use greenFrequencyWavenumber
	use errorMessage
	use pgPlotImage
	use pgPlotWindow
	implicit none
	type (cmdLine) :: cl
	type (error_message) :: errmsg
	type (green_frequency_wavenumber) :: gfk
	complex, dimension(:), pointer :: fk
	real, dimension(:,:), allocatable :: fkplot
	real :: zs,bg,fkmax,fmin,fmax,wmax
	integer :: cnt,nwe,if,nwn,isp,jsp
	logical :: vsflag,parflag
	character (len=132) :: dsvbasename,string
	character (len=2) :: cdsv
	type (pgplot_image) :: pgfk
	type (pgplot_window) :: pgwin
	real :: x1,x2,y1,y2
	character (len=1) :: ch,cjsp
	external printhelp
!-----------------------------------------------------------------
	data vsflag/.false./parflag/.false./zs/10000./bg/-4./
!-----------------------------------------------------------------
	call new(cl,5,1,'h z bg vs p','0 1 1 0 0',printhelp)
	dsvbasename = clManarg(cl,1)
	if(clOptset(cl,'z')) then; string=clOptarg(cl,'z'); read(string,*) zs; endif
	if(clOptset(cl,'bg')) then; string=clOptarg(cl,'bg'); read(string,*) bg; endif
	if(clOptset(cl,'vs')) then; vsflag = .true.; endif
	if(clOptset(cl,'p')) then; parflag = .true.; endif
	call dealloc(cl)
!-------------------------------------------------------
!  open frequency-wavenumber spectrum file
!  and read in group tree and header info
!
	if (.not. parflag) then
		errmsg = readHeaderGreenFrequencyWavenumber(gfk,1,2,dsvbasename,vsflag)
	else
		errmsg = readHeaderParallelGreenFrequencyWavenumber(gfk,7,dsvbasename,vsflag)
	endif
	if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  read fk-spectra for depth zs
!
	if (.not. parflag) then
		errmsg = readBodyGreenFrequencyWavenumber(gfk,zs,vsflag)
	else
		errmsg = readParallelBodyGreenFrequencyWavenumber(gfk,zs,vsflag)
	endif
	if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  start plotting with DSV-index isp = 1 and jump-index jsp = 1
!  should always be available
!
	fmin = (.nfa.gfk-1)*.df.gfk
	fmax = (.nfb.gfk-1)*.df.gfk
	wmax = (.nwnmax.gfk-1)*.dwn.gfk
	allocate(fkplot(.nf.gfk,.nwnmax.gfk))
	call new(pgwin)
!
!  print out which spectra are available
!
	call tableContentsGreenFrequencyWavenumber(gfk)
	print *,'Available interactive options: A D X z Z l r u d o s[isp][jsp] b[n]'
!
	isp = 1; jsp = 1
 1	if (isp > 9 .or. jsp > 4 .or. isp < 1 .or. jsp < 1) then
		print *,'Invalid values for isp,jsp. Choose from 1 <= isp <= 9 and 1 <= jsp <= 4'
		goto 3
	endif
	cdsv = getNameDsvFrequencyWavenumberMapping(isp)
	cnt = getLinearSpectralCountGreenFrequencyWavenumber(gfk,isp,jsp)
	if (cnt == 0) then
		print *,'Spectrum DSV = ',cdsv,' Jump = ',jsp,' is not available. Choose a different one'
		goto 3
	endif
 	fk => getSelectedSpectrumGreenFrequencyWavenumber(gfk,cnt)
	fkmax = maxval(cabs(fk))
	write(cjsp,'(i1)') jsp
!
!  unwrap fk-spectrum into matrix with x->f and y->k
!
 4	fkplot = bg
	nwe = 0
	do if = .nfa.gfk,.nfb.gfk
		nwn = getNwnSelectedGreenFrequencyWavenumber(gfk,if)
		fkplot(if-.nfa.gfk+1,1:nwn) = max(log10(cabs(fk(nwe+1:nwe+nwn)/fkmax)),bg)
		nwe = nwe+nwn;
	enddo
!
!  prepare plot image
!
	call new(pgfk,fkplot,x1 = fmin,xn = fmax,y1 = 0.0,yn = wmax,xlab = 'Frequency',ylab = 'Wavenumber', &
		  title = 'FK-Spectrum for DSV = '//cdsv//' and jump = '//cjsp)
!
!  plot
!
 2	call display(pgfk)
 3	call pgcurs(x1,y1,ch)
	select case (ch)
	case ('A','D','X','z','Z','l','r','u','d')
		call zoomBindingsPgPlotSelect(x1,y1,ch,x2,y2,.xmin.pgfk,.xmax.pgfk,.ymin.pgfk,.ymax.pgfk)
		call setExtremaPgPlotImage(pgfk,x1,x2,y1,y2)
		goto 2
	case ('o'); call setOriginalExtremaPgPlotImage(pgfk); goto 2
	case ('s'); call pgcurs(x1,y1,ch); read(ch,'(i1)') isp; call pgcurs(x1,y1,ch); read(ch,'(i1)') jsp; call dealloc(pgfk); goto 1
	case ('b'); call pgcurs(x1,y1,ch); read(ch,*) bg; bg = -bg; call dealloc(pgfk); goto 4
	end select
!
!  clean up
!
	deallocate(fkplot)
	call dealloc(pgfk)
	call dealloc(gfk)
!	
 end program
!-----------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   PLOTGreenFKSpectra'
	print *,''
	print *,'Usage: plotGreenFKSpectra [-h] [-z] [-bg] [-vs] [-p] dsvbasename'
	print *,''
	print *,'Plot frequency wavenumber spectra read from files with'
	Print *,'names starting with <dsvbasename>.'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-z    : depth of fk-spectra in meters '
	print *,'-bg   : log background value for spectrum (def = -4)'
	print *,'-vs   : verbosity flag for stream access input (def = false)'
	print *,'-p    : FK-spectra were calculated with a parallel run (def = .false.)'
	print '(50(1h-))'
	return
	end
