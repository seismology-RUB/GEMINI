!------------------------------------------------------------------------
!  Interactive browsing through modes:
!
!  f+n:        move on n frequencies
!  f-n:        go back n frequencies
!  b+n:        move on n branches
!  b-n:        go back n branches
!  c+n:        move on n DSV components
!  c-n:        go back n DSV components
!  p:          print values of current eigenfunction to terminal
!  i:          choose plot region shown by zooming and scrolling
!              using standard pgPlotBindings (A,D,X,l,r,u,d,z,Z)
!  x:          exit program
!------------------------------------------------------------------------
 program plotFlmodes
	use commandLine
	use flmodeHeader
	use flmodeBody
	use pgPlotXY
	implicit none
!
	type (cmdLine) :: cl
	type (flmode_header) :: mhd
	type (flmode_body) :: meig
	type (pgplot_xy) :: pleig
	integer :: nf,nf1,nf2,pgid,jrec,nrt,kf,frec,fjump,kbr,bjump,keig,eigjump,i
	integer :: pgopen
	real :: df,xref,yref
	real, dimension(:), pointer :: eifun
	real, dimension(:), allocatable :: z
	double precision :: rearth
	double precision, dimension(:), pointer :: rnod
	character (len=132) :: modesfile,eigbase,text
	character :: typ*1,ch*1,ch1*1,ch2*1,cmd*2
	external printhelp, pgopen
!
	data typ/'S'/
!-----------------------------------------------------------------------------------
	call new(cl,2,1,'h t','0 1',printhelp)
	eigbase = clManarg(cl,1)
	if(clOptset(cl,2)) typ = trim(clOptarg(cl,2))
	call dealloc(cl)
!
!  open plot window
!
	pgid = pgopen('/XSERVE')
!
!  name of modes file according to mode type
!
	if(typ == 'S') then
		modesfile = trim(eigbase)//'.sph'
	else if(typ == 'T') then
		modesfile = trim(eigbase)//'.tor'
	else
		print *,'Invalid mode type'; stop
	endif
!
!  open modes file
!
	open(1, file=modesfile, status='old',form='unformatted', access = 'direct', recl=24)
	read(1,rec=1) nf1,nf2,df; read(1,rec=2) rearth
	nf = nf2-nf1+1
	print *,'Modentyp: ',typ
	print *,'Frequency index range: ',nf1,nf2
	print *,'Lowest  frequency: ',(nf1-1)*df
	print *,'Highest frequency: ',(nf2-1)*df
!
!  get mode data for first frequency (kf=1) and eigenfunction U or W (keig=1)
!
	kf=1
	keig=1
	kbr=1
	call dealloc(meig)       ! Hack to avoid strange segmentation fault. meig%yeig is claimed to be associated.
	call newFrequency(1,typ,nf,kf,0,kbr,mhd,meig,jrec,nrt)
!
!  start loop for displaying various modes
!
	ch = ' '
	do while (ch /= 'x')
		eifun => getEigfunFlmodeBody(meig,keig)
		rnod => getNodesFlmodeHeader(mhd)
		allocate(z(size(rnod)))
		z = sngl(rnod-rearth)
		write(text,'(a,i3,a,i3,a,i3,a,i1)') 'Frequency:',kf+nf1-1,'; Branch:',kbr,' of',nrt,'; DSV:',keig
		call new(pleig,eifun,z,ylab = 'Depth',title = trim(text))
		call display(pleig,pgid)
		call pgcurs(xref,yref,ch)
		if(ch == 'x') exit
		if(ch == 'f') then
			call pgcurs(xref,yref,ch1); call pgcurs(xref,yref,ch2); cmd = ch1//ch2
			read(cmd,'(i2)') fjump
			call newFrequency(1,typ,nf,kf,fjump,kbr,mhd,meig,jrec,nrt)
		else if(ch == 'b') then
			call pgcurs(xref,yref,ch1);	call pgcurs(xref,yref,ch2); cmd = ch1//ch2
			read(cmd,'(i2)') bjump
			call newOvertone(1,jrec,nrt,kbr,bjump,mhd,meig)
		else if(ch == 'c') then
			call pgcurs(xref,yref,ch1);	call pgcurs(xref,yref,ch2); cmd = ch1//ch2
			read(cmd,'(i2)') eigjump
			call newEigenfunction(mhd,meig,keig,eigjump)
		else if(ch == 'p') then
			do i=1,size(eifun)
				write(6,'(f10.2, e15.4)') z(i),eifun(i)
			enddo
		else if(ch == 'i') then
			call interactivePgPlotXY(pleig,pgid)       !  allow interactive modification of plot
		endif
		call dealloc(pleig)
		deallocate(z)
	enddo
	call dealloc(mhd); call dealloc(meig)
!
 contains
!-----------------------------------------------------------------------------------
!  read in first overtone of new frequency
!
!  lu:       file unit from where mode data are read (input)
!  typ:      mode type (input)
!  nf:       number of frequencies (input) 
!  kf  :     frequency index
!  fjump:    number of frequencies we jump
!  kbr:      branch index (inout)
!  mhd:      flmode header object read out (output)
!  meig:     flmode body object read out (output)
!  jrec:     number of first record after just read in overtone (output)
!  nrt:      number of overtones for this frequency (output)
!
	subroutine newFrequency(lu,typ,nf,kf,fjump,kbr,mhd,meig,jrec,nrt)
	type (flmode_header) :: mhd
	type (flmode_body) :: meig
	integer :: lu,nf,jrec,nrt,frec,kf,fjump,jf,kbr
	character (len=1) :: typ
!
!  chack valid value of frec, if unvalid give warning and return without reading
!
	if(kf+fjump < 1 .or. kf+fjump > nf) then
		print *,'Warning: Mode for this frequency is not available. Keep previous mode'
		return
	else
		kf = kf+fjump
		frec = kf+2
		call dealloc(meig); call dealloc(mhd)
	endif
!
 10	read(lu,rec=frec) jf,jrec   ! frequency index and starting record number of new mode block
	read(lu,rec=jrec) nrt
!
	if(nrt > 0) then
		jrec=jrec+1
		call new(mhd,lu,typ,jrec)
		if(kbr > nrt) then
			print *,'Desired overtone not available at this frequency. I take kbr = 1'
			kbr = 1
		endif
		jrec = jrec + (kbr-1)*getNrecFlmodeHeader(mhd)
		call new(meig,lu,mhd,jrec)
	else
		print *,'No modes exist for the selected frequency'
		print *,'I move to the next one'
		kf=kf+1; frec=frec+1
		if(kf > nf) then
			print *,'Beyond highest frequency. I stop now'
			stop
		else
			goto 10
		endif
	endif 
	end subroutine newFrequency		
!------------------------------------------------------------------------------------
!  read in an overtone for same frequency
!
!  lu:      file unit (input)
!  jrec:    number of first record of next overtone (inout)
!  kbr:     overtone index (inout)
!  nrt:     total number of overtones (input)
!  bjump:   how many overtones we jump ahead (input)
!  mhd:     flmodeHeader object (input)
!  meig:    new overtone data (output)
!
	subroutine newOvertone(lu,jrec,nrt,kbr,bjump,mhd,meig)
	integer :: lu,jrec,nrt,kbr,bjump
	type (flmode_body) :: meig
	type (flmode_header) :: mhd
!
	if(kbr+bjump < 1 .or. kbr+bjump > nrt) then
		print *,'Warning: desired overtone not available. Keep previous one!'
		return
	else
		call dealloc(meig)
		kbr = kbr+bjump
		jrec = jrec + (bjump-1)*getNrecFlmodeHeader(mhd)
	endif
	call new(meig,lu,mhd,jrec)
	end subroutine newOvertone
!-----------------------------------------------------------------------------------
!  jump to new eigenfunction
!
!  mhd:      flmode_header object (input)
!  meig:     flmode_body object (input)
!  keig:     current eigenfunction index
!  eigjump:  how many eigenfunctions we jump
!
	subroutine newEigenfunction(mhd,meig,keig,eigjump)
	integer :: keig,eigjump
	type (flmode_header) :: mhd
	type (flmode_body) :: meig
	integer :: eiglim
!
	if(mhd%typ == 'S') then; eiglim = 6; else; eiglim = 3; endif
	if(keig+eigjump < 1 .or. keig+eigjump > eiglim) then
		print *,'Warning: desired eigenfunction not available. Keep previous one!'
		return
	else
		keig=keig+eigjump
	endif
	end subroutine newEigenfunction
!
 end program	
!------------------------------------------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   PLOTFLMODES'
	print *,''
	print *,'Usage: plotFlmodes [-h] [-t] eigbase'
	print *,''
	print *,'Interactive display of eigenfunctions. Data are read from modes file'
	print *,'which is either <eigbase.sph> or <eigbase.tor>.'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-t    : type of mode, either S or T (def = S)'
	print '(50(1h-))'
	return
	end
