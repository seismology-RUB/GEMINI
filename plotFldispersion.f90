 program plotFldispersion
	use commandLine
	use flmodeHeader
	use flmodeBody
	use pgPlotXY
	implicit none
!
	type (cmdLine) :: cl
	type (flmode_header) :: mhd
	type (flmode_body) :: meig
	type (pgplot_xy) :: pldisp
	integer :: nf1,nf2,pgid,jrec,nrt,kf,nbrtotal,n,jf,nbr,kd,nov
	integer :: pgopen
	real :: df,xref,yref,f,p,cg,q
	real, dimension(:,:), target, allocatable :: dispersion_matrix     ! columns: f,p,c,g,q
	real, dimension(:), pointer :: frequency, dispersion
	character (len=132) :: modesfile,eigbase,text
	character :: typ*1,ch*1
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
	read(1,rec=1) nf1,nf2,df
	print *,'Modentyp: ',typ
	print *,'Frequency index range: ',nf1,nf2
	print *,'Lowest  frequency: ',(nf1-1)*df
	print *,'Highest frequency: ',(nf2-1)*df
!
!  loop through frequencies and calculate total number of branches
!
	nbrtotal = 0
	do jf=nf1,nf2
		read(1,rec=jf-nf1+1+2) kf,jrec
		read(1,rec=jrec) nrt
		nbrtotal=nbrtotal+nrt
	enddo
	print *,'Total number of branches: ',nbrtotal
!
!  allocate space for dispersion_matrix
!
	allocate(dispersion_matrix(nbrtotal,5))
!
!  loop again over frequencies to read dispersion data
!  into dispersion_matrix
!
	nbr = 0
	do jf=nf1,nf2
		read(1,rec=jf-nf1+1+2) kf,jrec
		read(1,rec=jrec) nrt
		jrec=jrec+1
		call new(mhd,1,typ,jrec)
		do n=1,nrt
			read(1,rec=jrec) f,nov,p,cg,q
			dispersion_matrix(n+nbr,:) = (/f,p,1./p,cg,q/)
			jrec = jrec + getNrecFlmodeHeader(mhd)
		enddo
		call dealloc(mhd)
		nbr=nbr+nrt
	enddo
!
!  loop for displaying p,c,g,or q
!
	ch = ' '
	kd = 2
	text = 'Slowness'
	do while (ch /= 'x')
		frequency => dispersion_matrix(:,1)
		dispersion => dispersion_matrix(:,kd)
		call new(pldisp,frequency,dispersion,symflag = -4,ch = 1.0,xlab = 'Frequency',ylab = trim(text),&
	&	         title = 'Dispersion properties')
		call display(pldisp,pgid)
		call pgcurs(xref,yref,ch)
		select case (ch)
		case ('p'); kd = 2; text = 'Slowness'
		case ('c'); kd = 3; text = 'Phase velocity'
		case ('g'); kd = 4; text = 'Group velocity'
		case ('q'); kd = 5; text = 'Quality factor'
		case ('i'); call interactivePgPlotXY(pldisp,pgid)
		end select
		call dealloc(pldisp)
	enddo
!
	deallocate(dispersion_matrix)
!
 end program plotFldispersion		
!------------------------------------------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   PLOTFLDISPERSION'
	print *,''
	print *,'Usage: plotFldispersion [-h] [-t] eigbase'
	print *,''
	print *,'Interactive display of dispersion. Data are read from modes file'
	print *,'which is either <eigbase.sph> or <eigbase.tor>.'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-t    : type of mode, either S or T (def = S)'
	print '(50(1h-))'
	return
	end
