!-----------------------------------------------------------------
 program modeseis
	use commandLine
	use mathConstants
	use flmodeHeader
	use flmodeBody
	use sourceReceiver
	use sffHeader
	use sffDatablock
	use frequencyTime
	use timeSeries
	implicit none
!
	type (cmdLine) :: cl
	type (flmode_header) :: mhd
	type (flmode_body) :: meig
	type (sar) :: info
	type (sff_header) :: head
	type (sff_datablock) :: dbl
	real, parameter :: unconv_dis = 1.e-6, unconv_pres = 1.e-9, sigma = 0.
	integer :: nf,nf1,nf2,if,i,j,m,nrtmin,nrtmax,nordh,nordl,nr,nrt,nrt2,n,nstf,nsamp,ncv,nmom
	integer :: ifqtf,iprint,iprec,ifwin,jrec(2)
	real :: epimax,df,fmin,fmax,gv,tcut,fch,fcl,dtp,depsrc,dtstf,dt,dum
	real, dimension(6) :: rm
	double precision :: rearth,rs
	real, dimension(:), allocatable :: stf,urs,ftemp
	real, dimension(6,6) :: basemom
	complex, dimension(:,:,:), allocatable :: zdis
	complex, dimension(:), allocatable :: hflow,hfhigh
	complex :: psi(0:2)
	character (len=132) :: string,sfffile,infofile,eigbase,stfasc
	character :: csys*1,typ(2)*1
	logical last, bmflag
	external printhelp
!-----------------------------------------------------------
!  defaults for input parameters
!
	data nordh/0/sfffile/12hmodeseis.sff/ &
	     fch/10./nordl/0/fcl/100./epimax/100./ifwin/0/ &
	     gv/2./tcut/0./dtp/0.001/ifqtf/0/iprint/0/ &
	     iprec/1/nrtmin/1/nrtmax/1000/nmom/1/bmflag/.false./
!-----------------------------------------------------------
!  basic moment tensors, first index: tensor component, 
!  second index: basis moment index
!
	data basemom/ 1.e20, 1.e20, 1.e20,    0.,     0.,    0., &       !  isotropic part
     	           & -1.e20, 1.e20,    0.,    0.,     0.,    0., &       !  rt-45-double couple
     	           & -1.e20,    0., 1.e20,    0.,     0.,    0., &       !  rf-45-double couple
     	           &     0.,    0.,    0.,    0.,     0., 1.e20, &       !  tf-double couple
     	           &     0.,    0.,    0.,    0.,  1.e20,    0., &       !  rf-double couple
     	           &     0.,    0.,    0., 1.e20,     0.,    0./         !  rt-double couple
!------------------------------------------------------------------------
	call new(cl,12,2,'h o bm hp lp e gvw dt qtf p b c','0 1 0 1 1 1 1 1 1 1 1 1',printhelp)
	infofile = clManarg(cl,1)
	eigbase  = clManarg(cl,2)
	if(clOptset(cl,2))  sfffile = clOptarg(cl,2)
	if(clOptset(cl,3))  then; bmflag = .true.; nmom = 6; endif
	if(clOptset(cl,4))  then; string = clOptarg(cl,4); read(string,*) nordh,fch; endif  
	if(clOptset(cl,5))  then; string = clOptarg(cl,5); read(string,*) nordl,fcl; endif  
	if(clOptset(cl,6))  then; string = clOptarg(cl,6); read(string,*) epimax; endif
	if(clOptset(cl,7))  then; ifwin = 1; string = clOptarg(cl,7); read(string,*) gv,tcut; endif
	if(clOptset(cl,8))  then; string = clOptarg(cl,8); read(string,*) dtp; endif
	if(clOptset(cl,9))  then; stfasc=clOptarg(cl,9); ifqtf=1; endif
	if(clOptset(cl,10)) then; string=clOptarg(cl,10); read(string,*) iprint,iprec; endif
	if(clOptset(cl,11)) then; string = clOptarg(cl,11); read(string,*) nrtmin,nrtmax; endif
!
!  Read info file
!
	call new(info,1,infofile,epimax)
	if(sar_getNetFlag(info) .lt. 2) then
		if(.not.clOptset(cl,12)) then
			print *,'You need either an infofile with component specification'
			print *,'or you need to specify a component on the command line'
			stop
		else
			string = clOptarg(cl,12)
			call sar_setComps(info,string(1:1))
		endif
	endif
	call dealloc(cl)
!---------------------------------------------------------------
!  Prepare loop over frequencies
!---------------------------------------------------------------
	open(1,file=trim(eigbase)//'.sph',form='unformatted',access='direct',recl=24)
	open(2,file=trim(eigbase)//'.tor',form='unformatted',access='direct',recl=24)
	read(1,rec=1) nf1,nf2,df; read(1,rec=2) rearth
	read(2,rec=1) nf1,nf2,df; read(2,rec=2) rearth
	nf=nf2-nf1+1
	fmin=(nf1-1)*df
	fmax=(nf2-1)*df
	print *,'Frequencies: ',fmin,fmax,nf1,nf2,df,rearth
	nrtmin=max(nrtmin,1)
!
!  extract some information from info file
!
	call sar_rad2km(info,sngl(rearth))
	nr = sar_getNumberTraces(info)
	depsrc = sar_getSourceDepth(info)
	call sar_getMomentTensor(info,rm)     ! moment tensor from info file
	csys = sar_CoordinateSystem(info)
	if(csys == 'S') then                  ! depsrc in km for 'S' and im meters for 'C'
		rs = rearth-dble(depsrc)       ! rs always in km
	else
		rs = rearth-1.d-3*dble(depsrc)
	endif
!---------------------------------------------------------------
!  add up modal displacements for each frequency
!---------------------------------------------------------------------
	typ(1) = 'S'
	typ(2) = 'T'
	if (bmflag) then; allocate(zdis(nf2,nr,6)); else; allocate(zdis(nf2,nr,1)); endif
	do if=1,nf2
		zdis(if,:,:)=(0., 0.)
	enddo
!
!  spheroidal contributions: m=1; toroidal contributions: m=2
!
	do if=nf1,nf2
		do m=1,2
			read(m,rec=if-nf1+1+2) i,jrec(m)
			read(m,rec=jrec(m)) nrt; jrec(m)=jrec(m)+1
			if(nrt.eq.0) cycle
			call new(mhd,m,typ(m),jrec(m))
			jrec(m) = jrec(m) + (nrtmin-1)*getNrecFlmodeHeader(mhd)
			nrt2 = min(nrt,nrtmax)
			do n=nrtmin,nrt2
				call new(meig,m,mhd,jrec(m))
				do j=1,nmom
					if (bmflag) rm(1:6) = basemom(1:6,j)
					call excitationFlmodeBody(meig,rm,rs,psi)
					do i=1,nr
						zdis(if,i,j) = zdis(if,i,j) + &
 &					             displacementFlmodeBody(meig,csys,info.comp.i,info.dis.i, &
 &					                                    info.phi.i,info.propdir.i,psi)
					enddo
				enddo
				call dealloc(meig)
			enddo
			call dealloc(mhd)
		enddo
		do i=1,nr
			do j=1,nmom
				if((info.comp.i).eq.'H') then
					zdis(if,i,j)=zdis(if,i,j)*unconv_pres
				else
					zdis(if,i,j)=zdis(if,i,j)*unconv_dis
				endif
			enddo
		enddo
	enddo
	close(1); close(2)
	print *,'Modes read in'
!---------------------------------------------
!  compute filter responses, frequencies in Hz
!---------------------------------------------
	allocate(hflow(nf),hfhigh(nf))
	call hbutrw("low ",nordl,nf,2.*pi*fcl,2.*pi*fmin,2.*pi*df,hflow)
	call hbutrw("high",nordh,nf,2.*pi*fch,2.*pi*fmin,2.*pi*df,hfhigh)
!
!  apply all filters
!
	do i=1,nr
		do j=1,nmom
			zdis(nf1:nf2,i,j)=zdis(nf1:nf2,i,j)*hfhigh(1:nf)*hflow(1:nf)
		enddo
	enddo
!
	deallocate(hflow,hfhigh)
!
!  dt and nsamp
!
	call newNsampDtFrequencyTime(df,dtp,nsamp,dt)
!
!  read source time function if desired
!
	if(ifqtf.eq.1) then
		open(1,file=stfasc)
		read(1,*) nstf,dtstf
		allocate(stf(nstf))
		do i=1,nstf
			read(1,*) dum,stf(i)
		enddo
		close(1)
		if(abs(dtstf-dt)/dt.gt.1.e-6) then
			print *,'<seismogram>: dt of syns and stf inconsistent'
			print *,dtstf,dt
			stop
		endif
	endif
!
!  prepare writing sff file
!
	call new(head)
	call addSourceSFFHeader(head,'Beben',csys,.xs.info,.ys.info,depsrc,.sdate.info,.stime.info)
	call writeSFFHeader(head,1,sfffile)
!
!  back transform to time domain
!
	print *,'Back transform to time domain'
	allocate(urs(nsamp))
	last = .false.
	do i=1,nr
		do j=1,nmom
			if(i.eq.nr .and. j.eq.nmom) last = .true.
			call transformFrequencyTime(nf1,nf2,df,sigma,zdis(:,i,j),nsamp,dt,urs)
			if(ifqtf.eq.1) then
				allocate(ftemp(nsamp))
				call convolveTimeSeries(stf,urs,nsamp,ftemp)
				urs=ftemp
				deallocate(ftemp)
			endif
			if(ifwin == 1) then
				if(csys == 'S') ncv = min(nsamp,ceiling((tcut+(info.dis.i)/gv)/dt))
			else
				ncv = nsamp
			endif
			call new(dbl,ncv,dt,urs,info.netid.i,info.station.i,'XU'//(info.comp.i),last)
			call modifyDateSFFDataBlock(dbl,.syear.info,.smonth.info,.sday.info)
			call modifyTanfSFFDatablock(dbl,.therd.info)
			call addInfoSFFDatablock(dbl,csys,info.x.i,info.y.i,0.,1)
			call writeSFFDataBlock(dbl,1,last)
			call nullifyDataSFFDatablock(dbl)
		enddo
	enddo
	print *,'deallocate'
	deallocate(zdis)
	if (ifqtf == 1) deallocate(stf)
	call dealloc(head)
	deallocate(urs)
	call dealloc(info)
	print *,'Mode seismograms written to: ',trim(sfffile)
!
 end program
!------------------------------------------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   MODESEIS'
	print *,''
	print *,'Usage: modeseis [-h] [-o] [-bm] [-hp] [-lp] [-e] [-gvw] [-s] [-dt] '
	print *,'                [-qtf] [-p] [-b] [-c]'
	print *,'                infofile eigenfunction_base'
	print *,''
	print *,'Compute synthetic seismograms for given components'
	print *,'and for stations listed in an <infofile>'
	print *,'using spheroidal and/or toroidal modes read from files with'
	Print *,'names starting with <eigenfunction_base>.'
	print *,''
	print *,'<infofile> contains numeric station ids, names, locations,'
	print *,'source location, force or moment tensor components, and event time.'
	print *,'By default, the seismograms are NOT filtered.'
	print *,''
	print *,'Unit for displacement is nanometers/Nm'
	print *,'Unit for pressure is Pa/Nm'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-o    : name of output SFF file (def = modeseis.sff).'
	print *,'-bm   : work in basis moment tensor mode. Calculate for each trace in info file'
	print *,'        synthetic seismograms for six basis moment tensors and ignore moment tensor'
	print *,'        given in info file. Used for moment tensor inversion later on. (def = do not)'
	print *,'-hp   : order and corner frequency (Hz) of highpass filter'
	print *,'        such as def = <0 10.> '
	print *,'-lp   : order and corner frequency (Hz) of lowpass filter'
	print *,'        such as def = <0 100.> '
	print *,'-e    : maximum distance up to which stations'
	print *,'        are considered (def = 100 km), unit km,'
	print *,'        but if sar_csys in infofile is S, unit must be degrees.'
	print *,'-gvw  : cut length of seismogram accroding to group velocity <gv> and minimum length <tcut>:'
	print *,'        tmax = tcut + dis/gv'
	print *,'-dt   : desired sampling interval (def = 0.001 s). Actual sampling may be'
	print *,'        different because of th 2**n requirement of the FFT. It is written'
	print *,'        into the header of the seismogram file. Never use the value given here'
	print *,'        when plotting seismograms!'
	print *,'-qtf  : name of source time file'
	print *,'-p    : print spectrum for given receiver to screen'
	print *,'        Enter iprint, iprec (def = 0 1)'
	print *,'-b    : select range of branches used for seismogram calculation (def = all)'
	print *,'-c    : component if info file does not contain component information (values: Z N E R L T H)'
	print '(50(1h-))'
	return
	end
