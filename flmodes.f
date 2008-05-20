c----------------------------------------------------------------
c  Driver routine FLMODES
c
c  Description of output file in direact access mode with recl=24
c
c  Record number |   contents of record
c       1        |   nf1,nf2,df
c       2        |   rearth (double precision)
c    3 to nf+2   |   frequency index, number of record where overtone block for this frequency starts
c  --------------------------------------------------------------------------------------------------
c  After this follow the blocks with overtone data for given frequency. Record numbers are given
c  relative to the start of such a block.
c       1        |   nrt (number of overtones), if (frequency index)
c       2        |   rearth,nnod,jwd,jwu,nnod+1 (last one gives the number of records used for a single overtone eigenfunction)
c       3 ..     |   radial nodes (double precision) in groups of 3 values per record
c
c  Then follows the eigenfunction data for all overtones (record numbers relative to start of an eigenfunction block)
c
c       1        |   sngl(f),nov(n),sngl(p(n)),sngl(cg(n)),sngl(q(n))
c       2 ...    |   do j=1,nnod; write(2,rec = jsrec+j) (sngl(yeig(i,j,n)),i=1,6);	enddo
c
c  Reading the records 3 to nf+2 allows you to jump to the overtones of any frequency.
c  After reading the header information of an overtone block, you may jump to the eigenfunctions of any overtone
c  by using the information in record 2 of the header info.
c------------------------------------------------------------------------------------------------------------------------------
	include 'earthmodel.h'
	include 'flnm.h'
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'laydim.h'
	character*80 cmodel,eigfile,frechet
	double precision df,tlen,fmin,fmax,cmin,cmax,f,wlratio,rearth
	integer kernelswitch,iprint,nmin,if,i,j,k,n,ice,icf,nf1,nf2
	integer nrt,nov(nrtt),noeig,nlay
	integer n3,n3rest,n2,n2rest,n6,n6rest,jsrec,jtrec,jskrec,jtkrec
	double precision p(nrtt),cg(nrtt),q(nrtt)
	double precision yeig(nvx,nnd,nrtt),fkro(nnd,nrtt),fkvp(nnd,nrtt),fkvs(nnd,nrtt),fkd(nlayer,nrtt)
c
c  Declarations for tf_cmdline
c  Use data to provide optid and hasarg
c
	integer maxopt,manarg
	parameter(maxopt=12,manarg=1)
	logical optset(maxopt),hasarg(maxopt)
	integer lastarg,iargc
	character*6 optid(maxopt)
	character*72 optarg(maxopt)
	data optid/2h-h,2h-k,5h-fmin,5h-fmax,5h-cmin,
     1           5h-cmax,5h-tlen,2h-o,2h-p,5h-nmin,6h-noeig,4h-wlr/
	data hasarg/.false.,.true.,.true.,.true.,.true.,
     1            .true.,.true.,.true.,.true.,.true.,.false.,.true./
c
c  defaults for input parameters
c
	data iprint/0/kernelswitch/0/fmin/10.d0/nmin/5/noeig/0/
     1     fmax/100.d0/cmin/0.1d0/cmax/0.8d0/tlen/1.d0/
     1     eigfile/6heigfun/frechet/7hfrechet/wlratio/4.d0/
c---------------------------------------------------------------------
c  read commandline
c
	call tf_cmdline(1,lastarg,maxopt,optid,optarg,optset,hasarg)
c
c  help if no args or -h option
c
	if( (iargc().eq.0) .or. optset(1) .or. lastarg+manarg.ne.iargc()) then	
		call printhelp
		stop
	endif
	call getarg(lastarg+1,cmodel)
	if(cmodel(1:1).eq.'-') then
		call printhelp
		stop
	endif
c
c  set values from command lines
c
	if(optset(2)) then
		kernelswitch=1
		frechet=optarg(2)
	endif
	if(optset(3)) read(optarg(3),*) fmin
	if(optset(4)) read(optarg(4),*) fmax
	if(optset(5)) read(optarg(5),*) cmin
	if(optset(6)) read(optarg(6),*) cmax
	if(optset(7)) read(optarg(7),*) tlen
	if(optset(8)) eigfile=optarg(8)
	if(optset(9)) read(optarg(9),*) iprint
	if(optset(10)) read(optarg(10),*) nmin
	if(optset(11)) noeig=1
	if(optset(12)) read(optarg(12),*) wlratio
c--------------------------------------------------------------
c  Initialize flpsher: read model

	call flnm_read(1,cmodel)
c
c  Note: the earth model is for the reference frequency. If dispersion
c        is desired, elastic parameters should be recomputed for
c        each frequency.
c
	call earthmodel_fromflnm(flnm_fref,0)
	call earthmodel_layer
	call layer_printinfo
	call earthmodel_spline
	call earthmodel_getrearth(rearth)
c
c  Prepare loop over frequencies

	df=1.d0/tlen
	nf1=nint(fmin/df)+1
	nf2=nint(fmax/df)+1
	print *,'nf1= ',nf1,' nf2 = ',nf2,' df = ',df
	fmin=(nf1-1)*df
	fmax=(nf2-1)*df
c
c  Open files to store Green basis solutions
c  Kernels are split into separate files 
c
	call trimright(eigfile,ice)
	call trimright(frechet,icf)
c-
	open(2,file=eigfile(1:ice)//'.sph',form='unformatted', access='direct', recl=24)
	open(3,file=eigfile(1:ice)//'.tor',form='unformatted', access='direct', recl=24)
	if(kernelswitch.eq.1) then
		open(7,file=frechet(1:icf)//'.sph',form='unformatted', access='direct',recl=24)
		open(8,file=frechet(1:icf)//'.tor',form='unformatted', access='direct',recl=24)
	endif
	write(2,rec=1) nf1,nf2,sngl(df)
	write(2,rec=2) rearth
	write(3,rec=1) nf1,nf2,sngl(df)
	write(3,rec=2) rearth
	if(kernelswitch.eq.1) then
		write(7, rec=1) nf1,nf2,sngl(df)
		write(7, rec=2) rearth
		write(8, rec=1) nf1,nf2,sngl(df)
		write(8, rec=2) rearth
	endif
c
c  Loop over frequencies
c  records from 2 to nf+1 are reserved for later writing the number of the
c  first record of each frequency block
c-
	jsrec = nf2-nf1+1+3
	jtrec = nf2-nf1+1+3
	jskrec = nf2-nf1+1+3
	jtkrec = nf2-nf1+1+3
	do if=nf1,nf2
		f=(if-1)*df
		if(iprint.ge.0) then
			print *,'-------------------------------------------------'
			print *,'           FLMODES: f = ',f
			print *,'-------------------------------------------------'
		endif
		write(2,rec=2+if-nf1+1) if,jsrec     !  first record of frequency block if
		write(3,rec=2+if-nf1+1) if,jtrec
		if(kernelswitch.eq.1) then
			write(7,rec=2+if-nf1+1) if,jskrec
			write(8,rec=2+if-nf1+1) if,jtkrec
		endif
c
c  update model parameters to account for dispersion (seldamp=3)
c  if seldamp in the model file is set to 1, dispersion is suppressed
c
		if(earthmodel_seldamp.eq.3) then
			call earthmodel_fromflnm(f,if-nf1+1)
			call earthmodel_spline
		endif
c
c  compute spheroidal modes

		print *,'SPHEROIDAL MODES:'
		call modes_spheroidal(f,cmin,cmax,wlratio,nmin,iprint,noeig,
     1		                nrt,nov,p,cg,q,yeig,fkro,fkvp,fkvs,fkd)
c
c  write spheroidal eigenvalues and eigenfunctions to file
c-
		write(2,rec=jsrec) nrt,if
		if(nrt.gt.0) then
			n3 = nnod/3
			n3rest = mod(nnod,3)
			write(2, rec=jsrec+1) rearth,nnod+1
			write(2, rec=jsrec+2) nnod,jwd,jwu
			do k=1,n3
				write(2,rec=jsrec+2+k) (rnod((k-1)*3+i),i=1,3)
			enddo
			jsrec=jsrec+2+n3+1
			if(n3rest.gt.0) then
				write(2,rec=jsrec) (rnod(n3*3+i),i=1,n3rest)
				jsrec = jsrec+1
			endif
			do n=1,nrt
				write(2,rec = jsrec) sngl(f),nov(n),sngl(p(n)),sngl(cg(n)),sngl(q(n))
				do j=1,nnod
					write(2,rec = jsrec+j) (sngl(yeig(i,j,n)),i=1,6)
				enddo
				jsrec=jsrec+nnod+1
			enddo
		endif
		if(kernelswitch.eq.1) then
			call layer_getnlay(nlay)
			write(7, rec=jskrec) nrt,if
			if(nrt.gt.0) then
				n3 = nnod/3
				n3rest = mod(nnod,3)
				n2 = nnod/2
				n2rest = mod(nnod,2)
				n6 = nlay/6
				n6rest = mod(nlay,6)
				if(n6rest.eq.0) write(7, rec=jskrec+1) rearth,1+n2+n2rest+n6
				if(n6rest.gt.0) write(7, rec=jskrec+1) rearth,1+n2+n2rest+n6+1
				write(7, rec=jskrec+2) nnod,jwd,jwu
				do k=1,n3
					write(7,rec=jskrec+2+k) (rnod((k-1)*3+i),i=1,3)
				enddo
				jskrec=jskrec+2+n3+1
				if(n3rest.gt.0) then
					write(7,rec=jskrec) (rnod(n3*3+i),i=1,n3rest)
					jskrec = jskrec+1
				endif
				do n=1,nrt
					write(7, rec=jskrec) sngl(f),nov(n),sngl(p(n)),sngl(cg(n)),sngl(q(n))
					do j=1,2*n2,2
						jskrec=jskrec+1
						write(7, rec=jskrec) sngl(fkro(j,n)),sngl(fkvp(j,n)),sngl(fkvs(j,n)),
     1						                     sngl(fkro(j+1,n)),sngl(fkvp(j+1,n)),sngl(fkvs(j+1,n))
					enddo
					jskrec=jskrec+1
					if(n2rest.gt.0) then
						write(7,rec = jskrec) sngl(fkro(nnod,n)),sngl(fkvp(nnod,n)),sngl(fkvs(nnod,n))
						jskrec=jskrec+1
					endif
					write(7, rec=jskrec) nlay
					do k=1,n6
						write(7, rec=jskrec+k) (sngl(fkd((k-1)*6+i,n)),i=1,6)
					enddo
					jskrec=jskrec+n6+1
					if(n6rest.gt.0) then
						write(7, rec=jskrec) (sngl(fkd(n6*6+i,n)),i=1,n6rest)
						jskrec=jskrec+1
					endif
				enddo
			endif
		endif
c
c  compute toroidal modes

		print *,'TOROIDAL MODES:'
		call modes_toroidal(f,cmin,cmax,wlratio,nmin,iprint,noeig,
     1		              nrt,nov,p,cg,q,yeig,fkro,fkvs,fkd)
c
c  write toroidal eigenvalues and eigenfunctions to file
c-
		write(3,rec=jtrec) nrt,if
		if(nrt.gt.0) then
			n3 = nnod/3
			n3rest = mod(nnod,3)
			n2 = nnod/2
			n2rest = mod(nnod,2)
			write(3, rec=jtrec+1) rearth,n2+n2rest+1
			write(3, rec=jtrec+2) nnod,jwd,jwu
			do k=1,n3
				write(3,rec=jtrec+2+k) (rnod((k-1)*3+i),i=1,3)
			enddo
			jtrec=jtrec+2+n3+1
			if(n3rest.gt.0) then
				write(3,rec=jtrec) (rnod(n3*3+i),i=1,n3rest)
				jtrec = jtrec+1
			endif
			do n=1,nrt
				write(3,rec = jtrec) sngl(f),nov(n),sngl(p(n)),sngl(cg(n)),sngl(q(n))
!				write(6,*) sngl(f),n,sngl(p(n)),sngl(cg(n)),sngl(q(n))
				do j=1,2*n2,2
					jtrec=jtrec+1
					write(3,rec = jtrec) (sngl(yeig(i,j,n)),i=1,3),(sngl(yeig(i,j+1,n)),i=1,3)
				enddo
				jtrec=jtrec+1
!				if (if.eq.83) print *,n,nnod,n2,n2rest,jtrec,(sngl(yeig(i,nnod,n)),i=1,3)
				if(n2rest.gt.0) then
					write(3,rec = jtrec) (sngl(yeig(i,nnod,n)),i=1,3)
					jtrec=jtrec+1
				endif
			enddo
		endif
		if(kernelswitch.eq.1) then
			call layer_getnlay(nlay)
			write(8, rec=jtkrec) nrt,if
			if(nrt.gt.0) then
				n3 = nnod/3
				n3rest = mod(nnod,3)
				n2 = nnod/2
				n2rest = mod(nnod,2)
				n6 = nlay/6
				n6rest = mod(nlay,6)
				if(n6rest.eq.0) write(8, rec=jtkrec+1) rearth,1+n2+n2rest+n6
				if(n6rest.gt.0) write(8, rec=jtkrec+1) rearth,1+n2+n2rest+n6+1
				write(8, rec=jtkrec+2) nnod,jwd,jwu
				do k=1,n3
					write(8,rec=jtkrec+2+k) (rnod((k-1)*3+i),i=1,3)
				enddo
				jtkrec=jtkrec+2+n3+1
				if(n3rest.gt.0) then
					write(8,rec=jtkrec) (rnod(n3*3+i),i=1,n3rest)
					jtkrec = jtkrec+1
				endif
				do n=1,nrt
					write(8, rec=jtkrec) sngl(f),nov(n),sngl(p(n)),sngl(cg(n)),sngl(q(n))
					do j=1,2*n2,2
						jtkrec=jtkrec+1
						write(8, rec=jtkrec) sngl(fkro(j,n)),sngl(fkvs(j,n)),
     1						                     sngl(fkro(j+1,n)),sngl(fkvs(j+1,n))
					enddo
					jtkrec=jtkrec+1
					if(n2rest.gt.0) then
						write(8,rec = jtkrec) sngl(fkro(nnod,n)),sngl(fkvs(nnod,n))
						jtkrec=jtkrec+1
					endif
					write(8, rec=jtkrec) nlay
					do k=1,n6
						write(8, rec=jtkrec+k) (sngl(fkd((k-1)*6+i,n)),i=1,6)
					enddo
					jtkrec=jtkrec+n6+1
					if(n6rest.gt.0) then
						write(8, rec=jtkrec) (sngl(fkd(n6*6+i,n)),i=1,n6rest)
						jtkrec=jtkrec+1
					endif
				enddo
			endif
		endif
	enddo
	close(2)
	close(3)
	if(kernelswitch.eq.1) then
		close(7)
		close(8)
	endif
c-
	stop
	end
c------------------------------------------
c  printhelp
c-------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                     FLMODES'
	print *,''
	print *,'Usage: flmodes [-h] [-k] [-fmin] [-fmax] [-cmin] [-cmax]'
	print *,'               [-tlen] [-o] [-p] [-nmin] [-noeig] earthmodel'
	print *,''
	print *,'Compute spheroidal and toroidal l-modes for a spherically '
	print *,'symmetric earth model.'
	print *,''
	print *,'Options:'
	print *,'-h     : print help'
	print *,'-k     : compute also kernels and write to file with given basename.'
	print *,'         Files are named <basename>.sph or <basename>.tor'
	print *,'         (def = frechet)'
	print *,'-fmin  : minimum frequency in Hz (def = 5)'
	print *,'-fmax  : maximum frequency in Hz (def = 40)'
	print *,'-cmin  : minimum phase velocity in km/s (def = 0.01)'
	print *,'-cmax  : maximum phase velocity in km/s (def = 0.8)'
	print *,'-tlen  : time series length in seconds (def = 1.) '
	print *,'-o     : base name for eigenfunction output. Files are named'
	print *,'         either <basename>.sph or <basename>.tor'
	print *,'         (def = eigfun)'
	print *,'-p     : printlevel (def = 0)'
	print *,'-nmin  : minimum number of eigenfunction nodes per layer (def = 5)'
	print *,'-noeig : do not compute eigenfunctions'
	print *,'-wlr   : fraction of wavelength taken for node spacing'
	print '(50(1h-))'
	return
	end
