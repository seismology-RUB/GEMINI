c-----------------------------------------------
	include 'seismodim.h'
	include 'seismo.h'
	include 'sardim.h'
	include 'sar.h'
	include 'nodesdim.h'
	include 'eigdim.h'
	include 'pis.h'
	integer nff
	parameter(nff=1024)
c-
	integer icheig,nordh,nordl,nf1,nf2,if,ids,nrt_tor,nrt_sph,nf,i,j
	integer ifwin,ifsilent,ifqsp,ifqtf,iprint,iprec
	integer jsrec,jtrec,nrec_sph,nrec_tor,n6,n6rest,k,n2,n2rest,nrt2
	integer n,nov,ifh,nstf,hour,minute,nsamp,ncv,nnod,nrtmin,nrtmax,jwd,jwu
	real rearth,fch,fcl,epimax,gv,tcut,dt,df,fmin,fmax,f,zsbs,zsas,dtstf,vl
	real p,cg,qatt,us,ups,vs,vps,ws,wps,unconv_dis,unconv_pres,sec
	real dum,stfr,stfi,newhdur
	real yeig(6,nnd),ur,vr,wr,rr
	real stf(nsampp),urs(nsampp),znod(nnd)
	complex hflow(nff),hfhigh(nff),hf(nff)
	complex psi(0:2),zdis(nff,nstatt)
	character sdate*6,stime*10,channel*3
	character stfspec*80,stfasc*80,sfffile*80
	character eigbase*80,infofile*80,sfftoc*80
c--------------------------
c  for sff output
c-------------------------
	integer idata(nsampp),ierr
	character wid2line*132
	logical last
c----------------------------------------------------------
c  Declarations for tf_cmdline
c  Use data to provide optid and hasarg
c
	integer maxopt,manarg
	parameter(maxopt=14,manarg=2)
	logical optset(maxopt),hasarg(maxopt)
	integer lastarg,iargc
	character*6 optid(maxopt)
	character*132 optarg(maxopt)
	data optid/2h-h,2h-o,3h-hp,3h-lp,2h-e,4h-gvw,2h-s,3h-dt,4h-qsp,
     1            4h-qtf,2h-p,2h-b,2h-c,4h-toc/
	data hasarg/.false.,.true.,.true.,.true.,.true.,.true.,.false.,.true.,.true.,
     1             .true.,.true.,.true.,.true.,.true./
c-----------------------------------------------------------
c  defaults for input parameters
c-
	data nordh/0/sfffile/12hmodeseis.sff/sfftoc/12hmodeseis.toc/
     1     fch/10./nordl/0/fcl/100./epimax/100./ifwin/0/ifsilent/0/
     1     gv/2./tcut/0./dt/0.001/ifqsp/0/ifqtf/0/iprint/0/
     1     iprec/1/nrtmin/1/nrtmax/1000/unconv_dis/1.e-6/unconv_pres/1.e-9/
c---------------------------------------------------------------------
c  read commandline
c-
	call tf_cmdline(1,lastarg,maxopt,optid,optarg,optset,hasarg)
c
c  help if no args or -h option
c-
	if( (iargc().eq.0) .or. optset(1) .or.lastarg+manarg.ne.iargc()) then	
		call printhelp
		stop
	endif
	call getarg(lastarg+1,infofile)
	call getarg(lastarg+2,eigbase)
	if(infofile(1:1).eq.'-'.or.eigbase(1:1).eq.'-') then
		call printhelp
		stop
	endif
c
c  set values from command lines
c-
	if(optset(2)) sfffile=optarg(2)
	if(optset(3)) read(optarg(3),*) nordh,fch
	if(optset(4)) read(optarg(4),*) nordl,fcl
	if(optset(5)) read(optarg(5),*) epimax
	if(optset(6)) then
		ifwin=1
		read(optarg(6),*) gv,tcut
	endif
	if(optset(7)) ifsilent=1
	if(optset(8)) read(optarg(8),*) dt
	if(optset(9)) then
		stfspec=optarg(9)
		ifqsp=1
	endif
	if(optset(10)) then
		stfasc=optarg(10)
		ifqtf=1
	endif
	if(optset(11)) read(optarg(11),*) iprint,iprec
	if(optset(12)) read(optarg(12),*) nrtmin,nrtmax
	if(optset(14)) sfftoc=optarg(14)
c-------------------------------------------------------
	call trimright(eigbase,icheig)
c
c  Read info.eventid file
c  and convert station coordinates to epicentral
c  system -> dis(i), phi(i), dis(i) in km
c-
	call sar_readinfo(1,infofile,epimax)
	print *,sar_nr,' stations found in infofile'
	if(sar_nr.eq.0) then
		print *,'<seismogram:> No station within epimax ',epimax
		stop
	endif
	if(sar_netflag.lt.2) then
		if(.not.optset(13)) then
			print *,'You need either an infofile with component specification'
			print *,'or you need to specify a component on the command line'
			stop
		else
			do i=1,sar_nr
				sar_comp(i)=optarg(13)(1:1)
			enddo
		endif
	endif
c---------------------------------------------------------------
c  Prepare loop over frequencies
c---------------------------------------------------------------
	open(1,file=eigbase(1:icheig)//'.sph',form='unformatted',access='direct',recl=24)
	open(2,file=eigbase(1:icheig)//'.tor',form='unformatted',access='direct',recl=24)
	read(1,rec=1) nf1,nf2,df,rearth
	read(2,rec=1) nf1,nf2,df,rearth
	nf=nf2-nf1+1
	fmin=(nf1-1)*df
	fmax=(nf2-1)*df
	print *,'Frequencies: ',fmin,fmax,nf1,nf2,df
	nrtmin=max(nrtmin,1)
c
c  convert sar_dis from rad to km if spherical coordinates are used
c
	if(sar_csys.eq.'S') then
		call sar_rad2km(rearth)
	endif
c---------------------------------------------------------------
c  Loop over frequencies, now in Hz
c  znod contains depths!!
c  ids+1 belongs to point above source !
c  ups, vps and wps output of spline_interpol_2d is derivative
c  with respect to z -> take negative value which is derivative wrt r.
c---------------------------------------------------------------------
	jsrec = nf+2
	jtrec = nf+2
	do if=nf1,nf2
		do i=1,sar_nr
			zdis(if,i)=0.
		enddo
c
c  spheroidal stuff
c  go on to toroidal modes if nrt_sph=0
c  take displacement and stress eigenfunctions at top of model (nnod) 
c  if there is a water layer on top take both (U,V,R) at ocean bottom (jwd)!
c-
		read(1,rec=jsrec) nrt_sph
		if(nrt_sph.eq.0) then
			jsrec=jsrec+1
			goto 101
		endif
		nrt2=min(nrt_sph,nrtmax)
		read(1,rec=jsrec+1) f,nnod,jwd,jwu,nrec_sph
		n6=nnod/6
		n6rest=mod(nnod,6)
		do k=1,n6
			read(1,rec=jsrec+1+k) (znod((k-1)*6+i),i=1,6)
		enddo
		jsrec=jsrec+1+n6+1
		if(n6rest.gt.0) then
			read(1,rec=jsrec) (znod(n6*6+i),i=1,n6rest)
			jsrec=jsrec+1
		endif
		call locate(sar_depsrc,nnod,znod,ids)
		if(ids.eq.nnod) ids=ids-1
		zsbs = znod(ids)
		zsas = znod(ids+1)
		if(nnod.gt.nnd) then
			write(6,'(a)') '<modeseis>: There are more nodes than dimensioned'
			stop
		endif
		jsrec = jsrec+(nrtmin-1)*nrec_sph
		do 10 n=nrtmin,nrt2
			read(1,rec=jsrec) nov,p,cg,qatt
			do j=1,nnod
				read(1,rec=jsrec+j) (yeig(i,j),i=1,6)
			enddo
			jsrec=jsrec+nnod+1
			call spline_interpol_2p(sar_depsrc,zsas,zsbs,yeig(1,ids+1),yeig(1,ids),-yeig(5,ids+1),
     1		                               -yeig(5,ids),us,ups)
			call spline_interpol_2p(sar_depsrc,zsas,zsbs,yeig(3,ids+1),yeig(3,ids),-yeig(6,ids+1),
     1		                               -yeig(6,ids),vs,vps)
c			print *,n,ids,yeig(1,ids),us,yeig(1,ids+1),yeig(5,ids),-ups,yeig(5,ids+1)
c			print *,n,ids,yeig(3,ids),vs,yeig(3,ids+1),yeig(6,ids),-vps,yeig(6,ids+1)
			ur=yeig(1,nnod)
			rr=yeig(2,nnod)
			vr=yeig(3,nnod)
			if(jwd.gt.0) then
				ur=yeig(1,jwd)
				rr=yeig(2,jwd)
				vr=yeig(3,jwd)
			endif
			vl=2.*pi*f*rearth*p
			call excoso_sph(sar_mt,dble(rearth)-dble(sar_depsrc),vl,us,vs,-ups,-vps,psi)
			do i=1,sar_nr
				call add_sph_displ(i,f,p,cg,qatt,psi,rearth,ur,vr,rr,zdis(if,i))
			enddo
 10		continue
c
c  toroidal stuff
c  take displacement at the surface (nnod)
c  if a water layer is present, evaluate displacement at ocean bottom
c-
 101		read(2,rec=jtrec) nrt_tor
		if(nrt_tor.eq.0) then
			jtrec=jtrec+1
			goto 102
		endif
		nrt2=min(nrt_tor,nrtmax)
		read(2,rec=jtrec+1) f,nnod,jwd,jwu,nrec_tor
		n6 = nnod/6
		n6rest = mod(nnod,6)
		n2 = nnod/2
		n2rest = mod(nnod,2)
		do k=1,n6
			read(2,rec=jtrec+1+k) (znod((k-1)*6+i),i=1,6)
		enddo
		jtrec=jtrec+1+n6+1
		if(n6rest.gt.0) then
			read(2,rec=jtrec) (znod(n6*6+i),i=1,n6rest)
			jtrec = jtrec+1
		endif
		call locate(sar_depsrc,nnod,znod,ids)
		if(ids.eq.nnod) ids=ids-1
		zsbs = znod(ids)
		zsas = znod(ids+1)
		if(nnod.gt.nnd) then
			write(6,'(a)') '<modeseis>: There are more nodes than dimensioned'
			stop
		endif
		jtrec=jtrec+(nrtmin-1)*nrec_tor
		do 20 n=nrtmin,nrt2
			read(2,rec = jtrec) nov,p,cg,qatt
			do j=1,2*n2,2
				jtrec=jtrec+1
				read(2,rec = jtrec) (yeig(i,j),i=1,3),(yeig(i,j+1),i=1,3)
			enddo
			jtrec=jtrec+1
			if(n2rest.gt.0) then
				read(2,rec = jtrec) (yeig(i,nnod),i=1,3)
				jtrec=jtrec+1
			endif
			call spline_interpol_2p(sar_depsrc,zsas,zsbs,yeig(1,ids+1),yeig(1,ids),-yeig(3,ids+1),
     1		                               -yeig(3,ids),ws,wps)
c			print *,n,ids,yeig(1,ids),ws,yeig(1,ids+1),yeig(3,ids),-wps,yeig(3,ids+1)
			wr=yeig(1,nnod)
			if(jwd.gt.0) wr=yeig(1,jwd)
			vl=2.*pi*f*rearth*p
			call excoso_tor(sar_mt,dble(rearth)-dble(sar_depsrc),vl,ws,-wps,psi)
			do i=1,sar_nr
				call add_tor_displ(i,f,p,cg,qatt,psi,rearth,wr,zdis(if,i))
			enddo
 20		continue
c
c  convert to nanometers or Pa (for H)
c-
 102		do i=1,sar_nr
			if(sar_comp(i).eq.'H') then
				zdis(if,i)=zdis(if,i)*unconv_pres
			else
				zdis(if,i)=zdis(if,i)*unconv_dis
			endif
		enddo
	enddo
	close(1)
	close(2)
c---------------------------------------------
c  compute filter responses, frequencies in Hz
c---------------------------------------------
	call hbutrw("low ",nordl,nf,2.*pi*fcl,2.*pi*fmin,2.*pi*df,hflow)
	call hbutrw("high",nordh,nf,2.*pi*fch,2.*pi*fmin,2.*pi*df,hfhigh)
c
c  Source spectrum
c
	if(ifqsp.eq.1) then
		open(1,file=stfspec)
		do if=1,nf
			read(1,*) dum,dum,stfr,stfi
			hf(if)=cmplx(stfr,stfi)
		enddo
		close(1)
	else
		do if=1,nf
			hf(if)=(1.0,0.0)
		enddo
	endif
c
c  apply all filters
c-
	do i=1,sar_nr
		do if=nf1,nf2
			ifh=if-nf1+1
			zdis(if,i)=zdis(if,i)*hfhigh(ifh)*hflow(ifh)*hf(ifh)
		enddo
	enddo
c
c  read source time function if desired
c
	if(ifqtf.eq.1) then
		open(1,file=stfasc)
		read(1,*) nstf,dtstf
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
c----------------------------------------------------------------------
c  Computation of seismograms
c  determine seismo_dt consistent with desired dt and power of 2 requirement
c-
	call seismo_oversampling(df,dt)
c
c  use half duration to create a source time function if none else is given
c-
	if(ifqtf.eq.0.and.ifqsp.eq.0) then
		call seismo_stffromhdursinq(sar_hdur,nstf,stf,newhdur)
		print *,'Convolve impulse response using hdur = ',newhdur
	endif
c
c  prepare writing to SFF
c
	call sar_eventsdate(sdate)
	call sar_eventstime(stime)
	call sar_eventhms(hour,minute,sec)
	if(ifsilent.eq.0) print *,'Write seismogram to ',sfffile
	call sff_NewToc(1,2,sfffile,sfftoc,ierr)
	if(sar_csys.eq.'S') then
		call sff_WOpenTocS(1,2,sfffile,sfftoc,'Earthquake          ','S',
     1	                       90.-sar_xsrc,sar_ysrc,sar_depsrc,sdate,stime,ierr)
	else
		call sff_WOpenTocS(1,2,sfffile,sfftoc,'Sledge hammer       ','C',
     1	                       sar_xsrc,sar_ysrc,sar_depsrc,sdate,stime,ierr)
	endif
c
c  transform to time domain
c
	last=.false.
	do i=1,sar_nr
		call seismo_f2t(nf1,nf2,df,ifwin,sar_dis(i)/gv+tcut,0.0,zdis(1,i),nsamp,urs)
		ncv=nsamp
c
c  convolve with source time function
c  overwrites urs
c-
		if(ifqsp.ne.1) call causfalt(nsampp,nsamp,nstf,ncv,seismo_dt,urs,stf)

c  prepare wid2line
c-
		channel='XU'//sar_comp(i)
		call sff_PrepWid2(ncv,1./seismo_dt,sar_name(i),sar_year,sar_month,sar_day,hour,
     1		                  minute,channel,'NSP','NSP',sec,-1.,-1.,-1.,-1.,
     1		                  wid2line,ierr)
		if(ifsilent.eq.0) print *,wid2line(1:80)

		if(iprint.gt.0.and.iprec.eq.i) then
			print *,df,nf,sar_dis(i)/gv+tcut
			print *,nsamp,seismo_dt,seismo_nsamp,nstf,stf(1)
			print *,sar_xsrc,sar_ysrc
			print *,sar_comp(i),sar_xsta(i),sar_ysta(i),sar_dis(i)
			print *,(urs(j),j=1,ncv)
		endif
c
c  write datablock
c-
		if(i.eq.sar_nr) last=.true.
		if(sar_csys.eq.'S') then
			call sff_WTraceTocI(1,2,wid2line,ncv,urs,idata,last,'S',
     1		                           90.-sar_xsta(i),sar_ysta(i),0.,1,ierr)	
		else
			call sff_WTraceTocI(1,2,wid2line,ncv,urs,idata,last,'C',
     1		                           sar_xsta(i),sar_ysta(i),0.,1,ierr)	
		endif
 12		continue
	enddo
	stop
	end
c-----------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   MODESEIS'
	print *,''
	print *,'Usage: modeseis [-h] [-o] [-toc] [-hp] [-lp] [-e] [-gvw] [-s] [-dt] '
	print *,'                  [-qsp] [-qtf] [-p] [-b] [-c]'
	print *,'                  infofile eigenfunction_base'
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
	print *,'-toc  : name of sff table of contents file (def = modeseis.toc)'
	print *,'-hp   : order and corner frequency (Hz) of highpass filter'
	print *,'        such as def = <0 10.> '
	print *,'-lp   : order and corner frequency (Hz) of lowpass filter'
	print *,'        such as def = <0 100.> '
	print *,'-e    : maximum distance up to which stations'
	print *,'        are considered (def = 100 km), unit km,'
	print *,'        but if sar_csys in infofile is S, unit must be degrees.'
	print *,'-gvw  : cut length of seismogram accroding to group velocity <gv> and minimum length <tcut>:'
	print *,'        tmax = tcut + dis/gv'
	print *,'-s    : be silent'
	print *,'-dt   : desired sampling interval (def = 0.001 s). Actual sampling may be'
	print *,'        different because of th 2**n requirement of the FFT. It is written'
	print *,'        into the header of the seismogram file. Never use the value given here'
	print *,'        when plotting seismograms!'
	print *,'-qsp  : name of source spectrum file (either use -qsp or -qtf or nothing!)'
	print *,'-qtf  : name of source time file (either use -qsp or -qtf or nothing!)'
	print *,'-p    : print spectrum for given receiver to screen'
	print *,'        Enter iprint, iprec (def = 0 1)'
	print *,'-b    : select range of branches used for seismogram calculation (def = all)'
	print *,'-c    : component if info file does not contain component information (values: Z N E R L T H)'
	print '(50(1h-))'
	return
	end
