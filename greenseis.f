c---------------------------------------------------------------------
c  SEISMOGRAM
c---------------------------------------------------------------------
	include 'greendim.h'
	include 'green.h'
	include 'sardim.h'
	include 'sar.h'
	include 'seismodim.h'
	include 'seismo.h'
	include 'pis.h'
	include 'nkfmax.h'
	integer i,k,ifwin,ifqsp,ifqtf,ifh,kint,ico,kintmax
	integer ncv,nsamp,ierr,if,nf,nstf,precomputed_bessel,read_greenfile
	integer nordh,nordl,hour,minute
	integer iprint,iprec,ichdsv
	integer integral_jump_mom(14),integral_dsv_mom(14),integral_jump_force(7),integral_dsv_force(7)
	real zero,fch,fcl,dt,sec,dtstf,caz,saz,azi,gv,tcut
	real urs(nsampp),stf(nsampp),besselj(npp,0:2,nstatt)
	real epimax,tapfrac,dum,stfr,stfi,frequency_check,newhdur
	complex zdis(nff,nstatt,14),zoma,hflow(nff),hfhigh(nff)
	complex hfil(nff),greenkfsp(nkfmax),zsp(nff),hf(nff),zspn,zspe,zspt,zspl
	character channel*3,sdate*6,stime*10
	character integral_type_mom(14)*3,integral_type_force(7)*3
	character cjump*1,cdsv*1
	character greenfile*132,sfffile*132,infofile*80
	character stfspec*80,stfasc*80,dsvbasename*132
c
c  for sff output
c
	character wid2line*132
	integer idata(nsampp)
	logical last
c----------------------------------------------------------
c  Declarations for tf_cmdline
c  Use data to provide optid and hasarg
c
	integer maxopt,manarg
	parameter(maxopt=11,manarg=2)
	logical optset(maxopt),hasarg(maxopt)
	integer lastarg,iargc
	character*6 optid(maxopt)
	character*132 optarg(maxopt)
	data optid/2h-h,2h-o,3h-hp,3h-lp,2h-e,4h-gvw,3h-dt,4h-qsp,
     1            2h-t,4h-qtf,2h-p/
	data hasarg/.false.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,
     1             .true.,.true.,.true./
c
c  defaults for input parameters
c
	data sfffile/13hgreenseis.sff/nordh/0/
     1     fch/10./nordl/0/fcl/100./epimax/100./ifwin/0/
     1     gv/2.0/tcut/0./dt/0.001/ifqsp/0/tapfrac/0.2/
     1	    ifqtf/0/zero/1.e-6/iprint/0/iprec/1/
c
c  force: which jump basis solution is needed in the 7 integrals
c-
	data integral_jump_force/1,2,1,2,1,1,2/
c
c  moment: which jump basis solution is needed in the 14 integrals
c-
	data integral_jump_mom/1,2,3,4,1,2,3,4,1,2,1,2,3,4/
c
c  force: which motion type applies for each of the 7 integrals
c-
	data integral_type_force/'sph','sph','sph','sph','tor','tor','sph'/
c
c  moment: which motion type applies for each of the 14 integrals
c-
	data integral_type_mom/'sph','sph','sph','sph','sph','sph','sph','sph','tor','tor','tor','tor','sph','sph'/
c
c  force: which component of the displacement-stress-vector is needed in the 7 integrals
c-
	data integral_dsv_force/1,1,3,3,1,1,3/
c
c  moment: which component of the displacement-stress-vector is needed in the 14 integrals
c- 
	data integral_dsv_mom/1,1,1,1,3,3,3,3,1,1,1,1,3,3/
c---------------------------------------------------------------------
c  print dimensions

	call print_dimensions
c---------------------------------------------------------------------
c  read commandline
c
	call tf_cmdline(1,lastarg,maxopt,optid,optarg,optset,hasarg)
c
c  help if no args or -h option
c
	if( (iargc().eq.0) .or. optset(1) .or.lastarg+manarg.ne.iargc()) then	
		call printhelp
		stop
	endif
	call getarg(lastarg+1,infofile)
	call getarg(lastarg+2,dsvbasename)
	if(infofile(1:1).eq.'-') then
		call printhelp
		stop
	endif
c
c  set values from command lines
c
	if(optset(2)) sfffile=optarg(2)
	if(optset(3)) read(optarg(3),*) nordh,fch
	if(optset(4)) read(optarg(4),*) nordl,fcl
	if(optset(5)) read(optarg(5),*) epimax
	if(optset(6)) then
		ifwin=1
		read(optarg(6),*) gv,tcut
	endif
	if(optset(7)) read(optarg(7),*) dt
	if(optset(8)) then
		stfspec=optarg(8)
		ifqsp=1
	endif
	if(optset(9)) read(optarg(9),*) tapfrac
	if(optset(10)) then
		stfasc=optarg(10)
		ifqtf=1
	endif
	if(optset(11)) read(optarg(11),*) iprint,iprec
c-------------------------------------------------------
c  go
c-
	call trimright(dsvbasename,ichdsv)
c
c  Read info.eventid file
c  and convert station coordinates to epicentral
c  system -> dis(i), phi(i), dis(i) in km
c-
	call sar_readinfo(8,infofile,epimax)
	print *,sar_nr,' stations found in infofile'
	if(sar_nr.eq.0) then
		print *,'<seismogram:> No station within epimax ',epimax
		stop
	endif
	if(sar_netflag.lt.2) then
		print *,'You need an infofile with component specification'
		stop
	endif
c
c  Loop over integrals to be computed
c-
	precomputed_bessel=0
	read_greenfile=0
	if(sar_istyp.eq.1) then
		kintmax=14
	else if(sar_istyp.eq.0) then
		kintmax=7
	endif
	do kint=1,kintmax
c
c  set greenfile name
c-
		if(sar_istyp.eq.0) then
			write(cjump,'(i1)') integral_jump_force(kint)
			write(cdsv,'(i1)') integral_dsv_force(kint)
			greenfile = dsvbasename(1:ichdsv)//'.'//integral_type_force(kint)//
     1			            '.s'//'f'//'.j'//cjump//'.c'//cdsv
		else if (sar_istyp.eq.1) then
			write(cjump,'(i1)') integral_jump_mom(kint)
			write(cdsv,'(i1)') integral_dsv_mom(kint)
			greenfile = dsvbasename(1:ichdsv)//'.'//integral_type_mom(kint)//
     1			            '.s'//'m'//'.j'//cjump//'.c'//cdsv
		endif
		print *,greenfile(1:80)
c
c  Read header of Green kf-spectrum
c-
		call green_readkfhead(1,greenfile)
		call green_getnf(nf)
		if(nf*green_nwnmax.gt.2*nkfmax) then
			print *,'Dimension for greenkfsp is expected to be too small!'
			print *,'Expected: ',0.5*nf*green_nwnmax,' Dimensioned: ',nkfmax
			stop
		endif
		call green_readkfsp(1,greenkfsp)
c
c  convert sar_dis from rad to km if spherical coordinates are used, but only one time!
c-
		if(sar_csys.eq.'S'.and.read_greenfile.eq.0) then
			call sar_rad2km(green_rearth)
		endif
c
c  check frequency consistency
c-
		if(read_greenfile.eq.1) then
			if(abs((green_nf1+green_nf2)*green_dwn/frequency_check-1.).gt.zero) then
				print *, 'Frequency range and/or spacing is not uniform among DSV fk-spectra'
				print *, 'Error found with file: ',greenfile,' and integral number ',kint
				print *,(green_nf1+green_nf2)*green_dwn,frequency_check
			endif
		endif
c
c  compute kint-th wavenumber integral for all receivers
c 
		do i=1,sar_nr
			if(precomputed_bessel.eq.0) call green_bessel(sar_dis(i),besselj(1,0,i))
			if(sar_istyp.eq.1) call green_wnint_mom(kint,sar_dis(i),tapfrac,greenkfsp,besselj(1,0,i),zdis(1,i,kint))
			if(sar_istyp.eq.0) call green_wnint_force(kint,sar_dis(i),tapfrac,greenkfsp,besselj(1,0,i),zdis(1,i,kint))
		enddo
		precomputed_bessel=1
		read_greenfile=1
c
c  calculate frequency check number for current DSV fk-spectrum
c-
		frequency_check=(green_nf1+green_nf2)*green_dwn
	enddo
c
c  read response file
c-
	zoma=cmplx(2.*pi*(green_nf1-1)*green_df,-green_sigma)
	print *,'Frequency range: ',(green_nf1-1)*green_df,(green_nf2-1)*green_df,green_df
c
c  low and high-pass filter
c-
	call hbutcw("low ",nordl,nf,2.*pi*fcl,zoma,2.*pi*green_df,hflow)
	call hbutcw("high",nordh,nf,2.*pi*fch,zoma,2.*pi*green_df,hfhigh)	
c
c  Source spectrum
c-
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
c  multiply all filters
c
	do if=1,nf
		hfil(if)=hflow(if)*hfhigh(if)
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
c
c  Computation of seismograms
c  determine seismo_dt consistent with desired dt and power of 2 requirement
c-
	call seismo_oversampling(green_df,dt)
c
c  use half duration to create a source time function if none else is given
c-
	if(ifqtf.eq.0.and.ifqsp.eq.0) then
		call seismo_stffromhdursinq(sar_hdur,nstf,stf,newhdur)
		print *,'Convolve impulse response using hdur = ',newhdur
	endif
c
c  prepare writing to SFF, open file
c-
	call sar_eventsdate(sdate)
	call sar_eventstime(stime)
	call sar_eventhms(hour,minute,sec)
	call sff_New(1,sfffile,ierr)
	if(sar_csys.eq.'S') then
		call sff_WOpenS(1,sfffile,'Earthquake','S',
     1	                       90.-sar_xsrc,sar_ysrc,sar_depsrc,sdate,stime,ierr)
	else
		call sff_WOpenS(1,sfffile,'Sledge hammer','C',
     1	                      sar_xsrc,sar_ysrc,sar_depsrc,sdate,stime,ierr)
	endif
c
c  loop over receivers
c
	last=.false.
	do i=1,sar_nr
		if(i.eq.sar_nr) last=.true.
c
c  prepare rotation
c-
		call sar_propdir(sar_xsta(i),sar_ysta(i),azi)
		caz=cos(azi)
		saz=sin(azi)
c
c  displacement components for force
c-
		if(sar_istyp.eq.0) then
			if(sar_comp(i).eq.'Z') then
				do if=green_nf1,green_nf2
					zsp(if) = zdis(if,i,1)*sar_force(1)
     1						  +zdis(if,i,2)*(sar_force(2)*cos(sar_phi(i))+sar_force(3)*sin(sar_phi(i)))
				enddo
			else if(sar_comp(i).eq.'N' .or. sar_comp(i).eq.'E') then
				do if=green_nf1,green_nf2
					zspl = zdis(if,i,3)*sar_force(1)
     1					      +(zdis(if,i,4)+zdis(if,i,5))*(sar_force(2)*cos(sar_phi(i))+sar_force(3)*sin(sar_phi(i)))
					zspt = (zdis(if,i,6)+zdis(if,i,7))*(sar_force(2)*sin(sar_phi(i))-sar_force(3)*cos(sar_phi(i)))
					if(sar_csys.eq.'S') then
						zspn=-zspl*caz+zspt*saz
						zspe=+zspl*saz+zspt*caz
					else
						zspn=+zspl*caz-zspt*saz
						zspe=+zspl*saz+zspt*caz
					endif
					if(sar_comp(i).eq.'N') zsp(if)=zspn
					if(sar_comp(i).eq.'E') zsp(if)=zspe
				enddo
		   	endif
c
c  displacement components for moment tensor
c-
		else if(sar_istyp.eq.1) then
			if(sar_comp(i).eq.'Z') then
				do if=green_nf1,green_nf2
					zsp(if) = zdis(if,i,1)*sar_mt(1)
     1					         +zdis(if,i,2)*(sar_mt(2)+sar_mt(3))
     1					         +zdis(if,i,3)*(sar_mt(4)*cos(sar_phi(i))+sar_mt(5)*sin(sar_phi(i)))
     1					         +zdis(if,i,4)*((sar_mt(3)-sar_mt(2))*cos(2.*sar_phi(i))-2.*sar_mt(6)*sin(2.*sar_phi(i)))
				enddo
			else if(sar_comp(i).eq.'N' .or. sar_comp(i).eq.'E') then
				do if=green_nf1,green_nf2
					zspl = zdis(if,i,5)*sar_mt(1)
     1					      +zdis(if,i,6)*( sar_mt(2)+sar_mt(3) )
     1					      +zdis(if,i,7)*( sar_mt(4)*cos(sar_phi(i))+sar_mt(5)*sin(sar_phi(i)) )
     1					      +zdis(if,i,8)*( (sar_mt(3)-sar_mt(2))*cos(2.*sar_phi(i))-2.*sar_mt(6)*sin(2.*sar_phi(i)) )
     1					      +zdis(if,i,9)*( sar_mt(4)*cos(sar_phi(i))+sar_mt(5)*sin(sar_phi(i)) )
     1					      +zdis(if,i,10)*( (sar_mt(3)-sar_mt(2))*cos(2.*sar_phi(i))-2.*sar_mt(6)*sin(2.*sar_phi(i)) )
					zspt = zdis(if,i,11)*( sar_mt(5)*cos(sar_phi(i))-sar_mt(4)*sin(sar_phi(i)) )
     1					      +zdis(if,i,12)*( (sar_mt(3)-sar_mt(2))*sin(2.*sar_phi(i))+2.*sar_mt(6)*cos(2.*sar_phi(i)) )
     1					      +zdis(if,i,13)*( sar_mt(5)*cos(sar_phi(i))-sar_mt(4)*sin(sar_phi(i)) )
     1					      +zdis(if,i,14)*( (sar_mt(3)-sar_mt(2))*sin(2.*sar_phi(i))+2.*sar_mt(6)*cos(2.*sar_phi(i)) )
					if(sar_csys.eq.'S') then
						zspn=-zspl*caz+zspt*saz
						zspe=+zspl*saz+zspt*caz
					else
						zspn=+zspl*caz-zspt*saz
						zspe=+zspl*saz+zspt*caz
					endif
					if(sar_comp(i).eq.'N') zsp(if)=zspn
					if(sar_comp(i).eq.'E') zsp(if)=zspe
				enddo
			endif
		endif
c
c  print something
c-
		if(iprint.gt.0.and.iprec.eq.i) then
			print *,sar_xsrc,sar_ysrc
			print *,sar_xsta(i),sar_ysta(i)
			print *,(zsp(if),if=green_nf1,green_nf2)
		endif
c
c  apply all filters
c-
		do if=green_nf1,green_nf2
			ifh=if-green_nf1+1
			zsp(if)=zsp(if)*hfil(ifh)
		enddo
c
c  transform to time domain
c-
		call seismo_f2t(green_nf1,green_nf2,green_df,ifwin,sar_dis(i)/gv+tcut,
     1		                green_sigma,zsp,nsamp,urs)
		ncv=nsamp
c
c  convolve with source time function
c  overwrites urs
c-
		if(ifqsp.ne.1) call causfalt(nsampp,nsamp,nstf,ncv,seismo_dt,urs,stf)
c
c  prepare wid2line
c-
		channel='XU'//sar_comp(i)
		call sff_PrepWid2(ncv,1./seismo_dt,sar_name(i),sar_year,sar_month,sar_day,hour,
     1		                  minute,channel,'NSP','NSP',sec,-1.,-1.,-1.,-1.,
     1		                  wid2line,ierr)
		print *,wid2line(1:80)

		if(iprint.gt.0.and.iprec.eq.i) then
			print *,green_df,nf,green_sigma,sar_dis(i)/gv+tcut
			print *,nsamp,seismo_dt,seismo_nsamp,nstf,stf(1)
			print *,sar_xsrc,sar_ysrc
			print *,sar_xsta(i),sar_ysta(i)
			print *,(urs(k),k=1,ncv)
		endif
c
c  write datablock
c-
		if(sar_csys.eq.'S') then
			call sff_WTraceI(1+ico,wid2line,ncv,urs,idata,last,'S',
     1		                       90.-sar_xsta(i),sar_ysta(i),0.,1,ierr)	
		else
			call sff_WTraceI(1+ico,wid2line,ncv,urs,idata,last,'C',
     1		                       sar_xsta(i),sar_ysta(i),0.,1,ierr)	
		endif
	enddo
	stop
	end
c-----------------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                   GREENSEIS'
	print *,''
	print *,'Usage: greenseis  [-h] [-o] [-hp] [-lp] [-e] [-s] [-dt] '
	print *,'                  [-qsp] [-st] [-t] [-r] [-qtf] [-p]'
	print *,'                  infofile dsvbasename'
	print *,''
	print *,'Compute synthetic seismograms for receivers listed in an <infofile>'
	print *,'using DSV frequency-wavenumber spectra read from files with'
	Print *,'names starting with <dsvbasename>.'
	print *,''
	print *,'<infofile> contains numeric station ids, names, components, network, locations,'
	print *,'source location, force or moment tensor components, and event time.'
	print *,''
	print *,'Depending on the chosen components the required fk-spectra files are'
	print *,'automatically selected following the filename conventions in flgevask.'
	print *,'But make sure that frequency and wavenumber range and spacing are'
	print *,'identical for all fk-spectra needed.'
	print *,'' 
	print *,'By default, the seismograms are NOT filtered.'
	print *,''
	print *,'Unit for displacement is nanometers/Nm'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h    : print help'
	print *,'-o    : name of output SFF file (def = greenseis.sff)'
	print *,'-hp   : order and corner frequency (Hz) of highpass filter'
	print *,'        such as def = <0 10.> '
	print *,'-lp   : order and corner frequency (Hz) of lowpass filter'
	print *,'        such as def = <0 100.> '
	print *,'-e    : maximum distance up to which stations'
	print *,'        are considered (def = 100 km), unit km,'
	print *,'        but if sar_csys in infofile is S, unit must be degrees.'
	print *,'-gvw  : cut length of seismogram according to group velocity <gv> and minimum length <tcut>:'
	print *,'        tmax = tcut + dis/gv'
	print *,'-dt   : desired sampling interval (def = 0.001 s). Actual sampling may be'
	print *,'        different because of th 2**n requirement of the FFT. It is written'
	print *,'        into the header of the seismogram file. Never use the value given here'
	print *,'        when plotting seismograms!'
	print *,'-qsp  : name of source spectrum file (either use -qsp or -qtf or nothing!)'
	print *,'-qtf  : name of source time file (either use -qsp or -qtf or nothing!)'
	print *,'-t    : fraction of wavenumber range used for taper (def=0.2)'
	print *,'-p    : print spectrum for given receiver and component to screen'
	print *,'        Enter iprint, ipco, iprec (def = 0 1 1)'
	print '(50(1h-))'
	return
	end
c--------------------------------------------------
c  print dimensions

	subroutine print_dimensions
	include 'greendim.h'
	include 'sardim.h'
	include 'seismodim.h'
	include 'nkfmax.h'

	print *,'greendim.h:  nff    = ',nff
	print *,'greendim.h:  npp    = ',npp
	print *,'sardim.h:    nstatt = ',nstatt
	print *,'seismodim.h: nsampp = ',nsampp
	print *,'nkfmax.h:    nkfmax = ',nkfmax

	return
	end
