c---------------------------------------------------------------
c $Id: flgevask.f,v 1.7 2003/07/04 06:45:53 wolle Exp $
c---------------------------------------------------------------
	include 'pi.h'
	include 'earthmodel.h'
	include 'parameter.h'
	include 'greendim.h'
	double complex gfsrsph(4,4,npp),sksph(nkk,4,nparr,4,npp)
	double complex gfsrtor(2,2,npp),sktor(nkk,2,nparr,2,npp)
	double complex zsrcsph(4,4),zsrctor(2,2)
	double precision tlen,xlen,fmin,fmax,f,df,sigma,pmax,wn,dwn,fnode,om
	double precision zstore,zs,ze,wn2,wnmax,frac,fmax,fmin
	integer i,j,m,mm,if,nf1,nf2,iwn,nwn,icg,icp,k,iunit,nsrcsph,nsrctor,nwnmax
	integer iprint,ncomp,comps(6),nnmin,istyp,selflag,npar,pars(nparr)
	integer ncompsph,compssph(4),ncomptor,compstor(2),npartor,parstor(nparr)
	integer jump_spher(4), jump_tor(2)
	character*80 cmodel,greenfile,partials
	character ckomp*1,cpar*1,sourcetype*1,cjump*1
c----------------------------------------------------------
c  Declarations for tf_cmdline
c  Use data to provide optid and hasarg
c
	integer maxopt,manarg
	parameter(maxopt=19,manarg=1)
	logical optset(maxopt),hasarg(maxopt)
	integer lastarg,iargc
	character*6 optid(maxopt)
	character*72 optarg(maxopt)
	data optid/2h-h,2h-p,5h-tlen,5h-fmin,5h-fmax,5h-pmax,4h-sig,3h-zs,
     1	           3h-ze,2h-s,3h-od,3h-op,5h-xlen,5h-nmin,
     2	           3h-sf,2h-c,2h-m,5h-fnod,5h-wmrg/
	data hasarg/.false.,.true.,.true.,.true.,.true.,.true.,.true.,
     1            .true.,.true.,.true.,.true.,.true.,.true.,.true.,
     2	          .true.,.true.,.true.,.true.,.true./
c-----------------------------------------------------------
c  defaults for input parameters
c
	data iprint/0/tlen/1./sigma/0.d0/fmin/5./fmax/40./
     1	     pmax/5./zs/1./greenfile/5hgreen/
     2	     xlen/1./nnmin/5/istyp/2/ze/0./partials/8hpartials/ncomp/1/
     4	     comps/1,2,3,4,5,6/selflag/1/npar/1/pars/1,2,3,4,5/fnode/5./frac/0.1/
c---------------------------------------------------------------------
c  print dimensions to standard output

	call print_dimensions
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
	if(optset(2)) read(optarg(2),*) iprint
	if(optset(3)) read(optarg(3),*) tlen
	if(optset(4)) read(optarg(4),*) fmin
	if(optset(5)) read(optarg(5),*) fmax
	if(optset(6)) read(optarg(6),*) pmax
	if(optset(7)) read(optarg(7),*) sigma
	if(optset(8)) read(optarg(8),*) zs
	if(optset(9)) read(optarg(9),*) ze
	if(optset(10)) read(optarg(10),*) istyp
	if(optset(11)) greenfile=optarg(11)
	if(optset(12)) partials=optarg(12)
	if(optset(13)) read(optarg(13),*) xlen
	if(optset(14)) read(optarg(14),*) nnmin
	if(optset(15)) read(optarg(15),*) selflag
	if(optset(16)) then
		call isplitstring(optarg(16),' ',ncomp,comps)
	endif
	if(optset(17)) then
		call isplitstring(optarg(17),' ',npar,pars)
	endif
	if(optset(18)) read(optarg(18),*) fnode
	if(optset(19)) read(optarg(19),*) frac
c--------------------------------------------------------------
c  go
c
	sigma=sigma/tlen
	df=1.d0/tlen
	if(fmin.lt.df) fmin=df
	nf1=nint(fmin/df)+1
	nf2=nint(fmax/df)+1
	fmin=(nf1-1)*df
	fmax=(nf2-1)*df
	dwn=2.*pi/xlen
	nnmin=max(nnmin,2)
c
	if(istyp.eq.1.or.istyp.eq.5) sourcetype='m'
	if(istyp.eq.2.or.istyp.eq.3.or.istyp.eq.4) sourcetype='f'
c
c           Echo input
c
	print '(a,i2)','Control-output print level:',iprint
	print '(a,1p,2g15.5)','Seism. length, df : ',tlen,df
	print '(a,1p,g15.5)','Sigma: ',sigma
	print '(a,1p,2g13.3,2i5)','fmin,fmax,nf1,nf2 ',fmin,fmax,nf1,nf2
	print *,'pmax,Length of profile: ',pmax,xlen
	print '(2a)','Earth model: ',cmodel
	print '(a,g15.6)','Source depth (m): ',zs
	print *,'Green components: ',(comps(i),i=1,ncomp)
c
c  Initialize flgevas
c
	call flgevas_init(cmodel,fnode,zstore,zs,ze,nnmin,ncomp,comps,npar,pars,
     1	            ncompsph,compssph,ncomptor,compstor,npartor,parstor,
     1			istyp,nsrcsph,jump_spher,nsrctor,jump_tor,iprint)
c
c  print some information about files to be created
c
	print *,'Spheroidal basis solutions neded: ',(jump_spher(i),i=1,nsrcsph)
	print *,'Toroidal basis solutions needed:  ',(jump_tor(i),i=1,nsrctor)
	print *,'Spheroidal DSV-components needed: ',(compssph(i),i=1,ncompsph)
	print *,'Toroidal DSV-components needed:   ',(compstor(i),i=1,ncomptor)
c
c  Open files to store Green basis solutions
c  Kernels are split into separate files 
c
	call trimright(greenfile,icg)
	call trimright(partials,icp)
c
	wnmax=(nf2-1)*df*2.*pi*pmax
	nwnmax=int(wnmax/dwn)+1
	if(nwnmax.gt.npp) then
		print *,'More wavenumbers required than dimensioned'
		print *,nwnmax,npp
		stop
	endif
c
c  spheroidal results
c
	if(nsrcsph.gt.0) then
	do i=1,ncompsph
		write(ckomp,'(i1)') compssph(i)
		do m=1,nsrcsph
			write(cjump,'(i1)') jump_spher(m)
			iunit=6+(i-1)*nsrcsph*(npar+1)+(m-1)*(npar+1)+1
			open(iunit,file=greenfile(1:icg)//'.sph'//'.s'//sourcetype//'.j'//
     1		     cjump//'.c'//ckomp,form='unformatted')
			write(iunit) sngl(zs),sngl(ze),sngl(sigma),nf1,nf2,sngl(df),sngl(dwn),sngl(earthmodel_rearth),
     1		             jump_spher(m),compssph(i),sourcetype,nwnmax
			if(selflag.eq.1) then
				do j=1,npar
					iunit=iunit+1
					write(cpar,'(i1)') pars(j)
					open(iunit,file=partials(1:icp)//'.sph'//'.s'//sourcetype//'.j'//cjump//
     1				     '.c'//ckomp//'.p'//cpar,form='unformatted')
					write(iunit) sngl(zs),sngl(ze),sngl(sigma),nf1,nf2,sngl(df),
     1					       sngl(dwn),jump_spher(m),compssph(i),pars(j),sourcetype,nwnmax
					write(iunit) earthmodel_nk,(sngl((earthmodel_rearth-earthmodel_rk(k))*1.d3),
     1				                                k=1,earthmodel_nk)
				enddo
			endif
		enddo
	enddo
	endif
c
c  toroidal results
c
	if(nsrctor.gt.0) then
	do i=1,ncomptor
		write(ckomp,'(i1)') compstor(i)
		do m=1,nsrctor
			iunit=6+ncompsph*nsrcsph*(npar+1)+(i-1)*nsrctor*(npartor+1)+(m-1)*(npartor+1)+1
			write(cjump,'(i1)') jump_tor(m)	
			open(iunit,file=greenfile(1:icg)//'.tor'//'.s'//sourcetype//'.j'//cjump//
     1		     '.c'//ckomp,form='unformatted')
			write(iunit) sngl(zs),sngl(ze),sngl(sigma),nf1,nf2,sngl(df),sngl(dwn),sngl(earthmodel_rearth),
     1		         jump_tor(m),compstor(i),sourcetype,nwnmax
			if(selflag.eq.1) then
				do j=1,npartor
					iunit=iunit+1
					write(cpar,'(i1)') parstor(j)
					open(iunit,file=partials(1:icp)//'.tor'//'.s'//sourcetype//'.j'//cjump//
     1					     '.c'//ckomp//'.p'//cpar,form='unformatted')
					write(iunit) sngl(zs),sngl(ze),sngl(sigma),nf1,nf2,sngl(df),
     1				             sngl(dwn),jump_tor(m),compstor(i),parstor(j),sourcetype,nwnmax
					write(iunit) earthmodel_nk,(sngl((earthmodel_rearth-earthmodel_rk(k))*1.d3),
     1			                                k=1,earthmodel_nk)
				enddo
			endif
		enddo
	enddo
	endif
c
c  Frequency-slowness loop
c
	fmax=(nf2-1)*df
	fmin=(nf1-1)*df
	do if=nf1,nf2
		f=(if-1)*df
		om=2.*pi*f
c
c  update model parameters to account for dispersion (seldamp=3)
c  if seldamp in the model file is set to 1, dispersion is suppressed
c
		if(earthmodel_seldamp.eq.3) then
			call earthmodel_fromflnm(f,if-nf1)
			call earthmodel_spline
		endif
c
		if(iprint.gt.0) then
			print *,'-------------------------------------------------'
			print *,'           FLGEVASK: f = ',f
			print *,'-------------------------------------------------'
		endif
c
c  compute source terms
c  zeps is needed for output in case of moment tensor sources
c
		call sterms(sourcetype,zs,zsrcsph,zsrctor)
c
c  Loop over wavenumbers
c  Green functions and partial derivatives
c
		wn2=om*pmax+frac*wnmax*(fmax-f)/(fmax-fmin)
c		wn2=om*pmax*(1.+frac*(fmax-f)/(fmax-fmin))
		nwn=int(wn2/dwn)+1
		do iwn=1,nwn
			if(iprint.gt.0) then
				print *,'--------------------------------------------------------'
				print *,'        IWN = ',iwn
				print *,'--------------------------------------------------------'
			endif
			wn=(iwn-1)*dwn
			if(ncompsph.gt.0.and.nsrcsph.gt.0) then
				call greenkfsph(wn,f,sigma,zs,ze,ncompsph,compssph,npar,pars,
     1				          gfsrsph(1,1,iwn),sksph(1,1,1,1,iwn),selflag,
     1				          nsrcsph,zsrcsph,jump_spher,iprint)
			endif
			if(ncomptor.gt.0.and.nsrctor.gt.0) then
				call greenkftor(wn,f,sigma,zs,ze,ncomptor,compstor,npartor,parstor,
     1				          gfsrtor(1,1,iwn),sktor(1,1,1,1,iwn),selflag,
     1				          nsrctor,zsrctor,jump_tor,iprint)
			endif
		enddo
c
c  write spheroidal Green functions
c
		if(nsrcsph.gt.0) then
		do i=1,ncompsph
			do m=1,nsrcsph
				mm=jump_spher(m)
				iunit=6+(i-1)*nsrcsph*(npar+1)+(m-1)*(npar+1)+1
				write(iunit) nwn,(cmplx(gfsrsph(compssph(i),mm,iwn)),iwn=1,nwn)
				if(selflag.eq.1) then
					do j=1,npar
						iunit=iunit+1
						write(iunit) nwn
						do k=1,earthmodel_nk
							write(iunit) (cmplx(sksph(k,compssph(i),pars(j),mm,iwn)),iwn=1,nwn)
						enddo
					enddo
				endif
			enddo
		enddo
		endif
c
c  write toroidal Green functions
c
		if(nsrctor.gt.0) then
		do i=1,ncomptor
			do m=1,nsrctor
				mm=jump_tor(m)
				iunit=6+ncompsph*nsrcsph*(npar+1)+(i-1)*nsrctor*(npartor+1)+(m-1)*(npartor+1)+1
				write(iunit) nwn,(cmplx(gfsrtor(compstor(i),mm,iwn)),iwn=1,nwn)
				if(selflag.eq.1) then
					do j=1,npartor
						iunit=iunit+1
						write(iunit) nwn
						do k=1,earthmodel_nk
							write(iunit) (cmplx(sktor(k,compstor(i),parstor(j),mm,iwn)),iwn=1,nwn)
						enddo
					enddo
				endif
			enddo
		enddo
		endif
		if(if.lt.nf2) write(6,'(i6,i6,$)') if,nwn
		if(if.eq.nf2) write(6,'(i3,i6)') if,nwn
	enddo
c
c  close files
c
	do i=1,ncompsph*nsrcsph*(npar+1)+ncomptor*nsrctor*(npartor+1)
		close(6+i)
	enddo
c
	stop
	end
c------------------------------------------
c  printhelp
c-------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                     FLGEVASK'
	print *,''
	print *,'Usage: flgevask [-h] [-p] [-tlen] [-sig] [-fmin] [-fmax] [-pmax]'
	print *,'             [-zs] [-ze] [-od] [-op] [-xlen]'
	print *,'             [-nmin] [-sf] [-c] [-s] [-m] [-fnode] [-wmrg] earthmodel'
	print *,''
	print *,'This code is just a driver for greenkf and pardev. See doku there.'
	print *,''
	print *,'Options:'
	print *,'-h     : print help'
	print *,'-p     : printlevel (def = 0)'
	print *,'-tlen  : time series length in seconds (def = 1.) '
	print *,'-fmin  : minimum frequency in Hz (def = 5), if set to zero, it will be adjusted to 1/tlen'
	print *,'-fmax  : maximum frequency in Hz (def = 40)'
	print *,'-pmax  : maximum slowness in s/km (def = 5 s/km). Used to determine'
	print *,'         the wavenumber range.'
	print *,'-sig   : ratio of tlen and damping time for complex frequency (def = 5.)'
	print *,'-zs    : source depth in m (def = 1 m)'
	print *,'-ze    : receiver depth in m (def = 0 m)'
	print *,'-od    : basename of output file with gfs(comp,jump,k,f,re,rs) (green).'
	print *,'         Filename is composed as follows:'
	print *,'         basename.sph.s[f|m].j[jump].c[comp]  for spheroidal solutios'
	print *,'         basename.tor.s[f|m].j[jump].c[comp]  for toroidal solutions.'
	print *,'         [f|m] is either force or moment tensor source,'
	print *,'         [jump] takes values from 1 to 4 for spheroidal solutions and'
	print *,'         [jump] takes values from 1 to 2 for toroidal solutions.'
	print *,'         [comp] takes values from 1 to 4 (U,R,V,S) for spheroidal solutions'
	print *,'         [comp] takes values from 1 to 2 (W,T) for toroidal solutions'
	print *,'         see definition of jump vectors in gemini-displ.pdf'
	print *,'-op    : basename of output file with partials sk(nr,komp,par,jump,k,f)'
	print *,'         (def = partials). See filename conventions under option -od'
	print *,'         Index of parameter is appended to the filename.'
	print *,'-xlen  : length of profile in km (def = 1 km) '
	print *,'         (determines wavenumber stepping)'
	print *,'-nmin  : minimum number of solution nodes per layer (def=5)'
	print *,'         This does not affect the accuracy of the solution at the receiver'
	print *,'         but it affects the accuracy of the partials' 
	print *,'-sf    : select flag (def=1):'
	print *,'         (1) compute displacement-stress vector and Frechet derivatives'
	print *,'         (2) compute only displacement-stress vector (saves time)'
	print *,'-c     : a blank-separated list of displacement-stress vector'
	print *,'         components for which Green functions are desired (def = 1)'
	print *,'         (1) vertical displacement (U)'
	print *,'         (2) rr-stress ( = negative pressure in fluid regions) (R)'
	print *,'         (3) horizontal diplacement (V)'
	print *,'         (4) rx-stress (S)'
	print *,'         (5) horizontal displacement (W)'
	print *,'         (6) ry-stress (T)'
	print *,'-s     : type of source.'
	print *,'         (1) explosion'
	print *,'         (2) vertical force'
	print *,'         (3) horizontal force'
	print *,'         (4) general single force'
	print *,'         (5) moment tensor'
	print *,'-m     : a blank-separated list of model parameter indices for which'
	print *,'         Frechet kernels are desired (def = 3).'
	print *,'         (1) density'
	print *,'         (2) P-velocity'
	print *,'         (3) S-velocity'
	print *,'         (4) Q-kappa'
	print *,'         (5) Q-mue'
	print *,'-fnod  : frequency to determine node spacing in Hz (f=5)'
	print *,'-wmrg  : margin of wavenumbers beyond om*pmax at lowest frequency as'
	print *,'         a fraction of the maximum wavenumber (def=0.1)'
	print '(50(1h-))'
	return
	end
c-----------------------------------------------------------------------
c  print dimensions to standard output
c-----------------------------------------------------------------------
	subroutine print_dimensions
	include 'earthmodel.h'
	include 'flnm.h'
	include 'greendim.h'
	include 'nodesdim.h'
	include 'layer.h'
	include 'nvmax.h'
	include 'bsstore.h'
c
	print '(50(1h-))'
	print *,'The program was compiled with the following dimensions:'
	print *,''
	print *,'earthmodel.h: nkk     = ',nkk
	print *,'flnm.h:       nkk1    = ',nkk1
	print *,'greendim.h:   nff     = ',nff
	print *,'greendim.h:   npp     = ',npp
	print *,'nodesdim.h:   nnd     = ',nnd
	print *,'layer.h:      nlayer  = ',nlayer
	print *,'nvmax.h:      nvmax   = ',nvmax
	print *,'bsstore.h:    nnstore = ',nnstore
	print '(50(1h-))'
c
	return
	end
