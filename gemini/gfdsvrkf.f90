!---------------------------------------------------------------
 	program gfdsvrkf
	use commandLine
	use mathConstants
	use streamAccess
	use primitiveTypeEncoding
	implicit none
!
	include 'earthmodel.h'
	include 'flnm.h'
	include 'nodesdim.h'
	include 'nodes.h'
	type (file_stream_access) :: fda_sph,fda_tor
	type (group_stream_access) :: root_sph,root_tor,frequency,source_radius
	type (data_stream_access) :: header_real_sph,header_int_sph,rnodes_sph,dset
	type (data_stream_access) :: header_real_tor,header_int_tor,rnodes_tor
	double complex, dimension(:,:,:,:), allocatable :: gfesph,gfetor
	double complex, dimension(4,4,nnd) :: zsph
	double complex, dimension(2,2,nnd) :: ztor
	double precision :: tlen,xlen,fmin,fmax,f,df,sigma,wn,dwn,fnode,om
	double precision :: p1lim,p2lim
	double precision :: zstore,zs,ze,wn2,wnmax,frac,r1,r2
	real, dimension(:), allocatable :: d
	real :: max_source_depth
	integer, dimension(:), allocatable :: id
	integer, dimension(6) :: dsvmask_sph
	integer, dimension(3) :: dsvmask_tor
	integer :: i,m,if,ifh,nf1,nf2,iwn,nwn,ice,k,nsrcsph,nsrctor,nlay,nr
	integer :: jrstep,nrdsvsel,dsvstep,derivflag
	integer :: iprint,nnmin,istyp,nwnmax,lusph,lutor,ierr,jsmin,nrdsv,nkfmax,korg
	integer, dimension(:), allocatable :: dims
	character (len=132) :: string,cmodel,dsvbase,nodefile
	character (len=1) :: sourcetype
	logical :: rstflag,vsflag,alternative
	type (cmdLine) :: cl
	external printhelp
!-----------------------------------------------------------
!  defaults for input parameters
!
	data iprint/0/tlen/1./sigma/0.d0/fmin/5./fmax/40./derivflag/0/
     1     p1lim/0.5/p2lim/0.25/dsvbase/'gfdsvrfk'/rstflag/.false./vsflag/.false./
     2     xlen/1./nnmin/5/istyp/1/ze/0./max_source_depth/1000000./
     4     fnode/5./frac/0.1/nodefile/'nodes.out'/jrstep/1/alternative/.false./
!---------------------------------------------------------------------
!  print dimensions to standard output
!
	call print_dimensions
!---------------------------------------------------------------------
!  read commandline
!
	call new(cl,20,1,'h a d p tlen fmin fmax klim sig ze s o xlen nmin fnod wmrg on msd rst vs',
     1	         '0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0',printhelp)
	cmodel = clManarg(cl,1)
	if (clOptset(cl,'a')) then; alternative = .true.; endif
	if (clOptset(cl,'d')) then; derivflag = 1; endif
	if (clOptset(cl,'rst')) then; rstflag = .true.; endif
	if (clOptset(cl,'p')) then; string = clOptarg(cl,'p'); read(string,*) iprint; endif
	if (clOptset(cl,'tlen')) then; string = clOptarg(cl,'tlen'); read(string,*) tlen; endif
	if (clOptset(cl,'fmin')) then; string = clOptarg(cl,'fmin'); read(string,*) fmin; endif
	if (clOptset(cl,'fmax')) then; string = clOptarg(cl,'fmax'); read(string,*) fmax; endif
	if (clOptset(cl,'klim')) then; string = clOptarg(cl,'klim'); read(string,*) p1lim,p2lim; endif
	if (clOptset(cl,'sig')) then; string = clOptarg(cl,'sig'); read(string,*) sigma; endif
	if (clOptset(cl,'ze')) then; string = clOptarg(cl,'ze'); read(string,*) ze; endif
	if (clOptset(cl,'s')) then; string = clOptarg(cl,'s'); read(string,*) istyp; endif
	if (clOptset(cl,'o')) then; dsvbase = clOptarg(cl,'o'); endif
	if (clOptset(cl,'xlen')) then; string = clOptarg(cl,'xlen'); read(string,*) xlen; endif
	if (clOptset(cl,'nmin')) then; string = clOptarg(cl,'nmin'); read(string,*) nnmin; endif
	if (clOptset(cl,'fnod')) then; string = clOptarg(cl,'fnod'); read(string,*) fnode; endif
	if (clOptset(cl,'wmrg')) then; string = clOptarg(cl,'wmrg'); read(string,*) frac; endif
	if (clOptset(cl,'on')) then; nodefile = clOptarg(cl,'on'); endif
	if (clOptset(cl,'msd')) then; string = clOptarg(cl,'msd'); read(string,*) max_source_depth,jrstep; endif
	if (clOptset(cl,'vs')) then; vsflag = .true.; verboseStreamAccess = .true.; endif
	call clPrint(cl)
	call dealloc(cl)
!--------------------------------------------------------------
!  go
!
	sigma=sigma/tlen
	df=1.d0/tlen
	if(fmin.lt.df) fmin=df
	nf1=nint(fmin/df)+1
	nf2=nint(fmax/df)+1
	fmin=(nf1-1)*df
	fmax=(nf2-1)*df
	dwn=2.*mc_pi/xlen
	nnmin=max(nnmin,2)
!
!  Initialize
!
	call flnm_read(1,cmodel)
!
!  Note: the earth model is for the reference frequency. If dispersion
!        is desired, elastic parameters should be recomputed for
!        each frequency. This requires an additional routine.
!
	call earthmodel_fromflnm(flnm_fref,0)
	call earthmodel_layer
	call layer_printinfo
	call earthmodel_spline
	call layer_set_stepsize_wlfrac(0.25d0)
!
!  zstore is set here to depth of halfspace+1/4 of the thickness of the overburden.
!  The code assumes that zstore is in the halfspace !!! (nlstore=1)
!
	call layer_getrb(1,r1)
	call layer_getnlay(nlay)
	call layer_getrb(nlay,r2)
	zstore=(r2-r1)*1.25
!
!  if earth model has a water layer on top and ze is above ocean bottom, 
!  receiver is set to just below the ocean bottom
!  ze converted to km
!
	call layer_isfluid(nlay,i)
	call layer_getrb(nlay-1,r1)
	if (i.eq.1 .and. (r2-r1).gt.ze*1.d-3) then
		ze = r2-r1+1.d-6
	else
		ze = max(0.001d0,ze)*1.d-3
	endif
	if (alternative) then
		print *,'Source at depth: ',ze,' km'
	else
		print *,'Receiver at depth: ',ze,' km'
	endif
!
!  dummy for nodes_create, not used
!
	zs=ze
!
!  create DSV nodes
!
	call nodes_create(fnode,4.d0,zstore,zs,ze,nnmin,iprint)
	call nodes_print(6,'hallo')
	call nodes_print(2,nodefile)
!
!  set number of required jump solutions according to excitation
!
	call trimright(dsvbase,ice)
	if(istyp.eq.0) then
		sourcetype='f'
		nsrcsph=2
		nsrctor=1
	else if(istyp.eq.1) then
		sourcetype='m'
		nsrcsph=4
		nsrctor=2
	else
		print *,'Unknown source type: ',istyp
		stop
	endif
!
!  set masks of required DSV-spectra: Spheroidal: (U,R,V,S,UP,VP), Toroidal: (W,T,WP)
!
	if (rstflag) then; dsvstep = 1; else; dsvstep = 2; endif
	dsvmask_sph = (/ 1,0,1,0,0,0 /)
	dsvmask_tor = (/ 1,0,0 /)
	if (derivflag == 1) then
		dsvmask_sph(5:6) = 1
		dsvmask_tor(3) = 1
	endif
	if (dsvstep == 1) then
		dsvmask_sph(2:4:2) = 1
		dsvmask_tor(2) = 1
	endif
!
!  wavenumber limits
!
	wnmax = 2.*mc_pi*0.5*fmax*(p1lim+p2lim)/(1.-frac)
	nwnmax=int(wnmax/dwn)+1
	allocate(gfesph(nnd,6,4,nwnmax),gfetor(nnd,3,2,nwnmax))
	print *,'Maximum amount of wavenumbers: ',nwnmax
!
!  find DSV node index of maximum source OR receiver depth == minimum DSV radius
!
	call dlocate(earthmodel_rearth-max_source_depth*1.d-3,nnod,rnod,jsmin)
	if (jsmin == 0) jsmin = 1
	nrdsv = min(nnod-jsmin+1,nnod)
	nrdsvsel = (nrdsv-1)/jrstep+1         ! less nodes due to stepping
	if (alternative) then
		print *,'Receivers down to :',max_source_depth*1.e-3
		print *,'Number of receiver nodes in output: ',nrdsvsel
		print *,'Index of new bottom receiver node: ',jsmin
	else
		print *,'Sources down to :',max_source_depth*1.e-3
		print *,'Number of source nodes in output: ',nrdsvsel
		print *,'Index of new bottom source node: ',jsmin
	endif
!
!  open spheroidal dsvfile for stream access
!
	lusph=1
	ierr = createFileStreamAccess(fda_sph,lusph,dsvbase(1:ice)//'.s'//sourcetype//'.sph')    ! create file
	if (ierr /= 0) then; print *,'File open error'; stop; endif
!
!  open toroidal dsvfile for stream access
!
	lutor = 2
	ierr = createFileStreamAccess(fda_tor,lutor,dsvbase(1:ice)//'.s'//sourcetype//'.tor')    ! create file
	if (ierr /= 0) then; print *,'File open error'; stop; endif
!
!  root groups
!
	call createGroupStreamAccess(root_sph,'Root',0)                        ! create root group
	call createGroupStreamAccess(root_tor,'Root',0)                        ! create root group
!
!  header data for both files
!  write spheroidal and toroidal header data
!  integer header data written at end
!
	allocate(dims(1))
	dims(1) = 5
	allocate(d(5))
	d = (/ sngl(ze),sngl(sigma),sngl(df),sngl(dwn),sngl(earthmodel_rearth) /)
	call createDatasetStreamAccess(header_real_sph,1,dims,T_REAL)  ! real header data in root
	call createDatasetStreamAccess(header_real_tor,1,dims,T_REAL)  ! real header data in root
	call writeDatasetVectorStreamAccess(header_real_sph,fda_sph,d)
	call writeDatasetVectorStreamAccess(header_real_tor,fda_tor,d)
	call addDatasetStreamAccess(root_sph,header_real_sph)
	call addDatasetStreamAccess(root_tor,header_real_tor)
	deallocate(d)
!
!  values of node radii
!  write node radii to files
!
	dims(1) = nrdsvsel
	call createDatasetStreamAccess(rnodes_sph,1,dims,T_DOUBLE)    ! double data set with node radii
	call createDatasetStreamAccess(rnodes_tor,1,dims,T_DOUBLE)    ! double data set with node radii
	call writeDatasetVectorStreamAccess(rnodes_sph,fda_sph,rnod(jsmin:nnod:jrstep))
	call writeDatasetVectorStreamAccess(rnodes_tor,fda_tor,rnod(jsmin:nnod:jrstep))
	call addDatasetStreamAccess(root_sph,rnodes_sph)
	call addDatasetStreamAccess(root_tor,rnodes_tor)
!
!  Frequency-wavenumber loop
!
	fmax=(nf2-1)*df
	fmin=(nf1-1)*df
	nkfmax=0
	do if=nf1,nf2
		ifh = if-nf1+1
		f=(if-1)*df
		om=2.*mc_pi*f
!
!  update model parameters to account for dispersion (seldamp=3)
!  if seldamp in the model file is set to 1, dispersion is suppressed
!
		if(earthmodel_seldamp.eq.3) then
			call earthmodel_fromflnm(f,1)
			call earthmodel_spline
		endif
!
		if(iprint.gt.0) then
			print *,'-------------------------------------------------'
			print *,'           GFDSVRKF: f = ',f
			print *,'-------------------------------------------------'
		endif
!
!  compute jump vectors for each radial node
!  OR compute jump vector at ze only
!
		if (alternative) then
			call sterms(sourcetype,ze,zsph(1,1,1),ztor(1,1,1))       ! jump vector needed at one depth only
		else
			do nr=1,nnod
				call sterms(sourcetype,earthmodel_rearth-rnod(nr),zsph(1,1,nr),ztor(1,1,nr))
			enddo
		endif
!
!  Loop over wavenumbers
!  Green functions gfe(rnod,dsv,jump,wn)
!
		if (f <= 0.5*fmax) then
			wn2 = om*p1lim+frac*wnmax
		else
			wn2 = 2.*mc_pi*0.5*fmax*(p1lim+p2lim)+2.*mc_pi*(f-fmax)*p2lim+frac*wnmax
		endif
		nwn=int(wn2/dwn)+1
		do iwn=1,nwn
			if(iprint.gt.0) then
				print *,'--------------------------------------------------------'
				print *,'        IWN = ',iwn
				print *,'--------------------------------------------------------'
			endif
			wn=(iwn-1)*dwn
			call jump_solution_sph(wn,f,sigma,ze,gfesph(1,1,1,iwn),nsrcsph,zsph,iprint,alternative,derivflag)
			call jump_solution_tor(wn,f,sigma,ze,gfetor(1,1,1,iwn),nsrctor,ztor,iprint,alternative,derivflag)
		enddo
!
!  prepare writing
!
		dims(1) = nwn
		nkfmax=nkfmax+nwn    ! count total amount of wavenumbers
!
!  write spheroidal DSV at receiver for sources at nodes (gfesph)
!
		call createGroupStreamAccess(frequency,'Frequency',if)
		do k=1,nrdsvsel
			call createGroupStreamAccess(source_radius,'Radius',k)
			korg = 1+(k-1)*jrstep
			do i=1,6
				if (dsvmask_sph(i) == 0) cycle
				do m=1,nsrcsph
					call createDatasetStreamAccess(dset,1,dims,T_COMPLEX)
					call writeDatasetVectorStreamAccess(dset,fda_sph,cmplx(gfesph(korg+jsmin-1,i,m,1:nwn)))
					call addDatasetStreamAccess(source_radius,dset)
					call dealloc(dset)
				enddo
			enddo
			call addSubgroupStreamAccess(frequency,source_radius)
			call dealloc(source_radius)
		enddo
		call addSubgroupStreamAccess(root_sph,frequency)
		call dealloc(frequency)
!
!  write toroidal DSV at nodes from source (gfstor) and from receiver (gfetor)
!
		call createGroupStreamAccess(frequency,'Frequency',if)
		do k=1,nrdsvsel
			call createGroupStreamAccess(source_radius,'Radius',k)
			korg = 1+(k-1)*jrstep
			do i=1,3
				if (dsvmask_tor(i) == 0) cycle
				do m=1,nsrctor
					call createDatasetStreamAccess(dset,1,dims,T_COMPLEX)
					call writeDatasetVectorStreamAccess(dset,fda_tor,cmplx(gfetor(korg+jsmin-1,i,m,1:nwn)))
					call addDatasetStreamAccess(source_radius,dset)
					call dealloc(dset)
				enddo
			enddo
			call addSubgroupStreamAccess(frequency,source_radius)
			call dealloc(source_radius)
		enddo
		call addSubgroupStreamAccess(root_tor,frequency)
		call dealloc(frequency)
!
		if(if.lt.nf2) write(6,'(i6,i6,$)') if,nwn
		if(if.eq.nf2) write(6,'(i3,i6)') if,nwn
	enddo
!
!  write remaining integer header data
!
	allocate(id(9))
	dims(1) = 9
	id = (/ istyp,nf1,nf2,nkfmax,nwnmax,nrdsvsel,nsrcsph,dsvstep,derivflag /)
	call createDatasetStreamAccess(header_int_sph,1,dims,T_INTEGER)  ! integer header data in root
	call writeDatasetVectorStreamAccess(header_int_sph,fda_sph,id)
	call addDatasetStreamAccess(root_sph,header_int_sph)
	id = (/ istyp,nf1,nf2,nkfmax,nwnmax,nrdsvsel,nsrctor,dsvstep,derivflag /)
	call createDatasetStreamAccess(header_int_tor,1,dims,T_INTEGER)  ! integer header data in root
	call writeDatasetVectorStreamAccess(header_int_tor,fda_tor,id)
	call addDatasetStreamAccess(root_tor,header_int_tor)
	deallocate(id)
!
!  write group tree
!
	call writeGroupStreamAccess(root_sph,fda_sph)
	call writeGroupStreamAccess(root_tor,fda_tor)
!
!  deallocation
!
	call clearGroupTree(root_sph); call clearGroupTree(root_tor)
	call dealloc(fda_sph); call dealloc(fda_tor)
	deallocate(dims,gfesph,gfetor)
!
	end program
!------------------------------------------
!  printhelp
!-------------------------------------------
	subroutine printhelp
	print '(50(1h-))'
	print *,'                     GFDSVRKF'
	print *,''
	print *,'Usage: gfdsvrkf [-h] [-p] [-a] [-d] [-tlen] [-sig] [-fmin] [-fmax] [-klim] [-ze] '
	print *,'                [-o] [-xlen] [-nmin] [-s] [-fnod] [-wmrg] [-on] [-msd] earthmodel'
	print *,''
	print *,'Compute displacement-stress fk-spectra at a fixed receiver depth ze'
	print *,'for single force or moment tensor sources at all depth nodes.'
	print *,'OR alternatively,'
	print *,'compute displacement-stress fk-spectra and radial derivatives at receivers at all depth nodes'
	print *,'for a fixed single force or moment tensor source at depth ze.'
	print *,''
	print *,'Options:'
	print *,''
	print *,'-h     : print help'
	print *,'-a     : do the alternative calculation (spectra for one source at many receiver depths)'
	print *,'         Default is spectra at one receiver for many source depths.'
	print *,'-d     : compute also radial derivatives of displacment depctra (def = no)' 
	print *,'-rst   : also write fk-spectra of stress DSV: (def = no)'
	print *,'-p     : printlevel (def = 0)'
	print *,'-tlen  : time series length in seconds (def = 1.) '
	print *,'-fmin  : minimum frequency in Hz (def = 5)'
	print *,'-fmax  : maximum frequency in Hz (def = 40)'
	print *,'-klim  : Define upper limit of wavenumber versus frequency by two possibly different slopes p1 and p2'
	print *,'         First slope is valid up to the 0.5*fmax, second one beyond'
	print *,'-sig   : ratio of tlen and damping time for complex frequency (def = 5.)'
	print *,'-ze    : receiver depth in m (def = 0 m) OR source depth in m'
	print *,'-o     : basename of output file with gfdsv(rnod,k,f,comp,jump).'
	print *,'         Filename is composed as follows:'
	print *,'         basename.s[f|m].sph  for spheroidal solutios'
	print *,'         basename.s[f|m].tor  for toroidal solutions.'
	print *,'         [f|m] is either force or moment tensor source,'
	print *,'-xlen  : length of profile in km (def = 1 km) '
	print *,'         (determines wavenumber stepping)'
	print *,'-nmin  : minimum number of solution nodes per layer (def=5)'
	print *,'-s     : type of source.'
	print *,'         (0) general single force'
	print *,'         (1) moment tensor'
	print *,'-fnod  : frequency used for node spacing in Hz (5.)'
	print *,'-on    : name of output file for node information'
	print *,'-wmrg  : margin of wavenumbers beyond om*pmax at lowest frequency as'
	print *,'         a fraction of the maximum wavenumber (def=0.1)'
	print *,'-msd   : maximum source OR receiver depth in meter (def = bottom radial node)'
	print *,'         and stepping regarding radial nodes (def = 1)'
	print *,'-vs    : streamAccess file written with many comments for debugging'
	print '(50(1h-))'
	return
	end
!-----------------------------------------------------------------------
!  print dimensions to standard output
!-----------------------------------------------------------------------
	subroutine print_dimensions
	include 'earthmodel.h'
	include 'flnm.h'
	include 'nodesdim.h'
	include 'layer.h'
	include 'nvmax.h'
!
	print '(50(1h-))'
	print *,'The program was compiled with the following dimensions:'
	print *,''
	print *,'earthmodel.h: nkk     = ',nkk
	print *,'flnm.h:       nkk1    = ',nkk1
	print *,'nodesdim.h:   nnd     = ',nnd
	print *,'layer.h:      nlayer  = ',nlayer
	print *,'nvmax.h:      nvmax   = ',nvmax
	print '(50(1h-))'
!
	return
	end
