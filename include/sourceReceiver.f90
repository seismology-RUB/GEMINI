!---------------------------------------------
!  source and receiver information module (SAR)
!
!        subroutine sar_readinfo(this,iunit,infofile,epimax)
!        subroutine deallocSourceReceiver(this)
!        subroutine sar_read_rcv_bin(this,iunit,filename,epimax)
!        subroutine sar_epico(this,x,y,dis,phi)
!        subroutine sar_propdir(this,x,y,azi)
!        subroutine sar_getstaindex(this,selsta,k)
!        subroutine sar_swapstations(this,j,k)
!        subroutine sar_setnr(this,n)
!        subroutine sar_setComps(this,comp)
!        subroutine sar_rad2km(this,rearth)
!        subroutine sar_eventhms(this,hour,minute,sec)
!        subroutine sar_scalar_moment(this,rmsc)
!        subroutine sar_tm_epicar_globalcar(this)
!        subroutine sar_tm_epispher_globalcar(this,delta,azi,tm)
!        subroutine sar_geogr2geoce(this,cflag,j)
!        subroutine sar_sortPhi(this,indx,phisorted)
!        function   sar_getNumberTraces(this,nr)
!        function   sar_getEventid(this,eventid)
!        function   sar_CoordinateSystem(this)
!        function   sar_getDistanceSelected(this,ksta)
!        function   sar_getPhiSelected(this,ksta)
!        function   sar_getPropdirSelected(this,ksta)
!        function   sar_getNetIdSelected(this,ksta)
!        function   sar_getCompSelected(this,ksta)
!	  function   sar_getStationSelected(this,ksta)
!        function   sar_getXSelected(this,ksta)
!        function   sar_getYSelected(this,ksta)
!        function   sar_getStationDepth(this,ksta)
!        function   sar_getXSource(this)
!        function   sar_getYSource(this)
!        function   sar_getDateSource(this)
!        function   sar_getTimeSource(this)
!        function   sar_getTSAMSource(this)
!        function   sar_getYearSource(this)
!        function   sar_getMonthSource(this)
!        function   sar_getDaySource(this)
!        function   sar_getNetFlag(this)
!        function   sar_getSourceDepth(this)
!        function   sar_getSourceType(this)
!	  function   sar_getDistances(this)
!	  function   sar_getAzimuths(this)
!        subroutine sar_getMomentTensor(this,rm)
!        subroutine sar_getForce(this,f)
!-----------------------------------------------------------------------
!
!  cartesian coordinates: 
!	rxy with phi from x to y
!	forces: 1=fr, 2=fx, 3=fy
!	moments: 1=Mrr, 2=Mxx, 3=Myy, 4=Mrx, 5=Mry, 6=Mxy
!	displacement components: either RLT with T pointing left if one looks from source to receiver
!	                         or RXY
!	sar_depsrc = source depth in meters
!	sar_depsta = station depth in meters
!	sar_xsrc = x-coordinate of source in meters
!	sar_ysrc = y-coordinate of source in meters
!	sar_xsta = x-coordinate of station in meters
!	sar_ysta = y-coordinate of station in meters
!	sar_dis(i) in kilometers
!	sar_phi(i) in radians
!
!  spherical coordinates: 
!	rtf with phi and from south over east
!	forces: 1=fr, 2=ft, 3=ff using basis vectors at the source point!
!	moments: 1=Mrr, 2=Mtt, 3=Mff, 4=Mrt, 5=Mrf, 6=Mff using basis vectors at the source point!
!	displacement components: either RLT with T pointing left if one looks from source to receiver
!	                         or RNE
!
!	sar_depsrc = source depth in km
!	sar_depsta = station depth in km
!	sar_xsrc = Colatitide of source in degrees
!	sar_ysrc = Longitude of source in degrees
!	sar_xsta = Colatitude of station in degrees
!	syr_ysta = Longitude of station in degrees
!	sar_dis(i) in radians
!	sar_phi(i) in radians
!
!	netflag = 0: ID STA COLAT LON
!	netflag = 1: ID STA NET COLAT LON
!	netflag = 2: ID STA COMP NET COLAT LON
!	netflag = 3: ID STA COMP COLAT LON
!
!	gg2gcflag: convert from geographic to geocentric coordinates
!	           = 0:  do not transform at all
!	           = 1:  transform both source and receiver coordinates
!	           = 2:  transform receiver coordinates only
!
 module sourceReceiver
	use mathConstants
	use dayOfYear
	use realloc
	implicit none
	interface new
		module procedure sar_readinfo
	end interface
	interface dealloc
		module procedure deallocSourceReceiver
	end interface
	interface operator (.dis.);     module procedure sar_getDistanceSelected; end interface
	interface operator (.phi.);     module procedure sar_getPhiSelected;      end interface
	interface operator (.propdir.); module procedure sar_getPropdirSelected;  end interface
	interface operator (.netid.);   module procedure sar_getNetIdSelected;    end interface
	interface operator (.comp.);    module procedure sar_getCompSelected;     end interface
	interface operator (.station.); module procedure sar_getStationSelected;  end interface
	interface operator (.csys.);    module procedure sar_CoordinateSystem;    end interface
	interface operator (.x.);       module procedure sar_getXSelected;        end interface
	interface operator (.y.);       module procedure sar_getYSelected;        end interface
	interface operator (.ze.);      module procedure sar_getStationDepth;     end interface
	interface operator (.xs.);      module procedure sar_getXSource;          end interface
	interface operator (.ys.);      module procedure sar_getYSource;          end interface
	interface operator (.zs.);      module procedure sar_getSourceDepth;      end interface
	interface operator (.nr.);      module procedure sar_getNumberTraces;     end interface
	interface operator (.evid.);    module procedure sar_getEventId;          end interface
	interface operator (.sdate.);   module procedure sar_getDateSource;       end interface
	interface operator (.sdatefy.); module procedure sar_getFullYearDateSource; end interface
	interface operator (.stime.);   module procedure sar_getTimeSource;       end interface
	interface operator (.therd.);   module procedure sar_getTSAMSource;       end interface
	interface operator (.syear.);   module procedure sar_getYearSource;       end interface
	interface operator (.smonth.);  module procedure sar_getMonthSource;      end interface
	interface operator (.sday.);    module procedure sar_getDaySource;        end interface
	interface operator (.sdoy.);    module procedure sar_getDayOfYearSource;  end interface
	interface operator (.styp.);    module procedure sar_getSourceType;       end interface
	interface operator (.tfs.);     module procedure sar_getSourceTimeFullSeconds;       end interface
	interface operator (.tns.);     module procedure sar_getSourceTimeNanoSeconds;       end interface
	type sar
		private
!--------------------------------------------------------------------------------------------
		character (len=1) :: csys         ! coordinate system (C/S)
!--------------------------------------------------------------------------------------------
!  source part
!
		integer :: year                   !  time information about source
		integer :: month
		integer :: day
		integer :: hour
		integer :: minute
		real :: sec
		real :: cendt                      !  centroid time shift
		real :: xsrc                       !  x [m] / colatitude [deg] of source
		real :: ysrc                       !  y [m] / longitude  [deg] of source
		real :: depsrc                     !  depth of source (C:meters/S:converted to km)
		character (len=13) :: eventid      ! yymmdd_hhmmss
		integer :: istyp                   ! type of source (0: force, 1: moment tensor)
		real, dimension(6) :: mt           ! moment tensor
		real, dimension(3) :: force        ! force components
		real :: hdur                       ! half duration
		real, dimension(3,3) :: tm_ec_gc   ! transformation from epicentral cartesian to global cartesian
!--------------------------------------------------------------------------------------------------------------
!  receiver part
!
		integer :: netflag                 ! flag distinguishing info file types
		integer :: nr                      ! number of seismograms used
		integer :: nrorg                   ! original number of seismograms in info file
		integer, dimension(:), pointer :: lfdnr => null()  ! original number of seismogram
		character (len=5), dimension(:), pointer :: name => null()
		character (len=1), dimension(:), pointer :: comp => null()
		character (len=2), dimension(:), pointer :: net => null()
		real, dimension(:), pointer :: xsta => null()    ! x [m]/ colatitude [deg] of receiver
		real, dimension(:), pointer :: ysta => null()    ! y [m]/ longitude [deg] of receiver
		real, dimension(:), pointer :: depsta => null()  ! station depth (C:meters/S:converted to km)
		real, dimension(:), pointer :: dis => null()     ! distance from source in km / rad
		real, dimension(:), pointer :: phi => null()     ! azimuth of receiver counted from south over east in rad
!---------------------------------------------------------------------------------------------------------------
	end type sar
!
 contains
!----------------------------------------------------------------------------------------------
	subroutine sar_readinfo(this,iunit,infofile,epimax)
	type (sar) :: this
	integer :: iunit,i,j,staid,ascii_dez,gg2gcflag
	character :: infofile*(*),line*80
	real :: epimax
!
	open(iunit, file=infofile, status='old')
	read(iunit,'(a13)') this%eventid
	read(iunit,'(a1,1x,i1,1x,i1)') this%csys,this%netflag,gg2gcflag
	read(iunit,*) this%year,this%month,this%day,this%hour,this%minute,this%sec,this%cendt,this%hdur
	read(iunit,*) this%xsrc,this%ysrc,this%depsrc,this%istyp
	if(this%istyp.eq.1) then
		read(iunit,*) (this%mt(i),i=1,6)
	else if(this%istyp.eq.0) then
		read(iunit,*) (this%force(i),i=1,3)
	else
		print *,'Unknown source type'
		stop
	endif
	if(this%csys.eq.'S') this%depsrc=this%depsrc*1.e-3
	if(gg2gcflag.eq.1.and.this%csys.eq.'S') call sar_geogr2geoce(this,'S',0)
	call sar_tm_epicar_globalcar(this)
!
	read(iunit,'(a)') line
!
!  if the first character of line is not a number (48 <= ascii-code <= 57)
!  line contains the filename with receiver coordinates
!  else it contains the number of receivers lines that follow
! 
	ascii_dez=ichar(line(1:1))
!	print *,line(1:1),ascii_dez
	if(ascii_dez.lt.48.or.ascii_dez.gt.57) then
		close(iunit)
		call sar_read_rcv_bin(this,iunit,line,epimax)
		return
	else
		read(line,*) this%nr
	endif
!
!  allocate space for receiver data
!
	allocate(this%xsta(this%nr),this%ysta(this%nr),this%name(this%nr),this%net(this%nr),this%depsta(this%nr))
	allocate(this%comp(this%nr),this%dis(this%nr),this%phi(this%nr),this%lfdnr(this%nr))
!
	j=1
	do i=1,this%nr
		if(this%netflag.eq.0) then
			read(iunit,'(i4,1x,a5,f10.4,1x,f10.4,1x,f5.0)') &
     &			        staid,this%name(j),this%xsta(j),this%ysta(j),this%depsta(j)
			this%net(j) = 'XX'; this%comp(j) = 'X'
		else if(this%netflag.eq.1) then
			read(iunit,'(i4,1x,a5,1x,a2,3x,f10.4,1x,f10.4,1x,f5.0)') &
     &			        staid,this%name(j),this%net(j),this%xsta(j),this%ysta(j),this%depsta(j)
			this%comp(j) = 'X'
		else if(this%netflag.eq.2) then
			read(iunit,'(i4,1x,a5,1x,a1,1x,a2,3x,f10.4,1x,f10.4,1x,f5.0)') & 
     &			        staid,this%name(j),this%comp(j),this%net(j),this%xsta(j),this%ysta(j),this%depsta(j)
		else if(this%netflag.eq.3) then
			read(iunit,'(i4,1x,a5,1x,a1,f10.4,1x,f10.4,1x,f5.0)') &
     &			        staid,this%name(j),this%comp(j),this%xsta(j),this%ysta(j),this%depsta(j)
			this%net(j) = 'XX'
		endif
		if(gg2gcflag.ge.1.and.this%csys.eq.'S') call sar_geogr2geoce(this,'R',j)
		call sar_epico(this,this%xsta(j),this%ysta(j),this%dis(j),this%phi(j))
		this%lfdnr(j)=i
		if(this%csys.eq.'S') then
			this%depsta(j)=this%depsta(j)*1.e-3
			if(this%dis(j)*180./pi.le.epimax.and.staid.lt.1000) j=j+1
		else
			if(this%dis(j).le.epimax.and.staid.lt.1000) j=j+1
		endif
	enddo
	close(iunit)
	this%nr=j-1
	this%nrorg=this%nr
	end subroutine sar_readinfo
!---------------------------------------------------------------------------------------
!> \brief Create basic sar object
!
	subroutine sar_CreateBasic(this,ntraces)
	type (sar) :: this
	integer :: ntraces
!
	allocate(this%name(ntraces),this%comp(ntraces),this%net(ntraces))
	allocate(this%lfdnr(ntraces),this%xsta(ntraces),this%ysta(ntraces),this%depsta(ntraces))
	allocate(this%dis(ntraces),this%phi(ntraces))
	this%nr = ntraces
	this%nrorg = ntraces
	end subroutine sar_CreateBasic
!----------------------------------------------------------------------------------------
!  deallocate sourceReceiver object
!
	subroutine deallocSourceReceiver(this)
	type (sar) :: this
	if(associated(this%xsta)) deallocate(this%xsta)
	if(associated(this%ysta)) deallocate(this%ysta)
	if(associated(this%depsta)) deallocate(this%depsta)
	if(associated(this%name)) deallocate(this%name)
	if(associated(this%net)) deallocate(this%net)
	if(associated(this%comp)) deallocate(this%comp)
	if(associated(this%dis)) deallocate(this%dis)
	if(associated(this%phi)) deallocate(this%phi)
	if(associated(this%lfdnr)) deallocate(this%lfdnr)
	end subroutine deallocSourceReceiver
!----------------------------------------------------------------------------------------
!  read receiver coordinates from a binary file
!  used with bornspecforce and bornspecgreen, 
!  because there may be a lot of scattering points
!
	subroutine sar_read_rcv_bin(this,iunit,filename,epimax)
	type (sar) :: this
	integer iunit,i,j
	real epimax
	character*(*) filename
!
	open(iunit,file=filename,form='unformatted', status = 'old')
	read(iunit) this%nr
	print *,'Number of surface receiver points: ',this%nr
!
!  allocate space for receiver data
!
	allocate(this%xsta(this%nr),this%ysta(this%nr),this%name(this%nr),this%net(this%nr),this%depsta(this%nr),this%depsta(this%nr))
	allocate(this%comp(this%nr),this%dis(this%nr),this%phi(this%nr),this%lfdnr(this%nr))
!
	read(iunit) (this%xsta(i),this%ysta(i),i=1,this%nr)
	j=1
	do i=1,this%nr
		print *,this%xsta(j),this%ysta(j)
		call sar_epico(this,this%xsta(j),this%ysta(j),this%dis(j),this%phi(j))
		this%lfdnr(j)=i
		this%name(j)='BORNS'
		this%comp(j)='X'
		this%net(j) = 'XX'
		if(this%csys.eq.'S') then
			if(this%dis(j)*180./pi.le.epimax) j=j+1
		else
			if(this%dis(j).le.epimax) j=j+1
		endif
	enddo
	close(iunit)
	this%nr=j-1
	this%nrorg=this%nr
	end subroutine sar_read_rcv_bin
!-----------------------------------------------------------------------------------------
!  compute epicentral coordinates of a point
!
	subroutine sar_epico(this,x,y,dis,phi)
	type (sar) :: this
	real x,y,dis,phi,pythag,thgeo,phigeo,thpol,phipol
!
	if(this%csys.eq.'S') then
		thpol=this%xsrc*pi/180.
		phipol=this%ysrc*pi/180.
		thgeo=x*pi/180.
		phigeo=y*pi/180.
		call geo2epi(thgeo,phigeo,thpol,phipol,dis,phi)
	else
		dis=pythag(x-this%xsrc,y-this%ysrc)*1.e-3
		phi=atan2(y-this%ysrc,x-this%xsrc)
	endif
	return
	end subroutine sar_epico
!-----------------------------------------------------------------------------------------
!  compute direction of wave propagation at (x,y) 
!  counted from south over east in spherical coordinates
!  or counted from x over y in cartesian coordinates
!
	subroutine sar_propdir(this,x,y,azi)
	type (sar) :: this
	real bazi,azi,thpol,phipol,thgeo,phigeo,x,y,dis
!
	if(this%csys.eq.'S') then
		thpol=x*pi/180.
		phipol=y*pi/180.		
		thgeo=this%xsrc*pi/180.
		phigeo=this%ysrc*pi/180.
		call geo2epi(thgeo,phigeo,thpol,phipol,dis,bazi)
		azi=bazi-pi
	else
		azi=atan2(y-this%ysrc,x-this%xsrc)
	endif
	end subroutine sar_propdir
!-------------------------------------------------------------------------------------
!  find index of a selected station
!
	subroutine sar_getstaindex(this,selsta,k)
	type (sar) :: this
	character selsta*(*)
	integer k,j
!
	k=0
	do j=1,this%nr
		if(trim(this%name(j)).eq.trim(selsta)) k=j
	enddo
	return
	end subroutine sar_getstaindex
!------------------------------------------------------------------------------------------
!  swap two stations j and k
!
	subroutine sar_swapstations(this,j,k)
	type (sar) :: this
	integer j,k
	real x,y,dis,phi,dep
	character name*6,net*2,comp*1
	integer lfdnr
!
	net=this%net(k)
	name=this%name(k)
	comp=this%comp(k)
	x=this%xsta(k)
	y=this%ysta(k)
	dis=this%dis(k)
	phi=this%phi(k)
	dep=this%depsta(k)
	lfdnr=this%lfdnr(k)
!
	this%net(k)=this%net(j)
	this%name(k)=this%name(j)
	this%comp(k)=this%comp(j)
	this%xsta(k)=this%xsta(j)
	this%ysta(k)=this%ysta(j)
	this%dis(k)=this%dis(j)
	this%phi(k)=this%phi(j)
	this%depsta(k)=this%depsta(j)
	this%lfdnr(k)=this%lfdnr(j)
!	
	this%net(j)=net
	this%name(j)=name
	this%comp(j)=comp
	this%xsta(j)=x
	this%ysta(j)=y
	this%dis(j)=dis
	this%phi(j)=phi
	this%depsta(j)=dep
	this%lfdnr(j)=lfdnr
!
	end subroutine sar_swapstations
!--------------------------------------------------------------------------------------------
!  set number of receivers
!
	subroutine sar_setnr(this,n)
	type (sar) :: this
	integer n
!
	this%nr=n
	end subroutine sar_setnr
!-----------------------------------------------------------------------
!  set components to one common value
!
	subroutine sar_setComps(this,comp)
	type (sar) :: this
	character (len=*) :: comp
!
	this%comp = comp
	end subroutine sar_setComps
!--------------------------------------------------------------------------
!  convert sar_dis from rad to km
!
	subroutine sar_rad2km(this,rearth)
	type (sar) :: this
	integer j
	real rearth
!
	if(this%csys .eq. 'S') then
		do j=1,this%nr
			this%dis(j)=this%dis(j)*rearth
		enddo
	endif
	end subroutine sar_rad2km
!-----------------------------------------------------------------------
!  calculate hour, min, sec of rupture time from event time
!
	subroutine sar_eventhms(this,hour,minute,sec)
	type (sar) :: this
	integer hour,minute
	real sec,tanf
!
	tanf=this%sec+this%cendt+60.*(float(this%minute)+60.*float(this%hour))
	hour=int(tanf/3600.)
	minute=int((tanf-hour*3600.)/60.)
	sec=tanf-60.*(float(minute)+float(hour)*60.)
	end subroutine sar_eventhms
!---------------------------------------------------------------
!  compute scalar moment from moment tensor in Nm
!
	subroutine sar_scalar_moment(this,rmsc)
	type (sar) :: this
	real rmsc
	double precision rm2
!
	if(this%istyp.ne.1) then
		print *,'not reasonable for force source'
		stop
	endif
	rm2=dble(this%mt(1))**2+dble(this%mt(2))**2+dble(this%mt(3))**2+ &
     &	    2.*dble(this%mt(4))**2+2.*dble(this%mt(5))**2+2.*dble(this%mt(6))**2
	rmsc=sngl(dsqrt(rm2/2.))
	end subroutine sar_scalar_moment
!----------------------------------------------------------------
!  compute transformation matrix from epicentral cartesian coordinates
!  to global cartesian coordinates as described in born-scattering.tex,
!  Matrix T
!
	subroutine sar_tm_epicar_globalcar(this)
	type (sar) :: this
	real deg2rad,ct,st,cf,sf
	integer i,j
!
	deg2rad=pi/180.
	if(this%csys.eq.'S') then
		ct=cos(this%xsrc*deg2rad)
		st=sin(this%xsrc*deg2rad)
		cf=cos(this%ysrc*deg2rad)
		sf=sin(this%ysrc*deg2rad)
		this%tm_ec_gc(1,1)=ct*cf
		this%tm_ec_gc(1,2)=-sf
		this%tm_ec_gc(1,3)=st*cf
		this%tm_ec_gc(2,1)=ct*sf
		this%tm_ec_gc(2,2)=cf
		this%tm_ec_gc(2,3)=st*sf
		this%tm_ec_gc(3,1)=-st
		this%tm_ec_gc(3,2)=0.
		this%tm_ec_gc(3,3)=ct
	else
		do i=1,3
			do j=1,3
				this%tm_ec_gc(i,j)=0.
			enddo
			this%tm_ec_gc(i,i)=1.
		enddo
	endif
	end subroutine sar_tm_epicar_globalcar
!-----------------------------------------------------------------------
!  compute transformation matrix from epicentral spherical 
!  to global cartesian for one observation point
!  as described in born-scattering.tex (Matrix M)
!  delta: angular distance between Q,E and S
!  azi: azimuth (from south over east) of S as viewed from Q,E
!
	subroutine sar_tm_epispher_globalcar(this,delta,azi,tm)
	type (sar) :: this
	real delta,azi,tm(3,3),ca,sa,cd,sd
	integer i
!
	ca=cos(azi)
	sa=sin(azi)
	if(this%csys.eq.'S') then
		cd=cos(delta)
		sd=sin(delta)
		do i=1,3
			tm(i,1)=this%tm_ec_gc(i,1)*sd*ca + this%tm_ec_gc(i,2)*sd*sa + this%tm_ec_gc(i,3)*cd
			tm(i,2)=this%tm_ec_gc(i,1)*cd*ca + this%tm_ec_gc(i,2)*cd*sa - this%tm_ec_gc(i,3)*sd
			tm(i,3)=-this%tm_ec_gc(i,1)*sa   + this%tm_ec_gc(i,2)*ca
		enddo
	else
		tm(1,1)=0.
		tm(2,1)=0.
		tm(3,1)=1.
		tm(1,2)=ca
		tm(2,2)=sa
		tm(3,2)=0.
		tm(1,3)=-sa
		tm(2,3)=ca
		tm(3,3)=0.
	endif
	end subroutine sar_tm_epispher_globalcar
!--------------------------------------------------------
!   convert from geographic to geocentric colatitude
!   g: geographic colatitude (input)
!   c: geocentric colatitude (output)
!
!   tan c = (1+2*eps) * tan g
!   
!   Three remarks:
!   (1) If g > 90 deg, the tangents of g and c are negative.
!       To get a value between 0 and 180, add 180 deg to c.
!   (2) Near g = 90 deg, the tangent gets extremely large.
!       Apply a Taylor expansion about 90 deg and keep
!       terms linear in 90-g and 90-c and also terms (90-g)*eps.
!       Then  c = 90 -(90-g)/(1+2*eps) = 90-(90-g)*(1-2*eps)
!               = 90 -(90-g) + 2*eps*(90-g)
!               = g + 2*eps*(90-g)
!   (3) If g=180 deg then tan(g)=0 and c=0 instead of 180!
!---------------------------------------------------------------
	subroutine sar_geogr2geoce(this,cflag,j)
	type (sar) :: this
	character cflag*1
	integer j
	real g,c,pihmg,pih,eps
	data eps/3.3355570e-3/
!
	if(cflag.eq.'S') g = this%xsrc*pi/180.
	if(cflag.eq.'R') g = this%xsta(j)*pi/180.
!
	pih=.5*pi
	pihmg=pih-g
	if( abs(pihmg).lt.1.e-6 ) then
		c=g+2.*eps*pihmg
	else
		c=atan( (1.+2.*eps)*tan(g) )
		if(c.lt.0.) c=c+pi
	endif
	if(abs(c-g).gt..9*pi) c=c+pi
!
	if(cflag.eq.'S') this%xsrc = c*180./pi
	if(cflag.eq.'R') this%xsta(j) = c*180./pi
!
	end subroutine sar_geogr2geoce
!-----------------------------------------------------------------
!  sort traces according to phi
!
	subroutine sar_sortPhi(this,indx,phisorted)
	use sortArray
	type (sar) :: this
	integer, dimension(:) :: indx
	real, dimension(:) :: phisorted
!
	if(size(indx) < size(this%phi)) then
		print *,'<sar_sortPhi>: indx-array to small'
		stop
	endif
	phisorted = this%phi
	call sort(phisorted,indx)
	end subroutine sar_sortPhi
!-----------------------------------------------------------------
!  return number of traces
!
	integer function sar_getNumberTraces(this)
	type (sar), intent(in) :: this
!
	sar_getNumberTraces=this%nr
	end function sar_getNumberTraces
!-----------------------------------------------------------------
!  return eventid
!
	character (len=13) function sar_getEventid(this)
	type (sar), intent(in) :: this
!
	sar_getEventid=this%eventid
	end function sar_getEventid
!--------------------------------------------------------------------------
!  construct date string from event date (year with two digits)
!
	character (len=6) function sar_getDateSource(this)
	type (sar), intent(in) :: this
!
	write(sar_getDateSource,'(i2.2,i2.2,i2.2)') this%year,this%month,this%day
	end function sar_getDateSource
!--------------------------------------------------------------------------
!  construct date string from event date (year with 4 digits digits)
!  Beg your pardon! Fails after 2069 !
!
	character (len=8) function sar_getFullYearDateSource(this)
	type (sar), intent(in) :: this
	integer :: year
!
	if (this%year < 100) then
		if (this%year > 70) then
			year = this%year+1900
		else
			year = this%year+2000
		endif
	else
		year = this%year
	endif
	write(sar_getFullYearDateSource,'(i4.4,i2.2,i2.2)') year,this%month,this%day
	end function sar_getFullYearDateSource
!----------------------------------------------------------------------
!  construct time string from event time
!
	character (len=10) function sar_getTimeSource(this)
	type (sar), intent(in) :: this
!
	write(sar_getTimeSource,'(i2.2,i2.2,f6.3)') this%hour,this%minute,this%sec+this%cendt
	end function sar_getTimeSource
!---------------------------------------------------------------
!   calculate source time in seconds after midnight
!
	real function sar_getTSAMSource(this)
	type (sar), intent(in) :: this
!	
	sar_getTSAMSource = this%sec+this%cendt+60.*(float(this%minute)+60.*float(this%hour))
	end function sar_getTSAMSource
!-----------------------------------------------------------------
!  return coordinate system
!
	character (len=1) function sar_CoordinateSystem(this)
	type (sar), intent(in) :: this
	sar_CoordinateSystem = this%csys
	end function sar_CoordinateSystem
!--------------------------------------------------------------------
!  get distance of a selected trace by index
!
	real function sar_getDistanceSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getDistanceSelected = this%dis(ksta)
	end function sar_getDistanceSelected
!--------------------------------------------------------------------
!  get phi of a selected trace by index
!
	real function sar_getPhiSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getPhiSelected = this%phi(ksta)
	end function sar_getPhiSelected
!--------------------------------------------------------------------
!  get propagation direction for a selected trace by index
!
	real function sar_getPropdirSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	real :: azi,thpol,phipol,thgeo,phigeo,dis,bazi
!
	if(this%csys.eq.'S') then
		thpol =this%xsta(ksta)*pi/180.
		phipol=this%ysta(ksta)*pi/180.		
		thgeo =this%xsrc*pi/180.
		phigeo=this%ysrc*pi/180.
		call geo2epi(thgeo,phigeo,thpol,phipol,dis,bazi)
		azi=bazi-pi
	else
		azi=atan2(this%ysta(ksta)-this%ysrc,this%xsta(ksta)-this%xsrc)
	endif
	sar_getPropdirSelected = azi
	end function sar_getPropdirSelected
!--------------------------------------------------------------------
!  get net id of a selected trace by index
!
	character (len=2) function sar_getNetIdSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getNetIdSelected = this%net(ksta)
	end function sar_getNetIdSelected
!--------------------------------------------------------------------
!  get comp of a selected trace by index
!
	character (len=1) function sar_getCompSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getCompSelected = this%comp(ksta)
	end function sar_getCompSelected
!--------------------------------------------------------------------
!  get station of a selected trace by index
!
	character (len=5) function sar_getStationSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getStationSelected = this%name(ksta)
	end function sar_getStationSelected
!--------------------------------------------------------------------
!  get x-coordinate of a selected trace by index
!
	real function sar_getXSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getXSelected = this%xsta(ksta)
	end function sar_getXSelected
!--------------------------------------------------------------------
!  get y-coordinate of a selected trace by index
!
	real function sar_getYSelected(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getYSelected = this%ysta(ksta)
	end function sar_getYSelected
!--------------------------------------------------------------------
!  get pointer to x-coordinates of traces
!
	function sar_getX(this) result(x)
	type (sar), intent(in) :: this
	real, dimension(:), pointer :: x 
	x => this%xsta
	end function sar_getX
!--------------------------------------------------------------------
!  get pointer to y-coordinates of traces
!
	function sar_getY(this) result(y)
	type (sar), intent(in) :: this
	real, dimension(:), pointer :: y 
	y => this%ysta
	end function sar_getY
!--------------------------------------------------------------------
!  get station depth
!
	real function sar_getStationDepth(this,ksta)
	type (sar), intent(in) :: this
	integer, intent(in) :: ksta
	sar_getStationDepth = this%depsta(ksta)
	end function sar_getStationDepth
!--------------------------------------------------------------------
!  get x-coordinate of source
!
	real function sar_getXSource(this)
	type (sar), intent(in) :: this
	sar_getXSource = this%xsrc
	end function sar_getXSource
!--------------------------------------------------------------------
!  get y-coordinate of source
!
	real function sar_getYSource(this)
	type (sar), intent(in) :: this
	sar_getYSource = this%ysrc
	end function sar_getYSource
!--------------------------------------------------------------------
!  get year of source
!
	integer function sar_getYearSource(this)
	type (sar), intent(in) :: this
	sar_getYearSource = this%year
	end function sar_getYearSource
!--------------------------------------------------------------------
!  get month of source
!
	integer function sar_getMonthSource(this)
	type (sar), intent(in) :: this
	sar_getMonthSource = this%month
	end function sar_getMonthSource
!--------------------------------------------------------------------
!  get day of source
!
	integer function sar_getDaySource(this)
	type (sar), intent(in) :: this
	sar_getDaySource = this%day
	end function sar_getDaySource
!--------------------------------------------------------------------
!  get day of year of source
!
	integer function sar_getDayOfYearSource(this)
	type (sar), intent(in) :: this
	sar_getDayOfYearSource = getDayOfYear(this%year,this%month,this%day)
	end function sar_getDayOfYearSource
!-----------------------------------------------------------------
!  get netflag
!
	integer function sar_getNetFlag(this)
	type (sar), intent(in) :: this
	sar_getNetFlag = this%netflag
	end function sar_getNetFlag
!-------------------------------------------------------------------
!  get source depth
!
	real function sar_getSourceDepth(this)
	type (sar), intent(in) :: this
	sar_getSourceDepth = this%depsrc
	end function sar_getSourceDepth
!-------------------------------------------------------------------
!  get moment tensor
!
	subroutine sar_getMomentTensor(this,rm)
	type (sar), intent(in) :: this
	real, dimension(:) :: rm
	rm = this%mt
	end subroutine sar_getMomentTensor
!-------------------------------------------------------------------
!  get force vector
!
	subroutine sar_getForce(this,f)
	type (sar), intent(in) :: this
	real, dimension(:) :: f
	f = this%force
	end subroutine sar_getForce
!-----------------------------------------------------------------
!  get source type
!
	integer function sar_getSourceType(this)
	type (sar), intent(in) :: this
	sar_getSourceType = this%istyp
	end function sar_getSourceType
!-------------------------------------------------------------------
!  get pointer to distances
!
	function sar_getDistances(this) result(dis)
	type (sar), intent(in) :: this
	real, dimension(:), pointer :: dis
	dis => this%dis
	end function sar_getDistances
!-------------------------------------------------------------------
!  get pointer to azimuths
!
	function sar_getAzimuths(this) result(phi)
	type (sar), intent(in) :: this
	real, dimension(:), pointer :: phi
	phi => this%phi
	end function sar_getAzimuths
!------------------------------------------------------------------------
!  calculate time in full seconds after midnight from hour minute seconds
!
	function sar_getSourceTimeFullSeconds(this) result(tfs)
	type (sar), intent(in) :: this
	integer :: tfs
	tfs = (60*this%hour+this%minute)*60+int(this%sec)
	end function sar_getSourceTimeFullSeconds
!------------------------------------------------------------------------
!  calculate nanosecond part of source time after midnight from seconds
!
	function sar_getSourceTimeNanoSeconds(this) result(tns)
	type (sar), intent(in) :: this
	integer :: tns
	tns = nint((this%sec-int(this%sec))*1.e9)
	end function sar_getSourceTimeNanoSeconds
!----------------------------------------------------------------------
!> \brief set j-th trace name, component and network
!
	subroutine sar_setNameCompNet(this,j,name,comp,net)
	type (sar) :: this
	character (len=*) :: name,comp,net
	integer :: j
	this%lfdnr(j) = j
	this%name(j) = name
	this%comp(j) = comp
	this%net(j) = net
	this%netflag = 2
	end subroutine sar_setNameCompNet
!----------------------------------------------------------------------
!> \brief set j-th trace xsta,ysta,depsta,dis,phi
!
	subroutine sar_setStationCoordinates(this,j,x,y,z)
	type (sar) :: this
	real :: x,y,z
	integer :: j
	this%lfdnr(j) = j
	if (this%csys == 'S') then
		this%xsta(j) = 90.-x
		this%ysta(j) = y
		this%depsta(j) = z*1.e-3
	else
		this%xsta(j) = x
		this%ysta(j) = y
		this%depsta(j) = z
	endif
	call sar_epico(this,this%xsta(j),this%ysta(j),this%dis(j),this%phi(j))
	end subroutine sar_setStationCoordinates
!----------------------------------------------------------------------
!> \brief set source coordinates
!
	subroutine sar_setSourceCoordinates(this,x,y,z)
	type (sar) :: this
	real :: x,y,z
	if (this%csys == 'S') then
		this%xsrc = 90.-x
		this%ysrc = y
		this%depsrc = z*1.e-3
	else
		this%xsrc = x
		this%ysrc = y
		this%depsrc = z
	endif
	call sar_tm_epicar_globalcar(this)
	end subroutine sar_setSourceCoordinates
!----------------------------------------------------------------------
!> \brief set origin date and time
!
	subroutine sar_setOriginDateTime(this,timestring)
	type (sar) :: this
	character (len=*) :: timestring
	integer :: fs,ns
	read(timestring,'(i4.4,2i2.2,1x,3i2.2,1x,i9.9)') this%year,this%month,this%day,this%hour,this%minute,fs,ns
	this%sec = 1.0*fs+1.e-9*ns
	this%cendt = 0.0
	end subroutine sar_setOriginDateTime
!----------------------------------------------------------------------
!> \brief set source type
!
	subroutine sar_setSourceType(this,istyp)
	type (sar) :: this
	integer :: istyp
	this%istyp = istyp
	end subroutine sar_setSourceType
!----------------------------------------------------------------------
!> \brief set coordinate system
!
	subroutine sar_setCoordinateSystem(this,csys)
	type (sar) :: this
	character (len=*) :: csys
	this%csys = csys
	end subroutine sar_setCoordinateSystem
!----------------------------------------------------------------------
!> \brief set event id
!
	subroutine sar_setEventId(this,eventid)
	type (sar) :: this
	character (len=*) :: eventid
	integer :: n
	n = len_trim(eventid)
	this%eventid(1:n) = eventid
	this%eventid(n+1:len(this%eventid)) = ' '
	end subroutine sar_setEventId
!----------------------------------------------------------------------
!> \brief Extend info structure by new traces
!
	subroutine sar_extend(this,n)
	type (sar) :: this
	integer :: n
	this%nr = n
	this%nrorg = n
	this%lfdnr => reallocate(this%lfdnr,n)
	this%name => reallocate(this%name,n)
	this%comp => reallocate(this%comp,n)
	this%net => reallocate(this%net,n)
	this%xsta => reallocate(this%xsta,n)
	this%ysta => reallocate(this%ysta,n)
	this%depsta => reallocate(this%depsta,n)
	this%dis => reallocate(this%dis,n)
	this%phi => reallocate(this%phi,n)
	end subroutine sar_extend
!
 end module sourceReceiver
