c---------------------------------------------------------------------
c  Procedures related to source and receiver (sar) information
c  contained in info-Files
c
c  cartesian coordinates: 
c	rxy with phi from x to y
c	forces: 1=fr, 2=fx, 3=fy
c	moments: 1=Mrr, 2=Mxx, 3=Myy, 4=Mrx, 5=Mry, 6=Mxy
c	displacement components: either RLT with T pointing left if one looks from source to receiver
c	                         or RXY
c	sar_depsrc = source depth in meters
c	sar_xsrc = x-coordinate of source in meters
c	sar_ysrc = y-coordinate of source in meters
c	sar_xsta = x-coordinate of station in meters
c	sar_ysta = y-coordinate of station in meters
c	sar_dis(i) in kilometers
c	sar_phi(i) in radians
c
c  spherical coordinates: 
c	rtf with phi and from south over east
c	forces: 1=fr, 2=ft, 3=ff using basis vectors at the source point!
c	moments: 1=Mrr, 2=Mtt, 3=Mff, 4=Mrt, 5=Mrf, 6=Mff using basis vectors at the source point!
c	displacement components: either RLT with T pointing left if one looks from source to receiver
c	                         or RNE
c
c	sar_depsrc = source depth in km
c	sar_xsrc = Colatitide of source in degrees
c	sar_ysrc = Longitude of source in degrees
c	sar_xsta = Colatitude of station in degrees
c	syr_ysta = Longitude of station in degrees
c	sar_dis(i) in radians
c	sar_phi(i) in radians
c
c	netflag = 0: ID STA COLAT LON
c	netflag = 1: ID STA NET COLAT LON
c	netflag = 2: ID STA COMP NET COLAT LON
c	netflag = 3: ID STA COMP COLAT LON
c
c	gg2gcflag: convert from geographic to geocentric coordinates
c	           = 0:  do not transform at all
c	           = 1:  transform both source and receiver coordinates
c	           = 2:  transform receiver coordinates only
c-----------------------------------------------------------------------
	subroutine sar_readinfo(iunit,infofile,epimax)
	include 'sardim.h'
	include 'sar.h'
	include 'pi.h'
	integer iunit,i,j,staid,ascii_dez,gg2gcflag
	character infofile*(*),line*80
	real epimax
c
	open(iunit,file=infofile)
	read(iunit,'(a13)') sar_eventid
	read(iunit,'(a1,1x,i1,1x,i1)') sar_csys,sar_netflag,gg2gcflag
	read(iunit,*) sar_year,sar_month,sar_day,sar_hour,sar_minute,sar_sec,sar_cendt,sar_hdur
	read(iunit,*) sar_xsrc,sar_ysrc,sar_depsrc,sar_istyp
	if(sar_istyp.eq.1) then
		read(iunit,*) (sar_mt(i),i=1,6)
	else if(sar_istyp.eq.0) then
		read(iunit,*) (sar_force(i),i=1,3)
	else
		print *,'Unknown source type'
		stop
	endif
	if(sar_csys.eq.'S') sar_depsrc=sar_depsrc*1.e-3
	if(gg2gcflag.eq.1.and.sar_csys.eq.'S') call sar_geogr2geoce('S',0)
	call sar_tm_epicar_globalcar(sar_tm_ec_gc)
c
	read(iunit,'(a)') line
c
c  if the first character of line is not a number (48 <= ascii-code <= 57)
c  line contains the filename with receiver coordinates
c  else it contains the number of receivers lines that follow
c 
	ascii_dez=iachar(line(1:1))
c	print *,line(1:1),ascii_dez
	if(ascii_dez.lt.48.or.ascii_dez.gt.57) then
		call sar_zlimits(iunit)
		close(iunit)
		call sar_read_rcv_bin(iunit,line,epimax)
		return
	else
		read(line,*) sar_nr
	endif
	if(sar_nr.gt.nstatt) then
		print *,'More receivers than dimensioned'
		print *,'Increase NSTATT to at least ',sar_nr
		stop
	endif 
	j=1
	do i=1,sar_nr
		if(sar_netflag.eq.0) then
			read(iunit,'(i4,1x,a5,f10.4,1x,f10.4)') 
     1			        staid,sar_name(j),sar_xsta(j),sar_ysta(j)
		else if(sar_netflag.eq.1) then
			read(iunit,'(i4,1x,a5,1x,a2,3x,f10.4,1x,f10.4)') 
     1			        staid,sar_name(j),sar_net(j),sar_xsta(j),sar_ysta(j)
		else if(sar_netflag.eq.2) then
			read(iunit,'(i4,1x,a5,1x,a1,1x,a2,3x,f10.4,1x,f10.4)') 
     1			        staid,sar_name(j),sar_comp(j),sar_net(j),sar_xsta(j),sar_ysta(j)
		else if(sar_netflag.eq.3) then
			read(iunit,'(i4,1x,a5,1x,a1,f10.4,1x,f10.4)') 
     1			        staid,sar_name(j),sar_comp(j),sar_xsta(j),sar_ysta(j)
		endif
		if(gg2gcflag.ge.1.and.sar_csys.eq.'S') call sar_geogr2geoce('R',j)
		call sar_epico(sar_xsta(j),sar_ysta(j),sar_dis(j),sar_phi(j))
		sar_lfdnr(j)=i
		if(sar_csys.eq.'S') then
			if(sar_dis(j)*180./pi.le.epimax.and.staid.lt.1000) j=j+1
		else
			if(sar_dis(j).le.epimax.and.staid.lt.1000) j=j+1
		endif
	enddo
	call sar_zlimits(iunit)
	close(iunit)
	sar_nr=j-1
	sar_nrorg=sar_nr
c
	return
	end
c----------------------------------------------------------------------
c  read receiver coordinates from a binary file
c  used with bornspecforce and bornspecgreen, 
c  because there may be a lot of scattering points
c
	subroutine sar_read_rcv_bin(iunit,filename,epimax)
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
	integer iunit,i,j
	real epimax
	character*(*) filename
c
	open(iunit,file=filename,form='unformatted')
	read(iunit) sar_nr
	print *,'Number of surface receiver points: ',sar_nr
	if(sar_nr.gt.nstatt) then
		print *,'more receivers than dimensioned'
		print *,sar_nr,nstatt
		stop
	endif
	read(iunit) (sar_xsta(i),sar_ysta(i),i=1,sar_nr)
	j=1
	do i=1,sar_nr
		print *,sar_xsta(j),sar_ysta(j)
		call sar_epico(sar_xsta(j),sar_ysta(j),sar_dis(j),sar_phi(j))
		sar_lfdnr(j)=i
		sar_name(j)='BORNS'
		sar_comp(j)='Q'
		if(sar_csys.eq.'S') then
			if(sar_dis(j)*180./pi.le.epimax) j=j+1
		else
			if(sar_dis(j).le.epimax) j=j+1
		endif
	enddo
	close(iunit)
	sar_nr=j-1
	sar_nrorg=sar_nr
	return
	end
c----------------------------------------------------------------------
c  compute epicentral coordinates of a point
c
	subroutine sar_epico(x,y,dis,phi)
	real x,y,dis,phi,pythag,thgeo,phigeo,thpol,phipol
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
c
	if(sar_csys.eq.'S') then
		thpol=sar_xsrc*pi/180.
		phipol=sar_ysrc*pi/180.
		thgeo=x*pi/180.
		phigeo=y*pi/180.
		call geo2epi(thgeo,phigeo,thpol,phipol,dis,phi)
	else
		dis=pythag(x-sar_xsrc,y-sar_ysrc)*1.e-3
		phi=atan2(y-sar_ysrc,x-sar_xsrc)
	endif
	return
	end
c---------------------------------------------------------------------
c  compute direction of wave propagation at (x,y) 
c  counted from south over east in spherical coordinates
c  or counted from x over y in cartesian coordinates
c
	subroutine sar_propdir(x,y,azi)
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
	real bazi,azi,thpol,phipol,thgeo,phigeo,x,y,dis
c
	if(sar_csys.eq.'S') then
		thpol=x*pi/180.
		phipol=y*pi/180.		
		thgeo=sar_xsrc*pi/180.
		phigeo=sar_ysrc*pi/180.
		call geo2epi(thgeo,phigeo,thpol,phipol,dis,bazi)
		azi=bazi-pi
	else
		azi=atan2(y-sar_ysrc,x-sar_xsrc)
	endif
	return
	end
c------------------------------------------------------------------------
c  find index of a selected station
c
	subroutine sar_getstaindex(selsta,k)
	include 'sardim.h'
	include 'sar.h'
	character selsta*(*)
	integer k,j,ic1,ic2,jc1,jc2
c
	k=0
	call trimleft(selsta,ic1)
	call trimright(selsta,ic2)
	do j=1,sar_nr
		call trimleft(sar_name(j),jc1)
		call trimright(sar_name(j),jc2)
		if(sar_name(j)(jc1:jc2).eq.selsta(ic1:ic2)) k=j
	enddo
	return
	end
c------------------------------------------------------------------------
c  swap two stations j and k
c
	subroutine sar_swapstations(j,k)
	include 'sardim.h'
	include 'sar.h'
	integer j,k
	real x,y,dis,phi
	character*6 name
c
	name=sar_name(k)
	x=sar_xsta(k)
	y=sar_ysta(k)
	dis=sar_dis(k)
	phi=sar_phi(k)
c
	sar_name(k)=sar_name(j)
	sar_xsta(k)=sar_xsta(j)
	sar_ysta(k)=sar_ysta(j)
	sar_dis(k)=sar_dis(j)
	sar_phi(k)=sar_phi(j)
c	
	sar_name(j)=name
	sar_xsta(j)=x
	sar_ysta(j)=y
	sar_dis(j)=dis
	sar_phi(j)=phi
c
	return
	end
c-------------------------------------------------------------------------
c  set number of receivers
c
	subroutine sar_setnr(n)
	include 'sardim.h'
	include 'sar.h'
	integer n
c
	sar_nr=n
	return
	end
c--------------------------------------------------------------------------
c  convert sar_dis from rad to km
c
	subroutine sar_rad2km(rearth)
	include 'sardim.h'
	include 'sar.h'
	integer j
	real rearth
c
	do j=1,sar_nr
		sar_dis(j)=sar_dis(j)*rearth
	enddo
	return
	end
c--------------------------------------------------------------------------
c  construct date string from event date
c
	subroutine sar_eventsdate(sdate)
	include 'sardim.h'
	include 'sar.h'
	character sdate*(*)
c
	write(sdate,'(i2.2,i2.2,i2.2)') sar_year,sar_month,sar_day
	return
	end
c----------------------------------------------------------------------
c  construct time string from event time
c
	subroutine sar_eventstime(stime)
	include 'sardim.h'
	include 'sar.h'
	character stime*(*)
c
	write(stime,'(i2.2,i2.2,f6.3)') sar_hour,sar_minute,sar_sec+sar_cendt
	return
	end
c-----------------------------------------------------------------------
c  calculate hour, min, sec of rupture time from event time
c
	subroutine sar_eventhms(hour,minute,sec)
	include 'sardim.h'
	include 'sar.h'
	integer hour,minute
	real sec,tanf
c
	tanf=sar_sec+sar_cendt+60.*(float(sar_minute)+60.*float(sar_hour))
	hour=int(tanf/3600.)
	minute=int((tanf-hour*3600.)/60.)
	sec=tanf-60.*(float(minute)+float(hour)*60.)
	return
	end
c---------------------------------------------------------------
c   calculate source time in seconds after midnight
c
	subroutine sar_event_time(therd)
	include 'sardim.h'
	include 'sar.h'
	real therd
c	
	therd=sar_sec+sar_cendt+60.*(float(sar_minute)+60.*float(sar_hour))
	return
	end
c---------------------------------------------------------------
c  compute scalar moment from moment tensor in Nm
c
	subroutine sar_scalar_moment(rmsc)
	include 'sardim.h'
	include 'sar.h'
	real rmsc
	double precision rm2
c
	if(sar_istyp.ne.1) then
		print *,'not reasonable for force source'
		stop
	endif
	rm2=dble(sar_mt(1))**2+dble(sar_mt(2))**2+dble(sar_mt(3))**2+
     1	    2.*dble(sar_mt(4))**2+2.*dble(sar_mt(5))**2+2.*dble(sar_mt(6))**2
	rmsc=dsqrt(rm2/2.)
	return
	end
c----------------------------------------------------------------
c  compute transformation matrix from epicentral cartesian coordinates
c  to global cartesian coordinates as described in born-scattering.tex,
c  Matrix T
c
	subroutine sar_tm_epicar_globalcar(tm)
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
	real tm(3,3),deg2rad,ct,st,cf,sf
	integer i,j
c
	deg2rad=pi/180.
	if(sar_csys.eq.'S') then
		ct=cos(sar_xsrc*deg2rad)
		st=sin(sar_xsrc*deg2rad)
		cf=cos(sar_ysrc*deg2rad)
		sf=sin(sar_ysrc*deg2rad)
		tm(1,1)=ct*cf
		tm(1,2)=-sf
		tm(1,3)=st*cf
		tm(2,1)=ct*sf
		tm(2,2)=cf
		tm(2,3)=st*sf
		tm(3,1)=-st
		tm(3,2)=0.
		tm(3,3)=ct
	else
		do i=1,3
			do j=1,3
				tm(i,j)=0.
			enddo
			tm(i,i)=1.
		enddo
	endif
	return
	end
c-------------------------------------------------------
c  compute transformation matrix from epicentral spherical 
c  to global cartesian for one observation point
c  as described in born-scattering.tex (Matrix M)
c  delta: angular distance between Q,E and S
c  azi: azimuth (from south over east) of S as viewed from Q,E
c
	subroutine sar_tm_epispher_globalcar(delta,azi,tm)
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
	real delta,azi,tm(3,3),ca,sa,cd,sd
	integer i
c
	ca=cos(azi)
	sa=sin(azi)
	if(sar_csys.eq.'S') then
		cd=cos(delta)
		sd=sin(delta)
		do i=1,3
			tm(i,1)=sar_tm_ec_gc(i,1)*sd*ca + sar_tm_ec_gc(i,2)*sd*sa + sar_tm_ec_gc(i,3)*cd
			tm(i,2)=sar_tm_ec_gc(i,1)*cd*ca + sar_tm_ec_gc(i,2)*cd*sa - sar_tm_ec_gc(i,3)*sd
			tm(i,3)=-sar_tm_ec_gc(i,1)*sa   + sar_tm_ec_gc(i,2)*ca
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
	return
	end
c--------------------------------------------------------
c  read eventual depth limits for radial nodes

	subroutine sar_zlimits(iunit)
	include 'sardim.h'
	include 'sar.h'
	integer iunit
	character*7 line

	read(iunit,'(a7)',end=99) line
	if(line(1:7).eq.'ZLIMITS') then
		read(iunit,*) sar_zmin,sar_zmax
		sar_zlim_flag=1
		return
	endif
 99	sar_zlim_flag=0
	sar_zmin=-1.
	sar_zmax=-1.
	return
	end
c--------------------------------------------------------
c   convert from geographic to geocentric colatitude
c   g: geographic colatitude (input)
c   c: geocentric colatitude (output)
c
c   tan c = (1+2*eps) * tan g
c   
c   Three remarks:
c   (1) If g > 90 deg, the tangents of g and c are negative.
c       To get a value between 0 and 180, add 180 deg to c.
c   (2) Near g = 90 deg, the tangent gets extremely large.
c       Apply a Taylor expansion about 90 deg and keep
c       terms linear in 90-g and 90-c and also terms (90-g)*eps.
c       Then  c = 90 -(90-g)/(1+2*eps) = 90-(90-g)*(1-2*eps)
c               = 90 -(90-g) + 2*eps*(90-g)
c               = g + 2*eps*(90-g)
c   (3) If g=180 deg then tan(g)=0 and c=0 instead of 180!
c---------------------------------------------------------------
	subroutine sar_geogr2geoce(cflag,j)
	include 'sardim.h'
	include 'sar.h'
	include 'pis.h'
	character cflag*1
	integer j
	real g,c,pihmg,pih,eps
	data eps/3.3355570e-3/
c
	if(cflag.eq.'S') g = sar_xsrc*pi/180.
	if(cflag.eq.'R') g = sar_xsta(j)*pi/180.
c
	pih=.5*pi
	pihmg=pih-g
	if( abs(pihmg).lt.1.e-6 ) then
		c=g+2.*eps*pihmg
	else
		c=atan( (1.+2.*eps)*tan(g) )
		if(c.lt.0.) c=c+pi
	endif
	if(abs(c-g).gt..9*pi) c=c+pi
c
	if(cflag.eq.'S') sar_xsrc = c*180./pi
	if(cflag.eq.'R') sar_xsta(j) = c*180./pi
c
	return
	end
c-----------------------------------------------------------------
c  return number of traces
c
	subroutine sar_getNumberTraces(nr)
	include 'sardim.h'
	include 'sar.h'
	integer nr
c
	nr=sar_nr
	return
	end
c-----------------------------------------------------------------
c  return eventid
c
	subroutine sar_getEventid(eventid)
	include 'sardim.h'
	include 'sar.h'
	character*(*) eventid
c
	eventid=sar_eventid
	return
	end
c-----------------------------------------------------------------
c  return station names
c-
	subroutine sar_getStations(nr,staname)
	include 'sardim.h'
	include 'sar.h'
	integer i,nr
	character*(*) staname(nr)
c
	nr=sar_nr
	do i=1,sar_nr
		staname(i)=sar_name(i)
	enddo
	return
	end
c-----------------------------------------------------------------
c  return components
c-
	subroutine sar_getComps(nr,comp)
	include 'sardim.h'
	include 'sar.h'
	integer i,nr
	character*(*) comp(nr)
c
	nr=sar_nr
	do i=1,sar_nr
		comp(i)=sar_comp(i)
	enddo
	return
	end
c-----------------------------------------------------------------
c  return net id
c-
	subroutine sar_getNetIds(nr,net)
	include 'sardim.h'
	include 'sar.h'
	integer i,nr
	character*(*) net(nr)
c
	nr=sar_nr
	do i=1,sar_nr
		net(i)=sar_net(i)
	enddo
	return
	end
c-----------------------------------------------------------------
c  return azimuths and distances
c-
	subroutine sar_getDistanceAzimuth(nr,dist,azi)
	include 'sardim.h'
	include 'sar.h'
	integer i,nr
	real azi(nr),dist(nr)
c
	nr=sar_nr
	do i=1,sar_nr
		dist(i)=sar_dis(i)
		azi(i)=sar_phi(i)
	enddo
	return
	end
	
