c----------------------------------------------------------------------
c  $Id: sar.h,v 1.2 2003/04/03 12:09:38 wolle Exp $
c
c  Header file sar.h related to sar_procs.f
c  Declares variables associated with source-and-receiver (sar) information
c  Use requires include of sardim.h
c----------------------------------------------------------------------
	character sar_eventid*13,sar_name(nstatt)*5,sar_csys*1,sar_comp(nstatt)*1,sar_net(nstatt)*2
	integer sar_year,sar_month,sar_day,sar_hour,sar_minute
	integer sar_istyp,sar_nr,sar_nrorg,sar_lfdnr(nstatt)
	integer sar_zlim_flag,sar_netflag
	real sar_sec,sar_cendt,sar_xsrc,sar_ysrc,sar_depsrc
	real sar_zmin,sar_zmax
	real sar_mt(6),sar_force(3),sar_hdur,sar_tm_ec_gc(3,3)
	real sar_xsta(nstatt),sar_ysta(nstatt),sar_dis(nstatt),sar_phi(nstatt)
c
	common/sar/sar_xsta,sar_ysta,sar_mt,sar_force,sar_tm_ec_gc,
     1	     sar_sec,sar_cendt,sar_xsrc,sar_ysrc,sar_depsrc,sar_hdur,
     2	     sar_dis,sar_phi,sar_zmin,sar_zmax,
     2	     sar_istyp,sar_nr,sar_nrorg,sar_lfdnr,sar_zlim_flag,sar_netflag,
     2	     sar_year,sar_month,sar_day,sar_hour,sar_minute,
     3	     sar_eventid,sar_name,sar_csys,sar_comp,sar_net
