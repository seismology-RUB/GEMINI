c-------------------------------------------------------------------
c  $Id: earthmodel.h,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c  Declarations for class earthmodel that provides methods to
c  access earth model parameters
c
c  30/11/04: treat full-space models (earthmodel_fsflag: full-space flag)
c-------------------------------------------------------------------------
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	define(m4_elcon_type, `double complex')
c----------------------------------------------------------------------
	include 'nkk.h'
	integer earthmodel_nk,earthmodel_seldamp,earthmodel_aniflag,earthmodel_fsflag
	m4_elcon_type earthmodel_za(nkk),earthmodel_zc(nkk),earthmodel_zf(nkk)
	m4_elcon_type earthmodel_zl(nkk),earthmodel_zn(nkk)
	m4_elcon_type earthmodel_zkap(nkk),earthmodel_zmue(nkk)
	double precision earthmodel_rk(nkk),earthmodel_ro(nkk),earthmodel_fref
	double precision earthmodel_qk(nkk),earthmodel_qm(nkk)
	double precision earthmodel_qk2(nkk),earthmodel_qm2(nkk)
	double precision earthmodel_rearth
c
	m4_elcon_type earthmodel_za2(nkk),earthmodel_zc2(nkk),earthmodel_zf2(nkk)
	m4_elcon_type earthmodel_zl2(nkk),earthmodel_zn2(nkk)
	m4_elcon_type earthmodel_zkap2(nkk),earthmodel_zmue2(nkk)
	double precision earthmodel_ro2(nkk)
c
	common/earthmodel/earthmodel_rk,earthmodel_ro,earthmodel_qk,earthmodel_qm,
     1	   earthmodel_za,earthmodel_zc,
     1	   earthmodel_zf,earthmodel_zl,earthmodel_zn,earthmodel_zkap,earthmodel_zmue,
     2	   earthmodel_fref,earthmodel_rearth,earthmodel_nk,
     3	   earthmodel_seldamp,earthmodel_aniflag,earthmodel_fsflag
c
	common/emspline/earthmodel_ro2,earthmodel_qk2,earthmodel_qm2,
     1	   earthmodel_za2,earthmodel_zc2,
     1	   earthmodel_zf2,earthmodel_zl2,earthmodel_zn2,earthmodel_zkap2,earthmodel_zmue2
