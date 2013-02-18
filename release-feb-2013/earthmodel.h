c--------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
c-------------------------------------------------------------------
c  Declarations for class earthmodel that provides methods to
c  access earth model parameters
c
c  30/11/04: treat full-space models (earthmodel_fsflag: full-space flag)
c-------------------------------------------------------------------------
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	
c----------------------------------------------------------------------
	include 'nkk.h'
	integer earthmodel_nk,earthmodel_seldamp,earthmodel_aniflag,earthmodel_fsflag
	double complex earthmodel_za(nkk),earthmodel_zc(nkk),earthmodel_zf(nkk)
	double complex earthmodel_zl(nkk),earthmodel_zn(nkk)
	double complex earthmodel_zkap(nkk),earthmodel_zmue(nkk)
	double precision earthmodel_rk(nkk),earthmodel_ro(nkk),earthmodel_fref
	double precision earthmodel_qk(nkk),earthmodel_qm(nkk)
	double precision earthmodel_qk2(nkk),earthmodel_qm2(nkk)
	double precision earthmodel_rearth
c
	double complex earthmodel_za2(nkk),earthmodel_zc2(nkk),earthmodel_zf2(nkk)
	double complex earthmodel_zl2(nkk),earthmodel_zn2(nkk)
	double complex earthmodel_zkap2(nkk),earthmodel_zmue2(nkk)
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
