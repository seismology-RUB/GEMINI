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
c  Declarations for class flnm that provides methods to
c  access earth model parameters
c
c  Model files have the following form
c  1st line: FLNODE
c  2nd line: some comment string
c  3rd line: flnm_fref flnm_seldamp flnm_aniflag
c            where flnm_aniflag: (0)=isotropic model, (1)=transversely anisotropic model
c            where fref is the reference frequency and seldamp is
c            a switch allowing to select among different modi how
c            attenuation is treated:
c            1 = change only imaginary part, no dispersion
c            2 = fully real moduli (not implemented)
c            3 = change of real and imaginary part of elastic moduli (not implemented)
c  4th line: number of the following entries
c  5th -   : entries giving constants with increasing depth
c            depth (m), density (g/cm^3), vp (km/s), vs (km/s), qkap, qmu
c
c  Fluid regions get negative qmu and beta=0.
c
c  30/11/04: treat full-space models (flnm_fsflag: full-space flag)
c-------------------------------------------------------------------------
	integer nkk1
	parameter(nkk1=100)
	integer flnm_nk,flnm_seldamp,flnm_aniflag,flnm_fsflag
	double precision flnm_qk(nkk1),flnm_qm(nkk1),flnm_fref,flnm_rearth
	double precision flnm_rk(nkk1),flnm_ro(nkk1),flnm_vpv(nkk1),flnm_vsv(nkk1)
	double precision flnm_vph(nkk1),flnm_vsh(nkk1),flnm_eta(nkk1)
	common/flnm/flnm_rk,flnm_ro,flnm_vpv,flnm_vsv,flnm_vph,flnm_vsh,flnm_eta,
     1	            flnm_qk,flnm_qm,flnm_fref,flnm_rearth,flnm_nk,flnm_seldamp,
     1	            flnm_aniflag,flnm_fsflag

