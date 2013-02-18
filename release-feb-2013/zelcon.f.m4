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
c-------------------------------------------------------------
c  compute complex elastic constants from real values
c  of ro,vpv,vph,vsv,vsh,eta,qk,qm
c--------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_elcon_type,`double complex')
	define(m4_zqk,ifelse(m4_elcon_type,`double complex',`dcmplx(1.d0,qk)',`1.d0'))
	define(m4_zqm,ifelse(m4_elcon_type,`double complex',`dcmplx(1.d0,qm)',`1.d0'))
c---------------------------------------------------------------
	subroutine zelcon(fratio,ro,alfv,alfh,betv,beth,eta,qk,qm,
     1		zpa,zpc,zpf,zpl,zpn,zkap,zmue)
	include 'pi.h'
	include 'zero.h'
	double precision fratio,ro,alfv,alfh,betv,beth,eta,qk,qm
	double precision epa,epc,epn,epl,epf,kappa,mue
	double precision da,dc,df,dl,dn,lamp2mue,lambda
	m4_elcon_type zpa,zpc,zpn,zpl,zpf,zkap,zmue,zlamp2mue,zlambda
c
	if(betv.gt.zero) then
		epa=ro*alfh*alfh
		epc=ro*alfv*alfv
		epn=ro*beth*beth
		epl=ro*betv*betv
		epf=eta*(epa-2.*epl)
c
c  isotropic equivalents
c
		kappa=(epc+4.*epa-4.*epn+4.*epf)/9.d0
		mue=(epc+epa+6.*epl+5.*epn-2.*epf)/15.d0
c		print *,'zelcon: ',ro,alfh,beth,alfv,betv,eta
c		print *,'zelcon: ',epa,epc,epf,epl,epn,kappa,mue
c
c  purely anisotropic parts
c
		lamp2mue=kappa+4.*mue/3.
		lambda=kappa-2.*mue/3.
		da=epa-lamp2mue
		dc=epc-lamp2mue
		df=epf-lambda
		dn=epn-mue
		dl=epl-mue
c
c  complex isotropic constants
c
		zkap=kappa*m4_zqk*dexp(dlog(fratio)*2.*datan(qk)/pi)
		zmue=  mue*m4_zqm*dexp(dlog(fratio)*2.*datan(qm)/pi)
		zlamp2mue=zkap+4.*zmue/3.
		zlambda=zkap-2.*zmue/3.
c
		zpa=zlamp2mue+da
		zpc=zlamp2mue+dc
		zpn=zmue+dn
		zpl=zmue+dl
		zpf=zlambda+df
	else
		kappa=ro*alfv*alfv
		zkap=kappa*m4_zqk
		zmue=0.d0
		zpa=zkap
		zpc=zkap
		zpf=zkap
		zpl=0.d0
		zpn=0.d0
	endif
c
	return
	end
