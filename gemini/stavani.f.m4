c---------------------------------------------------------------
c  Calculates STARTING VALUES for integration of the spheroidal
c  motion upwards to the source level.
c  The values are derived from a homogeneous isotropic Earth model
c  assumed beneath the starting radius r1. The results are
c  spherical Bessel functions. To avoid explicit calculation of
c  those functions, we use ratios of them. These ratios can 
c  be computed from a recurrence relation. See subroutine ZETL().
c  Author: J.R. Dalkolmo, Modified: W. Friederich
c---------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_function_type,`double complex')
	define(m4_elcon_type,`double complex')
	define(m4_abs,ifelse(m4_function_type,`double complex',`zabs($*)',`dabs($*)'))
	define(m4_sqrt,ifelse(m4_function_type,`double complex',`zsqrt($*)',`dsqrt($*)'))
c---------------------------------------------------------------
	subroutine stavani(zyst,nvar,rstart,nla,csys,iprint)
 	integer nla,ifliq,nvar,iprint
	m4_function_type zyst(nvar),zom,zxa2,zeta,zxb2,zetb,za,zb
	m4_elcon_type zpa,zpc,zpf,zpl,zpn,zkap,zmu
	double precision om,elp1,rstart,ro,el,rsq
	character csys*1
	common/degree/elp1
	common/omega/zom,om
c
c  Get equivalent isotropic moduli at the starting radius
c
	call layer_setnlactive(nla)
	call earthmodel_elpar(rstart,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c
	call layer_isfluid(nla,ifliq)
	el=-0.5d0+dsqrt(0.25d0+elp1)
	rsq=rstart*rstart
c	print *,'stavani: rstart,nla: ',rstart,nla
c
c  Spheroidal motion
c
	if(csys.eq.'M') then
		zxa2=ro*zom*zom*rsq/(zkap+4.d0*zmu/3.d0)
c
c  Ratio of spherical Bessel functions
c
		if(m4_abs(zxa2/elp1).le.1.d0) then
			call zetl(zxa2,el,zeta,iprint)
		else
		ifelse(m4_function_type,`double complex',`
			call sphancf(m4_sqrt(zxa2),el,zeta,iprint)
		',`
			call zetl(zxa2,el,zeta,iprint)
		')
		endif
c		print *, '<stavani>: ',zxa2,ro,zkap+4d0*zmu/3.d0,rsq,zeta
		if(ifliq.eq.0) then
			zxb2=ro*zom*zom*rsq/zmu
			if(m4_abs(zxb2/elp1).le.1.d0) then
				call zetl(zxb2,el,zetb,iprint)
			else
			ifelse(m4_function_type,`double complex',`
				call sphancf(m4_sqrt(zxb2),el,zetb,iprint)
			',`
				call zetl(zxb2,el,zetb,iprint)
			')
			endif
c			print *, '<stavani>: ',zxb2,ro,zmu,rsq,zetb
		endif
c
c Starting values for integration UPWARDS in liquid or solid layer.
c
		if(ifliq.eq.1)then
			zyst(1) = -(1.d0-zeta/el)/(zom*zom*ro*rstart)
			zyst(2) = 1.d0/el
			return 
		else
			za=zeta/el
			zb=zetb/(el+1.d0)
			zyst(2)=(-za+zb*(za-1.d0))/el
			zyst(1)=zmu/rstart*(zxb2/el+2.d0*elp1*zyst(2))
			zyst(3)=zmu/rstart*(-2.d0*zyst(2)+zxb2/elp1*(za-1.d0))
			zyst(4)=-zmu/rstart*(zxb2/(el*el)*(1.d0-zb)+4.d0*zyst(2))
			zyst(5)=zmu*zmu/rsq*(-4.d0*(el-1.d0)*(el+2.d0)*zyst(2)
     1			        +zxb2/el*(zxb2/elp1-2.d0*(el-1.d0)*(2.d0*el+1.d0)/elp1
     2			        -4.d0/(el+1.d0)*za-2.d0/el*zb))
c			print *,(zyst(i),i=1,5)
		endif
c
c  Toroidal motion
c
	else if(csys.eq.'T') then
		zxb2=ro*zom*zom*rsq/zmu
		if(m4_abs(zxb2/elp1).le.1.d0) then
			call zetl(zxb2,el,zetb,iprint)
		else
		ifelse(m4_function_type,`double complex',`
			call sphancf(m4_sqrt(zxb2),el,zetb,iprint)
		',`
			call zetl(zxb2,el,zetb,iprint)
		')
		endif
		zyst(1) = 1.d0/el
		zyst(2) = zmu/(rstart*el)*(el-1.d0-zetb)
	endif
c
	return
	end 
