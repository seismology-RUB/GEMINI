c---------------------------------------------------------------
c  $Id:$
c
c  Calculate starting values for minor integration at the
c  surface (z=0). Treats both free surface and full-space models
c  with an homomgeneous spherical shell above the layer stack.
c
c  Never put a half space on top of a water layer!
c---------------------------------------------------------------
	subroutine stavsurf(zyst,nvar,rstart,csys,iprint)
 	integer nvar,iprint,nlay,fsflag
	double complex zyst(nvar),zom,zxa2,zeta,zxb2,zetb,za,zb
	double complex zpa,zpc,zpf,zpl,zpn,zkap,zmu
	double precision om,elp1,ro,el,rsq,rstart,vela,velb
	character csys*1
	common/degree/elp1
	common/omega/zom,om
c
	el=-0.5d0+dsqrt(0.25d0+elp1)
c
c  CASE 1: outer sphere exists
c          start integration downwards at rearth
c
	call earthmodel_getfsflag(fsflag)
	call earthmodel_getrearth(rstart)
	if(fsflag.eq.1) then
		call layer_getnlay(nlay)
		call layer_setnlactive(nlay+1)
		call earthmodel_elpar(rstart,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c
c  Spheroidal motion
c
		if(csys.eq.'M') then
			vela=dsqrt( real(zkap+4.d0*zmu/3.d0)/ro )
			velb=dsqrt( real(zmu)/ro )
			if(iprint.gt.0) then
				print *,'<stavsurf>: rturn-P = ',vela*dsqrt(elp1)/om
				print *,'<stavsurf>: rturn-S = ',velb*dsqrt(elp1)/om
			endif
			rsq=rstart*rstart
c
c  Ratio of spherical Hankel of the second kind functions
c
			zxa2=ro*zom*zom*rsq/(zkap+4.d0*zmu/3.d0)
			call sphan2cf(zsqrt(zxa2),el,zeta,iprint)
			zxb2=ro*zom*zom*rsq/zmu
			call sphan2cf(zsqrt(zxb2),el,zetb,iprint)
c
c Starting values for integration DOWNWARDS.
c
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
c
c  Toroidal motion
c
		else if(csys.eq.'T') then
			velb=dsqrt( real(zmu)/ro )
			rsq=rstart*rstart
			zxb2=ro*zom*zom*rsq/zmu
			call sphan2cf(zsqrt(zxb2),el,zetb,iprint)
			zyst(1) = 1.d0/el
			zyst(2) = zmu/(rstart*el)*(el-1.d0-zetb)
		endif
c
c  CASE 2: free surface
c          if top layer is a water layer:
c             for spheroidal case starting values are set separately in minorint
c             for toroidal case, starting values apply for solid boundary below water
c
	else
		if(csys.eq.'M') then
			zyst(1)=0.d0
			zyst(2)=1.d0
			zyst(3)=0.d0
			zyst(4)=0.d0
			zyst(5)=0.d0
		else if(csys.eq.'T') then
			zyst(1)=1.d0
			zyst(2)=0.d0
		endif
	endif
c
	return
	end 
