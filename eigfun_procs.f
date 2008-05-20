c-----------------------------------------------------------------
c   $Id: eigfun_procs.f,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c  Routines associated with class eigfun
c-----------------------------------------------------------------
c  print overview of results to screen 
c
	subroutine eigfun_printeigval(nrt,nov,p,cg,q,anorm,raylquo,m5root)
	integer nrt,n,nov(nrt)
	double precision p(nrt),cg(nrt),q(nrt),anorm(nrt),raylquo(nrt),m5root(nrt)
c
	write(6,'(a8,4a10,2a10,a10)') 'Overtone','Slowness','Phase','Group','Q','Norm','Wdiff','M5'
	do n=1,nrt
		write(6,'(i8,4f10.4,2e10.2,d10.2)') nov(n),p(n),1.d0/p(n),cg(n),
     1	                                    q(n),dsqrt(anorm(n)),raylquo(n),m5root(n)
	enddo
c
	return
	end
c---------------------------------------------------------------------
c  construct eigenfunction at a node from minors up to a sign
c  m(6,2): minors according to Woodhouse's definition (see savemin)
c  y: U,R,V,S,dU/dr,dV/dr
c  om: angular frequency in rad/s
c---------------------------------------------------------------------
	subroutine eigfun_create(elp1,nla,js,m1,m2,yeig)
	include 'nodesdim.h'
	include 'eigdim.h'
	include 'nodes.h'
	integer nrbot,nrtop,nlay,nl,nla,idmax,j,id,js
	double precision m1(6,nnd),m2(6,nnd),xmd(4)
	double precision yeig(nvx,nnd)
	double precision absxmd,x12,x13,x14,x23,x24,x34
	double precision div,sfl1,rr,s,om,elp1
	double precision zpa,zpc,zpf,zpl,zpn,zkap,zmu,ro,rop
	common/omega/om
c
	sfl1=dsqrt(elp1)
	call layer_getnlay(nlay)
	do 10 nl=nla,nlay
		call layer_setnlactive(nl)
		call layer_getbotnode(nl,nrbot)
		call layer_gettopnode(nl,nrtop)
		call layer_isfluid(nl,j)
		if(j.eq.1) goto 15
		do j=max(js,nrbot),nrtop
c
c  diagonal elements of X=M^(1) M^(2T) Sigma
c
			xmd(1)=-m1(2,j)*m2(3,j)+m1(3,j)*m2(2,j)
			xmd(2)=-m1(4,j)*m2(5,j)+m1(5,j)*m2(4,j)
			xmd(3)=-m1(2,j)*m2(4,j)+m1(4,j)*m2(2,j)
			xmd(4)=-m1(3,j)*m2(5,j)+m1(5,j)*m2(3,j)
			absxmd=dabs(xmd(1))
			div=absxmd
			idmax=1
			do id=2,4
				absxmd=dabs(xmd(id))
				if(absxmd.gt.div) then
					idmax=id
					div=absxmd
				endif
			enddo
c
c  now that maximum diagonal element is found compute
c  eigenfunction from it and its off-diagonals
c
			x12= m1(1,j)*m2(6,j)-m1(2,j)*m2(5,j)+m1(3,j)*m2(4,j)
			x13=-m1(1,j)*m2(2,j)+m1(2,j)*m2(1,j)
			x14=-m1(1,j)*m2(3,j)+m1(3,j)*m2(1,j)
			x23=-m1(1,j)*m2(4,j)+m1(4,j)*m2(1,j)
			x24=-m1(1,j)*m2(5,j)+m1(5,j)*m2(1,j)
			x34=-m1(2,j)*m2(5,j)+m1(4,j)*m2(3,j)+m1(6,j)*m2(1,j)
c
			s=sign(1.d0,xmd(idmax))
			div=dsqrt(div)
			yeig(idmax,j)=div
			if(idmax.eq.1) then
				yeig(2,j)=s*x12/div
				yeig(3,j)=s*x13/div
				yeig(4,j)=s*x14/div
			else if(idmax.eq.2) then
				yeig(1,j)=s*x12/div
				yeig(3,j)=s*x23/div
				yeig(4,j)=s*x24/div
			else if(idmax.eq.3) then
				yeig(1,j)=s*x13/div
				yeig(2,j)=s*x23/div
				yeig(4,j)=s*x34/div
			else if(idmax.eq.4) then
				yeig(1,j)=s*x14/div
				yeig(2,j)=s*x24/div
				yeig(3,j)=s*x34/div
			endif
c
c  change to usual definition in terms of U,R,V,S
c
			yeig(3,j)=yeig(3,j)/sfl1
			yeig(4,j)=yeig(4,j)/sfl1
c
c  fill y(5) with du/dr and y(6) with dv/dr
c
			call earthmodel_elpar(rnod(j),ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			rr=1.d0/rnod(j)
			yeig(5,j)=rr*zpf/zpc*(-2.d0*yeig(1,j)+sfl1**2*yeig(3,j))+yeig(2,j)/zpc
			yeig(6,j)=rr*(-yeig(1,j)+yeig(3,j))+yeig(4,j)/zpl

c			print *,nl,j,rnod(j),ro,zpa,zpl,div,
c			print*,'<ef_create:> ',nl,j,(yeig(id,j),id=1,nvx)
		enddo
		goto 10
c
c  fluid case, usual definition of eigenfunctions
c  take downward solution to construct eigenfunction
c  in the water layer.
c
 15		continue 
		call earthmodel_rholiq(rnod(j),ro,rop)
		do j=max(js,nrbot),nrtop
			if(nl.ne.nlay) then
				yeig(1,j)=m1(1,j)
				yeig(2,j)=m1(2,j)
			else
 				yeig(1,j)=m2(1,j)
				yeig(2,j)=m2(2,j)
			endif
			yeig(4,j)=0.d0
			call earthmodel_elpar(rnod(j),ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			rr=1.d0/rnod(j)
			yeig(3,j)=-yeig(2,j)/(rnod(j)*om*om*ro)
c
c  fill y(5) with dU/dr and y(6) with dV/dr
c
			yeig(5,j)=-2.d0*rr*yeig(1,j)+rr*sfl1**2*yeig(3,j)+yeig(2,j)/zpc
			yeig(6,j)=rr*yeig(1,j)-yeig(3,j)*(rr+rop/ro)
c			print *,nl,j,rnod(j),ro,zpa,zpl,(yeig(id,j),id=1,nvx)
		enddo
 10	continue
c
	return
	end
c-------------------------------------------------------------------
c  fix sign of eigenfunctions by looking for
c  the smoothest cubic spline interpolation of y_1
c  it works as follows:
c  Given y, yp at point i and y and yp at point i+1
c  up to one sign
c  we assume a cubic polynomial for y(r) and look which
c  sign of y(i+1) gives the smallest integrated curvature.
c-------------------------------------------------------------------
	subroutine eigfun_fixsign(nlstrt,js,yeig)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	double precision yeig(nvx,nnd)
	double precision s2,su,sv,cubl,cuwe,cvbl,cvwe,h,fsc
	integer nlstrt,nrbot,nl,nlay,j,nrbot,nrtop,i,jj,js,n
c
	call layer_getnlay(nlay)
	do nl=nlay,nlstrt,-1
		call layer_setnlactive(nl)
		call layer_gettopnode(nl,nrtop)
		call layer_getbotnode(nl,nrbot)
c
c  first ensure that y(1) has same sign above and below any interface
c
		if(nl.ne.nlay.and.yeig(1,nrtop+1)*yeig(1,nrtop).lt.0.d0) then
			do i=1,nvx
				yeig(i,nrtop)=-yeig(i,nrtop)
			enddo
		endif
c
c  make eigenfunction y(1) continuous across liquid-solid
c  interfaces. Note that we must multiply all eigenfunction
c  values down to rstart with the factor.
c
		i=0
		if(nl.lt.nlay) call layer_isfluid(nl+1,i)
		if(i.eq.1) then
			if(yeig(1,nrtop).eq.0.d0) then
				print *,'<eifsign>: eigenfunction U at top of water layer vanishes!'
				print *,' n = ',n
				stop
			endif
			fsc=yeig(1,nrtop+1)/yeig(1,nrtop)
			do jj=nrtop,js,-1
				do i=1,nvx
					yeig(i,jj)=fsc*yeig(i,jj)
				enddo
			enddo
		endif
c
c  jump to next layer if layer is liquid
c  because eigenfunction is known in this case
c
		call layer_isfluid(nl,i)
		if(i.eq.0) then
			do j=nrtop,max(js+1,nrbot+1),-1
c				h=rnod(j)-rnod(j-1)
				h=rnod(j-1)-rnod(j)
 				call eigfun_curv(yeig(1,j),yeig(5,j),yeig(1,j-1),yeig(5,j-1),
     1				                 h,cubl,cuwe)
				call eigfun_curv(yeig(3,j),yeig(6,j),yeig(3,j-1),yeig(6,j-1),
     1				                 h,cvbl,cvwe)
c
c  use curvatures with larger contrast to determine sign
c
				su=sign(1.d0,cuwe-cubl)
				sv=sign(1.d0,cvwe-cvbl)
				if(dabs(cubl-cuwe)/(cubl+cuwe).gt.dabs(cvbl-cvwe)/(cvbl+cvwe)) then
					s2=su
				else
					s2=sv
				endif
c
c  if cwe is less than cbl
c  change sign of eigenfunctions
c
				if(s2.lt.0.d0) then
					do i=1,nvx
						yeig(i,j-1)=-yeig(i,j-1)
					enddo
				endif
c				print *,'<eifsign:> ',nl,j,(yeig(i,j),i=1,6)
			enddo
		endif
	enddo
c
c  zero eigenfunctions below js
c
	do j=1,js-1
		do i=1,nvx
			yeig(i,j)=0.d0
		enddo
	enddo
c
	return
	end
c---------------------------------------------------------------------
c  compute integrated squared curvature of cubic polynomial
c  in interval with length h, given function values and derivatives
c  at ends of interval. Sign of y2 and y2p is unknown. Curvature
c  is evaluated for both cases
c  cbl: curvature with signs as is
c  cwe: curvature with changed signs
c---------------------------------------------------------------------
	subroutine eigfun_curv(y1,y1p,y2,y2p,h,cbl,cwe)
	double precision y1,y1p,y2,y2p,h,cbl,cwe
	double precision p,q,h2
c
	h2=h*h
c
c  integrated square of second derivative for signs as is
c
	p=y2-y1-h*y1p
	q=y2p-y1p
	cbl=3.d0*p*p-3.d0*h*p*q+h2*q*q
c
c  integrated square of second derivative for changed signs
c
	p=-y2-y1-h*y1p
	q=-y2p-y1p
	cwe=3.d0*p*p-3.d0*h*p*q+h2*q*q
c
	return
	end
c-----------------------------------------------------------------------
c  compute norm of eigenfunctions
c  group velocity, Q and Rayleigh quotient
c
	subroutine eigfun_integrals(elp1,nlstrt,js,yeig,ctyp,anorm,cg,qatt,raylquo)
	include 'eigfun.h'
	include 'nodesdim.h'
	integer i,j,nlstrt,js,nlay,nl,nrtop,nrbot
	double precision fint,anorm,cg,qatt,raylquo,om,rearth
	double precision elp1,yeig(nvx,nnd)
	character ctyp*1
	common/omega/om
c
	ef_typ=ctyp
	ef_elp1=elp1
	anorm=0.d0
	cg=0.d0
	qatt=0.d0
	raylquo=0.d0
	call layer_getnlay(nlay)
	do nl=nlstrt,nlay
		call layer_setnlactive(nl)
		call layer_gettopnode(nl,nrtop)
		call layer_getbotnode(nl,nrbot)
		do j=max(js+1,nrbot+1),nrtop
			do i=1,nvx
				ef_y1(i)=yeig(i,j-1)
				ef_y2(i)=yeig(i,j)
			enddo
			ef_node=j
			call eigfun_gauslv(fint,1)
			anorm=anorm+fint
			call eigfun_gauslv(fint,2)
			cg=cg+fint
			call eigfun_gauslv(fint,3)
			qatt=qatt+fint
			call eigfun_gauslv(fint,4)
			raylquo=raylquo+fint
		enddo
	enddo
c	print *,'integrals: ',anorm,cg,qatt,raylquo,' om = ',om
	call earthmodel_getrearth(rearth)
c
c  convert group velocity cg to km/s
c
	cg=cg/anorm*dsqrt(elp1)*rearth/om
	qatt=om*om/(qatt/anorm)
	raylquo=1.-raylquo/(om*om*anorm)

	return
	end
c------------------------------------------------------
c  compute the Frechet kernels for vp and vs and ro
c  at radial nodes
c  Use formulae of Dahlen & Tromp(p. 305) but
c  note that V(Dahlen)=sqelp1*V
c
c  requires change for an anisotropic medium !
c-------------------------------------------------------
	subroutine eigfun_frechet_sph(elp1,nlstrt,js,yeig,fkro,fkvp,fkvs,fkd)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'laydim.h'
	integer j,nl,js,nlstrt,nlay,nrtop,nrbot
	double precision om,elp1,akk,akm,akro,akd,alf,bet,ro,third
	double precision zpa,zpc,zpf,zpl,zpn,zkap,zmu
	double precision r,u,up,v,vp
	double precision yeig(nvx,nnd),fkro(nnd),fkvp(nnd),fkvs(nnd)
	double precision fkd(nlayer),fkdbot(nlayer),fkdtop(nlayer)
	common/omega/om
c
	third=0.333333333333333333d0
c
c  Loop over layers and nodes
c
	call layer_getnlay(nlay)
	do nl=nlstrt,nlay
		call layer_setnlactive(nl)
		call layer_gettopnode(nl,nrtop)
		call layer_getbotnode(nl,nrbot)
		do j=max(js,nrbot),nrtop
			r=rnod(j)
			u=yeig(1,j)
			v=yeig(3,j)
			up=yeig(5,j)
			vp=yeig(6,j)
			akk=(r*up+2.*u-elp1*v)**2
			akm=third*(2.*r*up-2.*u+elp1*v)**2
     1        	    +elp1*(r*vp-v+u)**2+elp1*(elp1-2.)*v*v
			akro=-(om*r)**2*(u*u+elp1*v*v)
			call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			alf=dsqrt((zkap+4.*third*zmu)/ro)
			bet=dsqrt(zmu/ro)
			fkro(j)=(zkap*akk+zmu*akm+ro*akro)/(2.*om*ro)
			fkvp(j)=2.*ro*alf*akk/(2.*om)
			fkvs(j)=2.*ro*bet*(akm-4.*third*akk)/(2.*om)
			akd = -zkap*akk-zmu*akm-ro*akro
     1			      +2.*zkap*r*up*(r*up-2.*u-elp1*v)
     1			      +4.*third*zmu*r*up*(2.*r*up-2.*u+elp1*v)
     2			      +2.*zmu*r*elp1*vp*(r*vp-v+u)
			if(j.eq.nrtop) fkdtop(nl) = akd
			if(j.eq.nrbot) fkdbot(nl) = akd
		enddo
	enddo
c
	do j=1,js-1
		fkro(j)=0.d0
		fkvp(j)=0.d0
		fkvs(j)=0.d0
	enddo
c
c  calculate frechet kernel for discontinuities
c  index nl refers to the discontinuity at the top of the layer !!
c-
	do nl=1,nlstrt-1
		fkd(nl)=0.d0
	enddo
	do nl=nlstrt,nlay-1
		fkd(nl)=(fkdbot(nl+1)-fkdtop(nl))/(2.*om)
	enddo
	fkd(nlay)=-fkdtop(nlay)/(2.*om)
c
	return
	end
c------------------------------------------------------
c  compute the Frechet kernels for vs
c  at radial nodes
c  Use formulae of Dahlen & Tromp(p. 305) but
c  note that V(Dahlen)=sqelp1*V
c
c  requires change for an anisotropic medium !
c-------------------------------------------------------
	subroutine eigfun_frechet_tor(elp1,nlstrt,js,yeig,fkro,fkvs,fkd)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'laydim.h'
	integer j,nl,js,nlstrt,nlay,nrtop,nrbot,ifliq
	double precision om,elp1,akro,akm,akd,bet,ro,third
	double precision zpa,zpc,zpf,zpl,zpn,zkap,zmu
	double precision r,w,wp
	double precision yeig(nvx,nnd),fkro(nnd),fkvs(nnd)
	double precision fkd(nlayer),fkdbot(nlayer),fkdtop(nlayer)
	common/omega/om
c
	third=0.333333333333333333d0
c
c  Loop over layers and nodes
c-
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	if(ifliq.eq.1) then
		call layer_gettopnode(nlay,nrtop)
		call layer_getbotnode(nlay,nrbot)
		do j=nrbot,nrtop
			fkro(j)=0.d0
			fkvs(j)=0.d0
		enddo
		fkd(nlay)=0.d0
		nlay=nlay-1
	endif	
	do nl=nlstrt,nlay
		call layer_setnlactive(nl)
		call layer_gettopnode(nl,nrtop)
		call layer_getbotnode(nl,nrbot)
		do j=max(js,nrbot),nrtop
			r=rnod(j)
			w=yeig(1,j)
			wp=yeig(3,j)
			akro=-(om*r)**2*elp1*w*w
			akm=elp1*(r*wp-w)**2+elp1*(elp1-2.)*w*w
			call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
			bet=dsqrt(zmu/ro)
			fkro(j)=(zmu*akm+ro*akro)/(2.*om*ro)
			fkvs(j)=2.*ro*bet*akm/(2.*om)
			akd = -zmu*akm-ro*akro+2.*zmu*r*elp1*wp*(r*wp-w)
			if(j.eq.nrtop) fkdtop(nl) = akd
			if(j.eq.nrbot) fkdbot(nl) = akd
		enddo
	enddo
c
	do j=1,js-1
		fkro(j)=0.d0
		fkvs(j)=0.d0
	enddo
c
c  calculate frechet kernel for discontinuities
c  index nl refers to the discontinuity at the top of the layer !!
c-
	do nl=1,nlstrt-1
		fkd(nl)=0.d0
	enddo
	do nl=nlstrt,nlay-1
		fkd(nl)=(fkdbot(nl+1)-fkdtop(nl))/(2.*om)
	enddo
	fkd(nlay)=-fkdtop(nlay)/(2.*om)
c
	return
	end
c------------------------------------------------------------------------
c  do a Gauss-Legendre integration over one node interval
c  isel selects from different integrands evaluated in intgds
c
	subroutine eigfun_gauslv(fint,isel)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'eigfun.h'
	integer isel,j,i
	double precision fint,val,val1,w(2),x(2)
	double precision rm,h,t1
	data w,x/.478628670499366d0,.236926885056189d0,
     +	         .538469310105683d0,.906179845938664d0/

	j=ef_node
	rm=.5d0*(rnod(j)+rnod(j-1))
	h=rnod(j)-rnod(j-1)
	call eigfun_intgds(isel,rm,val)
	fint=.568888888888889d0*val
	do i=1,2
		t1=x(i)*0.5d0*h
		call eigfun_intgds(isel,rm+t1,val)
		call eigfun_intgds(isel,rm-t1,val1)
		fint=fint+w(i)*(val+val1)
	enddo
	fint=fint*0.5d0*h
	return
	end
c------------------------------------------------------------------------
c  evaluates integrands for norm, group velocity
c  q and Rayleigh quotient at given radius
c  val(1): norm
c  val(2): group velocity
c  val(3): Q
c  val(4): Rayleigh quotient
c--------------------------------------------------
	subroutine eigfun_intgds(isel,r,val)
	include 'eigfun.h'
	integer isel
	double precision r,r2,val
	double precision u,up,v,vp,w,wp,ff,xf
	double precision zpa,zpc,zpf,zpl,zpn,ro,zkap,zmu,qk,qm
c
	r2=r*r
c
c  spheroidal integrals
c-
	if(ef_typ.eq.'S') then
c
c  first evaluate U,UP,V,VP at r by spline interpolation
c-
		call eigfun_interpol_sph(r,u,up,v,vp)
		ff=(2.d0*u-ef_elp1*v)/r
		xf=vp+(u-v)/r
c
c  get material parameters for r
c
		call earthmodel_q(r,qk,qm)
		call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c
c  normalization integrand = T
c
		if(isel.eq.1) then
			val=ro*(u*u+ef_elp1*v*v)*r2
c
c  group velocity integrand -> cI_1*cg
c
		else if(isel.eq.2) then
			val=zpl*u*u+zpa*ef_elp1*v*v-zpf*r*up*v-2.d0*(zpa-zpn)*u*v+zpl*u*(r*vp-v)
c
c  Q-integrand -> om**2*T*1/Q
c
		else if(isel.eq.3) then
			val=zkap*(r*up+r*ff)**2*qk+zmu*( (2.d0*r*up-r*ff)**2/3.d0+ef_elp1*r2*xf*xf
     1		                                     +ef_elp1*(ef_elp1-2.d0)*v*v )*qm
c
c  Rayleigh quotient, elastic energy
c
		else if(isel.eq.4) then
			val=zpc*up*up*r2+2.d0*zpf*r2*up*ff+(zpa-zpn)*r2*ff*ff
     1		        +zpl*ef_elp1*r2*xf*xf+ef_elp1*(ef_elp1-2.d0)*zpn*v*v
		endif
c
c  toroidal mode integrals
c
	else if(ef_typ.eq.'T') then
c
c  first evaluate W,WP at r by spline interpolation
c-
		call eigfun_interpol_tor(r,w,wp)
c
c  get material parameters for r
c
		call earthmodel_q(r,qk,qm)
		call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
c
c  normalization integrand = T
c-
		if(isel.eq.1) then
			val=ro*ef_elp1*w*w*r2
c
c  group velocity integrand -> cI_1*cg
c-
		else if(isel.eq.2) then
			val=ef_elp1*zpn*w*w
c
c  Q-integrand -> om**2*T*1/Q
c-
		else if(isel.eq.3) then
			xf=r*wp-w
			val=zmu*ef_elp1*( xf*xf+(ef_elp1-2.d0)*w*w )*qm
c
c  Rayleigh quotient, elastic energy
c-
		else if(isel.eq.4) then
			val=ef_elp1*(zpl*xf*xf+(ef_elp1-2.d0)*zpn*w*w)
		endif
	endif
c-
	return
	end
c------------------------------------------------------------------------
c  interpolate U,dU/dr,V,dV/dr betwen two nodes using a cubic polynomial
c  j is the index of the upper node
c  yh contains U,dU/dr,V,dV/dr
c
	subroutine eigfun_interpol_sph(r,u,up,v,vp)
	include 'eigfun.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'zero.h'
	integer j
	double precision h,h2,p,q,a2(2),a3(2),u,up,v,vp,t,r
c
	j=ef_node
	h=rnod(j)-rnod(j-1)
	if(dabs(h).lt.zero) then
		u=ef_y1(1)
		v=ef_y1(3)
		up=ef_y1(5)
		vp=ef_y1(6)
		return
	endif
	h2=h*h
c
c  quadratic and cubic coefficients for
c  spline interpolation of eigenfunctions
c
	p=ef_y2(1)-ef_y1(1)-h*ef_y1(5)
	q=ef_y2(5)-ef_y1(5)
	a2(1)=(3.d0*p-h*q)/h2
	a3(1)=(h*q-2.d0*p)/(h2*h)
	p=ef_y2(3)-ef_y1(3)-h*ef_y1(6)
	q=ef_y2(6)-ef_y1(6)
	a2(2)=(3.d0*p-h*q)/h2
	a3(2)=(h*q-2.d0*p)/(h2*h)
c
	t=r-rnod(j-1)
	u=ef_y1(1)+t*(ef_y1(5)+t*(a2(1)+t*a3(1)))
	v=ef_y1(3)+t*(ef_y1(6)+t*(a2(2)+t*a3(2)))
	up=ef_y1(5)+t*(2.d0*a2(1)+t*3.d0*a3(1))
	vp=ef_y1(6)+t*(2.d0*a2(2)+t*3.d0*a3(2))
c
	return
	end
c------------------------------------------------------------------------
c  interpolate W,dW/dr between two nodes using a cubic polynomial
c  j is the index of the upper node
c  yh contains W,dW/dr
c
	subroutine eigfun_interpol_tor(r,w,wp)
	include 'eigfun.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'zero.h'
	integer j
	double precision h,h2,p,q,a2,a3,w,wp,t,r
c
	j=ef_node
	h=rnod(j)-rnod(j-1)
	if(dabs(h).lt.zero) then
		w=ef_y1(1)
		wp=ef_y1(3)
		return
	endif
	h2=h*h
c
c  quadratic and cubic coefficients for
c  spline interpolation of eigenfunctions
c-
	p=ef_y2(1)-ef_y1(1)-h*ef_y1(3)
	q=ef_y2(3)-ef_y1(3)
	a2=(3.d0*p-h*q)/h2
	a3=(h*q-2.d0*p)/(h2*h)
c
	t=r-rnod(j-1)
	w=ef_y1(1)+t*(ef_y1(3)+t*(a2+t*a3))
	wp=ef_y1(3)+t*(2.d0*a2+t*3.d0*a3)
c
	return
	end
c------------------------------------------------------------------------
c   normalize eigenfunctions and zero below js
c
	subroutine eigfun_normalize(js,yeig,anorm)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	integer j,i,js
	double precision anorm,yeig(nvx,nnd)
c
	do j=1,js-1
		do i=1,nvx
			yeig(i,j)=0.d0
		enddo
	enddo
	do j=js,nnod
		do i=1,nvx
			yeig(i,j)=yeig(i,j)/dsqrt(anorm)
		enddo
	enddo
	return
	end
c-------------------------------------------------------------
c  fill dW/dr into third component of yeig
c
	subroutine eigfun_torderiv(nlstrt,js,yeig)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	integer nl,nlay,ifliq,nlstrt,js,nrbot,nrtop,j
	double precision yeig(nvx,nnd)
	double precision ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu
c
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	do nl=nlstrt,nlay
		call layer_setnlactive(nl)
		call layer_gettopnode(nl,nrtop)
		call layer_getbotnode(nl,nrbot)
		if(ifliq.eq.1.and.nl.eq.nlay) then
			do j=max(js,nrbot),nrtop
				yeig(3,j)=0.d0
			enddo
		else
			do j=max(js,nrbot),nrtop
				call earthmodel_elpar(rnod(j),ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
				yeig(3,j)=yeig(1,j)/rnod(j)+yeig(2,j)/zpl
			enddo
		endif
	enddo
	return
	end
