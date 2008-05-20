c---------------------------------------------------------------------
c  $Id: frechettor.f,v 1.2 2003/06/01 07:44:50 wolle Exp $
c
c                           FRECHETTOR
c
c  Compute partial derivatives of solution yout at re for perturbations of structure.
c	delta y(i,r_e,r_s)=int [G(i,n,r_e,r) (delta A(n,j,r)) y(j,r,r_s)] dr.
c	Here we compute the integrand in the above expression.
c
c  ipar:      index of parameter for which derivative is computed
c  ncomp:     number of components desired
c  comps:     components
c  green:     G(i,n,re,r)=green(n,i,nr)
c  yout:      y(j,r,rs)
c  pdev:      = G(i,n,r_e,r) (delta A(n,j,r)) y(j,r,r_s) (elpar/del)
c
c  Specification of inversion parameters is done in isoptbelpar.f
c  NOTE: kernels are only computed for the solidpart of the medium
c        pdev is set to zero in a top water layer
c---------------------------------------------------------------------------
	subroutine frechettor(ipar,ncomp,comps,green,yout,pdev,skd,iprint)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'layer.h'
	include 'zero.h'
	double complex green(2,2,nnd),yout(2,nnd),yh(2)
	double complex a(2,2),dela(2,2),pdev(nnd,2),skd(nlayer,2)
	double complex zpa,zpc,zpf,zpl,zpn,zkap,zmu
	double complex zpan,zpcn,zpfn,zpln,zpnn
	double precision r,ro,ron,ptb
	integer i,j,ico,ncomp,comps(ncomp),ipar,top,iprint
	integer jend,nr,nl,nlend
c
c  zero pdev at all nodes and comps
c
	do i=1,ncomp
		do nr=1,nnod
			pdev(nr,comps(i))=(0.d0,0.d0)
		enddo
	enddo
c
c  perturbation in solid part of the medium
c
	ptb=0.02
	jend=nnod
	if(layer_iflso(layer_nlay).eq.1) jend=jwd
	do nr=jstu,jend
		r=rnod(nr)
		call layer_getindex(r,nl,top)
		if(top.eq.1.and.dabs(rnod(nr-1)-rnod(nr)).lt.zero) nl=nl+1
		call layer_setnlactive(nl)
		call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
		call sysmattor(r,ro,zpa,zpc,zpf,zpl,zpn,a)
		call isoptbelpar(ipar,ro,zkap,zmu,ron,zpan,zpcn,zpfn,zpln,zpnn,ptb)
		call sysmattor(r,ron,zpan,zpcn,zpfn,zpln,zpnn,dela)
c
c  form (dela-a)/del
c
		do j=1,2
			do i=1,2
				dela(i,j)=(dela(i,j)-a(i,j))/ptb
			enddo
		enddo
c
c  delta A*yout
c
		do i=1,2
			yh(i)=0.d0
			do j=1,2
				yh(i)=yh(i)+dela(i,j)*yout(j,nr)
			enddo
		enddo
c
c  Frechet kernel for each component
c
		do i=1,ncomp
			ico=comps(i)
			pdev(nr,ico)=0.d0
			do j=1,2
				pdev(nr,ico)=pdev(nr,ico)+green(j,ico,nr)*yh(j)
			enddo
		enddo
		if(iprint.gt.0) then
			write(6,'(a,i5,4e15.5)') 'yout:  ',nr,(zabs(yout(i,nr)),i=1,2)
			write(6,'(a,i5,4e15.5)') 'yh:    ',nr,(zabs(yh(i)),i=1,2)
			write(6,'(a,i5,4e15.5)') 'green: ',nr,(zabs(green(i,comps(1),nr)),i=1,2)
			write(6,'(a,i5,e15.5)') 'pdev:  ',nr,zabs(pdev(nr,comps(1)))
		endif
	enddo
c
c  discontinuity kernels skd (do not consider bottom of a water layer and surface)
c  the layer index of skd refers to the top boundary of the layer
c  do not worry if layer boundary is below starting radius because
c  yout and green will have been set to zero there
c
	if(layer_iflso(layer_nlay).eq.1) then
		nlend=layer_nlay-2
	else
		nlend=layer_nlay-1
	endif
	do nl=1,nlend
		call layer_setnlactive(nl)
		call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
		call sysmattor(r,ro,zpa,zpc,zpf,zpl,zpn,a)
		call layer_setnlactive(nl+1)
		call earthmodel_elpar(r,ro,zpa,zpc,zpf,zpl,zpn,zkap,zmu)
		call sysmattor(r,ro,zpa,zpc,zpf,zpl,zpn,dela)
c
c  form jump of A
c
		do j=1,2
			do i=1,2
				dela(i,j)=dela(i,j)-a(i,j)
			enddo
		enddo
c
c  delta A*yout, take value of yout at layer top
c  should not matter because yout is continuous across discontinuities
c
		do i=1,2
			yh(i)=0.d0
			do j=1,2
				yh(i)=yh(i)+dela(i,j)*yout(j,layer_nrne(nl))
			enddo
		enddo
c
c  Frechet kernel for each component, skd is the negative jump of G*dA*yout
c
		do i=1,ncomp
			ico=comps(i)
			skd(nl,ico)=(0.d0,0.d0)
			do j=1,2
				skd(nl,ico)=skd(nl,ico)-green(j,ico,layer_nrne(nl))*yh(j)
			enddo
		enddo
	enddo
c
	return
	end
c---------------------------------------------------------------------
c                       SYSMATTOR
c
c  sytem matrix
c---------------------------------------------------------------------
	subroutine sysmattor(r,ro,zpa,zpc,zpf,zpl,zpn,a)
	double complex zom,zpa,zpc,zpf,zpl,zpn,a(2,2)
	double precision r,ro,elp1,om
	common/omega/zom,om
	common/degree/elp1
c
	a(1,1)=1.d0/r
	a(1,2)=1.d0/zpl
	a(2,1)=-zom*zom*ro-zpn/(r*r)*(2.d0-elp1)
	a(2,2)=-3.d0*a(1,1)
c
	return
	end
