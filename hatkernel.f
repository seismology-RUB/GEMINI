c-------------------------------------------------------------------------------
c  $Id: hatkernel.f,v 1.3 2003/06/01 07:44:50 wolle Exp $
c
c  compute change of DSV for each node of the earth model. Assume a hat function
c  perturbation for each node with relative perturbation of 1. Let rk(k) be
c  the target earthmodel node. Then we integrate pdev(r)*hat(r) from rk(k-1)
c  to rk(k+1) and hat(rk(k-1))=hat(rk(k+1))=0. If rk(k) is a on discontinuity
c  a half hat is used on the appropriate side.
c  We return sk(model-node,comp) containing the hat kernels. Multiplication
c  with the relative perturbations at the earthmodel nodes and summing over
c  the nodes results in the change of the DSV for given kf.
c  Note pdev(r) has jumps at source and receiver!
c--------------------------------------------------------------------------------
	subroutine hatkernel(ncomp,comps,pdev,sk,iprint)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'earthmodel.h'
	include 'layer.h'
	include 'zero.h'
	integer ncomp,comps(ncomp)
	integer nl,k,i,nr1,nr2,nr,l,ico,iprint
	double precision r1,r2,rm,h1,h2,w(2),dr
	double complex pdev(nnd,4),sk(nkk,4)
	double precision hat
c
	if(iprint.gt.0) print *,'            HATKERNEL' 
c
c  Loop over layers
c
	do nl=1,layer_nlay
c
c  Loop over model nodes in layer
c  if we are in the halfspace, r1 = deepest Green node
c  if we are at the layer bottom rk(k-1) has the same value as rk(k)
c  and does not contribute to the integral
c  if rk(k) is the surface we set rk(k+1) to rb(nlay)
c
		do k=layer_iktop(nl-1)+1,layer_iktop(nl)
c
c  zero sk for this hat
c
			do i=1,ncomp
				sk(k,comps(i))=(0.d0,0.d0)
			enddo
c
c  set limits for integration over the hat
c
			if(nl.eq.1.and.k.eq.1) then
				r1=rnod(1)
			else
				r1=earthmodel_rk(k-1)
			endif
			rm=earthmodel_rk(k)
			if(k.eq.layer_iktop(layer_nlay)) then
				r2=layer_rb(layer_nlay)
			else
				r2=earthmodel_rk(k+1)
			endif
c
c  find Green nodes inside the hat but use Green nodes
c  at bottom and top boundary of layer
c  check for double green node
c
			if(k.eq.layer_iktop(nl-1)+1) then
				nr1=layer_nrna(nl)
			else
				call dlocate(r1,nnod,rnod,nr1)
				nr1=nr1+1
 1				if(rnod(nr1).le.r1) then
					nr1=nr1+1
					goto 1
				endif
			endif
c
c  make sure that we get an inner node nr2 if r2=rb(nl) but rm<rb(nl)
c
			if(k.eq.layer_iktop(nl)) then
				nr2=layer_nrne(nl)
			else
				if(dabs(r2-layer_rb(nl)).lt.zero) then
					nr2=layer_nrne(nl)-1
				else				
					call dlocate(r2,nnod,rnod,nr2)
				endif
			endif
c
			if(iprint.gt.0) then
				print *,'nl,k,r1,rm,r2= ',nl,k,r1,rm,r2
				print *,'nr1,nr2=',nr1,nr2
				print '(a,2d20.12)','rnod(nr1),rnod(nr2)=',rnod(nr1),rnod(nr2)
			endif
c
c  we need at least two green nodes in each hat
c
			if(nr2-nr1+1.lt.2) then
				print *,'Not enough Green nodes in layer: ',nl,k,nr1,nr2
				stop
			endif
c
			h1=rm-r1
			h2=r2-rm
c
c  sum up hat*pdev in the hat
c  double Green nodes do not hurt because dr=0
c
			do nr=nr1,nr2-1
				do l=nr,nr+1
					w(l-nr+1)=hat(rnod(l),r1,rm,r2,h1,h2)
				enddo
				dr=rnod(nr+1)-rnod(nr)
				do i=1,ncomp
					ico=comps(i)
					sk(k,ico)=sk(k,ico)+0.5*dr*(w(1)*pdev(nr,ico)+w(2)*pdev(nr+1,ico))
				enddo
			enddo
c
c  add rand terms: nr=nr1 and nr=nr2
c  do not add left rand term if rm=rb(nl-1)
c
			if(k.ne.layer_iktop(nl-1)+1) then
				w(1)=hat(rnod(nr1),r1,rm,r2,h1,h2)
			else
				w(1)=0.
			endif
c
c  do not add right rand term if rm=rb(nl)
c
			if(k.ne.layer_iktop(nl)) then
				w(2)=hat(rnod(nr2),r1,rm,r2,h1,h2)
			else
				w(2)=0.
			endif
c
			do i=1,ncomp
				ico=comps(i)
				sk(k,ico)=sk(k,ico)+0.5*(rnod(nr1)-r1)*w(1)*pdev(nr1,ico)
     1				                   +0.5*(r2-rnod(nr2))*w(2)*pdev(nr2,ico)
			enddo
			if(iprint.gt.0) print *,nl,k,(sk(k,comps(i)),i=1,ncomp)
		enddo
	enddo
	return
	end
c-------------------------------------------------------------------
c   evaluate hat function
c
	double precision function hat(r,r1,rm,r2,h1,h2)
	double precision r,r1,rm,r2,h1,h2
c
	if(r.lt.rm) then
		hat=(r-r1)/h1
	else if(r.gt.rm) then
		hat=(r2-r)/h2
	else
		hat=1.d0
	endif
	return
	end
