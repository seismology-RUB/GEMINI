c----------------------------------------------------------
c  root finding procedures:
c
c----------------------------------------------------------
c  bracket roots using root counting algorithm
c-------------------------------------------------------
	subroutine root_bracket(lmin,lmax,nval2,ctyp,iprint)
	include 'root.h'
	include 'nvmax.h'
	double precision lmin,lmax,wlfrac,vertadd,prm5,rstart,wlfrac_factor
	double precision ymin(nvmax)
	double precision elnew(nvall),elh(nvall),m5h(nvall),m5new(nvall)
	integer ifliq,nlay,irep,iprint,ibis,ilast,nbis,nlstrt
	integer ibismax,nval2,j,nrbr,krval
	integer knew(nvall),kh(nvall)
	character ctyp*1
c
	wlfrac=0.25d0
	vertadd=0.d0
	wlfrac_factor=1.d0
	irep=0
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
c
 103	call layer_set_stepsize_wlfrac(wlfrac*wlfrac_factor)
	if(ctyp.eq.'S') then
		call sphmot(lmin,nlstrt,rstart,ymin,kbr(1),vertadd,iprint,0,1)
		m5(1)=ymin(5)
		if(ifliq.eq.1) m5(1)=ymin(2)
		if(iprint.gt.1) print *,'Starting radius for lmin: ',rstart
		call sphmot(lmax,nlstrt,rstart,ymin,kbr(2),vertadd,iprint,0,1)
		m5(2)=ymin(5)
		if(ifliq.eq.1) m5(2)=ymin(2)
		if(iprint.gt.1) print *,'Starting radius for lmax: ',rstart
	else if(ctyp.eq.'T') then
		call tormot(lmin,nlstrt,rstart,ymin,kbr(1),vertadd,iprint,0,1)
		m5(1)=ymin(2)
		if(iprint.gt.1) print *,'Starting radius for lmin: ',rstart
		call tormot(lmax,nlstrt,rstart,ymin,kbr(2),vertadd,iprint,0,1)
		m5(2)=ymin(2)
		if(iprint.gt.1) print *,'Starting radius for lmax: ',rstart
	endif
c
	elbr(1)=lmin
	elbr(2)=lmax
	nrbr=kbr(2)-kbr(1)
	print *,'Bracket: ',nrbr,' roots counted in ',lmin,lmax,kbr(1),kbr(2)
c
c  case that no root was found
c-
	if(nrbr.eq.0) then
		nval2=0
		print *,'<root_bracket>: Return because no root was found'
		return
	endif
c
	nbis=nint(dlog(dble(nrbr))/dlog(2.d0))+5
	if(irep.eq.0) ibismax=2*nbis-1
	nval=1
	ibis=1
c
c  loop over bisections
c
c  bisect each existing interval if its root count > 1
c
 101	ilast=1
c
c  clean up once after enough bisections
c
	if(ibis.eq.nbis) then
		call root_cleanup(iprint)
	endif
c
 	do 10 j=1,nval
		krval=kbr(j+1)-kbr(j)
		elh(j)=0.5*(elbr(j)+elbr(j+1))
		prm5=m5(j+1)*m5(j)
c
c  do not bisect interval with no roots
c  or only one root and correct prm5
c
		if((krval.eq.0.and.prm5.gt.0.d0).or.(krval.eq.1.and.prm5.lt.0.d0)) then
			elh(j)=-1.d0
			goto 10
		endif
		if((krval.eq.0.and.prm5.lt.0.d0).or.
     1		   (krval.eq.1.and.prm5.gt.0.d0).or.(krval.lt.0)) then
			vertadd=vertadd+2.d0
			wlfrac_factor=wlfrac_factor*0.5
			print *,'Incorrect root count: use wlfrac = ',wlfrac*wlfrac_factor,
     1		        ' ,vertadd: ',vertadd
			goto 103
		endif
c
		if(ctyp.eq.'S') then
			call sphmot(elh(j),nlstrt,rstart,ymin,kh(j),vertadd,iprint,0,1)
			m5h(j)=ymin(5)
			if(ifliq.eq.1) m5h(j)=ymin(2)
		else if(ctyp.eq.'T') then
			call tormot(elh(j),nlstrt,rstart,ymin,kh(j),vertadd,iprint,0,1)
			m5h(j)=ymin(2)
		endif
 10	continue
c
c
c  reorder elbr and kbr
c
	elnew(1)=elbr(1)
	knew(1)=kbr(1)
	m5new(1)=m5(1)
	do 20 j=1,nval
	if(elh(j).gt.0.d0) then
		ilast=ilast+1
		elnew(ilast)=elh(j)
		knew(ilast)=kh(j)
		m5new(ilast)=m5h(j)
		ilast=ilast+1
		elnew(ilast)=elbr(j+1)
		knew(ilast)=kbr(j+1)
		m5new(ilast)=m5(j+1)
	else
		ilast=ilast+1
		elnew(ilast)=elbr(j+1)
		knew(ilast)=kbr(j+1)
		m5new(ilast)=m5(j+1)
	endif
 20	continue
c
	nval2=ilast-1
	if(nval2+1.gt.nvall) stop 'Bracket: Too many bisections for dimension'
c
c  copy knew and elnew onto kbr and elbr
c
	do 30 j=1,nval2+1
	m5(j)=m5new(j)
	elbr(j)=elnew(j)
 30	kbr(j)=knew(j)
c
	if(iprint.gt.0) then
	  print *,'Bisection: ',ibis,nbis,ibismax,nval
	  do j=1,nval2
	     write(6,'(f15.6,e16.6,2i7)') elbr(j),m5(j),kbr(j),kbr(j+1)-kbr(j)
	  enddo
	     write(6,'(f15.6,e16.6,i7)') elbr(nval2+1),m5(nval2+1),kbr(nval2+1)
	endif
c
c  check if we are done, that is no additional intervals required
c  reset layer_stepsize_wlfrac
c
	if(nval2.eq.nval) then
		call root_cleanup(iprint)
		call layer_set_stepsize_wlfrac(wlfrac)
		return
	endif
c
c  else continue with bisection
c
	nval=nval2
	ibis=ibis+1
	if(ibis.eq.ibismax) then
		wlfrac_factor=wlfrac_factor*0.5
		call layer_set_stepsize_wlfrac(wlfrac*wlfrac_factor)
		print *,'Bracket: redo bracketing with wlfrac = ',wlfrac*wlfrac_factor
		ibismax=ibismax+nbis/2
		irep=irep+1
		goto 103
	endif
	goto 101
c
	end
c-------------------------------------------------------
c  clean up list of roots
c-------------------------------------------------------
	subroutine root_cleanup(iprint)
	include 'root.h'
	integer j,je1,je2,krval,iprint,i
c
	j=2
	je1=1
	je2=0
 1	krval=kbr(j+1)-kbr(j)
	if(krval.eq.0.and.je1.eq.1.and.j.lt.nval) then
		j=j+1
		je1=j
		goto 1
	else if(krval.eq.0.and.je1.gt.1) then
		if(j.eq.nval) then
			je2=j
			goto 2
		else 
			j=j+1
			goto 1
		endif
	else if(krval.eq.1.and.je1.gt.1) then
		je2=j
		goto 2
	else if(krval.eq.2.and.je1.gt.1) then
		if(j.eq.nval) then
			je2=j
			goto 2
		else 
			j=j+1
			goto 1
		endif
	endif
c
c  normal case, next interval or done
c
	if(j.lt.nval) then
		j=j+1
		goto 1
	else
		goto 99
	endif
c
c  reorder list
c
 2	do i=0,nval-je2
		kbr(je1+i)=kbr(je2+i+1)
		elbr(je1+i)=elbr(je2+i+1)
		m5(je1+i)=m5(je2+i+1)
	enddo
c	print *,'Cleanup: List reordered',je1,je2
c
	if(je2.le.nval) then
		j=je1
		nval=nval-(je2-je1+1)
		je1=1
		je2=0
		goto 1
	endif		
 99	if(iprint.gt.1) print *,'Cleanup: ',(kbr(i),i=1,nval+1)
 	return
	end
c---------------------------------------------------------
c  refine bracketed eigendegree using bisection
c-------------------------------------------------------
	subroutine root_pindown(nrt,nov,elroot,m5root,ctyp,iprint,x1s)
	include 'nvmax.h'
	include 'root.h'
	integer maxit
	double precision eps
	parameter(MAXIT=50, eps=1.d-8)
	integer j,nrt,iprint,ifliq,i,nlay,knt,nov(nrt),nlstrt
	character ctyp*1
	double precision ymin(nvmax)
	double precision elroot(nrt),m5root(nrt),x1s(nrt)
	double precision xacc,dx,xmid,m5a,m5b,ela,elb,rstart,wlfrac,hdiv,vertadd,m5ref
c
	wlfrac=0.25d0
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	hdiv=1.d0
	vertadd=0.d0
	nrt=0
	do 10 j=nval,1,-1
		if(kbr(j+1)-kbr(j).eq.1) then
			call layer_set_stepsize_wlfrac(wlfrac*hdiv)
			nrt=nrt+1
 			ela=elbr(j)
			elb=elbr(j+1)
			m5a=m5(j)
			m5b=m5(j+1)	
			xacc=0.5d0*(ela+elb)*eps
			if(m5a*m5b.gt.0.d0) stop 'root must be bracketed in root_bracket'
			if(m5a.lt.0.d0) then
				elroot(nrt)=ela
				dx=elb-ela
			else 
				elroot(nrt)=elb
				dx=ela-elb
			endif
			m5ref=eps
			do i=1,MAXIT
				dx=.5d0*dx
				xmid=elroot(nrt)+dx
				if(ctyp.eq.'S') then
					call sphmot(xmid,nlstrt,rstart,ymin,knt,vertadd,iprint,0,0)
					if(ifliq.eq.0) then
						m5root(nrt)=ymin(5)
					else
						m5root(nrt)=ymin(2)
					endif
				else if(ctyp.eq.'T') then
					call tormot(xmid,nlstrt,rstart,ymin,knt,vertadd,iprint,0,0)
					m5root(nrt)=ymin(2)
				endif
				m5ref=max(m5ref,dabs(m5root(nrt)))
				if(iprint.gt.1) then
					print *,'Pindown: j,el,m5',j,xmid,m5root(nrt)
				endif
				if(m5root(nrt).le.0.d0) elroot(nrt)=xmid
				if(dabs(dx).lt.xacc.or.m5root(nrt).eq.0.d0) then
					m5root(nrt)=m5root(nrt)/m5ref
					x1s(nrt)=rstart
					nov(nrt)=-kbr(j+1)
					elroot(nrt)=xmid
					goto 10
				endif
			enddo
		stop 'too many bisections in pindown'
 10	continue
	end
c------------------------------------------------------------------
c   print result of root bracketing
c
	subroutine root_print
	include 'root.h'
	integer j
c
	print *,'RESULT OF ROOT BRACKETING'
	do j=1,nval
		write(6,'(f15.5,d15.6,2i7)') elbr(j),m5(j),kbr(j),kbr(j+1)-kbr(j)
	enddo
	j=nval+1
	write(6,'(f15.5,d15.6,i7)') elbr(j),m5(j),kbr(j)
	return
	end
