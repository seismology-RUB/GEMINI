c-----------------------------------------------------------------------
c  Routines of class bsint. Used for integration of
c  a system of ordinary differential equations
c
c  NOTES:
c     tiny increased to 1.e-10 (larger than precision) to avoid
c     a series of identical xstore's leading to NANs in nodes_getsolution.
c     This happened when x-x2 was just smaller than the previously set 1.e-12
c     and another step was begun.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c  Bulirsch-Stoer driver routine
c
c  ystart:     initial values of function vector
c  nvar:       number of components of function vector
c  x1:         start of integration interval
c  x2:         end of integration interval
c  h1:         hint for stepsize to be used
c  sysmat:     name of subroutine which evaluates system matrix
c----------------------------------------------------------------------
c  m4 macros to set type of function vector and correct abs-intrinsic
c
	define(m4_function_type, `double complex')
	define(m4_abs,ifelse(m4_function_type,`double complex',`zabs($*)',`dabs($*)'))
c----------------------------------------------------------------------
	subroutine bsint_driver(ystart,nvar,x1,x2,h1,sysmat)
	include 'bsint.h'
	include 'nvmax.h'
	include 'bsstore.h'
	m4_function_type ystart(nvar),y(nvmax)
	double precision x,x1,x2,h,h1,hdid,hnext
	integer i,nvar,nstp,nsub
	external sysmat
c
	x=x1
	h=sign(h1,x2-x1)
	do i=1,nvar
		y(i)=ystart(i)
	enddo
c
c  integration loop
c  jstore is the index of the ystore-array at which new values will
c  be stored (one after the last element). Here initialized to 1.
c
	jstore=1
	do nstp=1,maxstep
c
c If step can overshoot end, cut down stepsize
c
		if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
		if(dabs(h).lt.tiny) return
		if(jstore+maxsub.gt.nnstore) then
			print *,'<bsint_driver>: More stored function values than dimensioned'
			stop
		endif
		call bsint_step(y,nvar,x,h,hdid,hnext,nsub,
     1		                ystore(1,jstore),dystore(1,jstore),sysmat)
c
c  Store values from xm(0) to xm(nsub-1). Omit xm(nsub) because it is identical
c  to xm(0) of the next step. Exception is made for the last step.
c  jstore is increased by nsub and points to one element after the last stored.
c
		do i=0,nsub-1
			xstore(i+jstore)=x-hdid+i*hdid/nsub
		enddo
		jstore=jstore+nsub
c
c  Are we done ?
c  jstore points to one element after the last stored and that's where
c  we store the solution at the end of the interval.
c  On return, jstore points to the last element.
c
		if(dabs(x-x2).lt.tiny) then
			xstore(jstore)=x2
			do i=1,nvar
				ystart(i)=y(i)
				ystore(i,jstore)=y(i)
			enddo
			return
		endif
c
c  Check against minimum stepsize
c
		if(dabs(hnext).lt.hmin) then
			print *,'nvar, h, x = ', nvar,h,x
			stop ' stepsize smaller than min'
		endif
		h=hnext
	enddo
c
c  something went wrong
c
	stop ' too many steps '
	end
c-----------------------------------------------------------------
c             Bulirsch-Stoer STEP
c
c  Taken and slightly modified from:
c  Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c  1986, Cambridge University Press.
c------------------------------------------------------------------------
	subroutine bsint_step(y,nvar,x,htry,hdid,hnext,nsub,ym,dym,sysmat)
	include 'bsint.h'
	include 'nvmax.h'
	external sysmat
	m4_function_type y(nvar),yerr(nvmax),ysav(nvmax),yseq(nvmax)
	m4_function_type ym(nvmax,0:maxsub),dym(nvmax,0:maxsub)
	double precision h,htry,x,xsav,hdid,hnext,xest,errmax
	integer nseq(imax),nvar,i,j,nsub
	data nseq/2,4,6,8,12,16,24/
	save nseq
c
	h=htry
	xsav=x
	do i=1,nvar
		ysav(i)=y(i)
	enddo
c
c  loop over subdivisions 
c  get accuracy for component with greatest modulus not for all of them
c
 1	do i=1,imax
		call bsint_mmid(ysav,nvar,xsav,h,nseq(i),i,yseq,ym,dym,sysmat)
		xest=(h/nseq(i))**2
		call bsint_rzextr(i,xest,yseq,y,yerr,nvar)
c		if(i.le.imin) write(6,'(a,i5,7d11.4)') 'bsint_step: ',i,xsav,h/nseq(i),(y(j),j=1,nvar)
c
c  Guard against spurious early convergence.
c
		if(i.gt.imin) then
	      	errmax=0.d0
			do j=1,nvar
				errmax=max( errmax,m4_abs(yerr(j))/m4_abs(y(j)) )
			enddo
c			write(6,'(a,i5,8d11.4)') 'bsint_step: ',i,xsav,h/nseq(i),(y(j),j=1,nvar),errmax
			errmax=errmax/eps
c
c  Cancel and new step if accuracy is sufficient
c
			if(errmax.lt.1.d0) then
				x=x+h
				hdid=h
				if(i.eq.nuse) then
					hnext=h*shrink
				else if(i.eq.nuse-1) then
					hnext=h*grow
				else
					hnext=(h*nseq(nuse-1))/nseq(i)
				endif
				nsub=nseq(i)
				return
			endif
		endif
	enddo
c
c  subdivisions were not sufficient
c
	h=0.25d0*h/2**((imax-nuse)/2)
c	print *,' step size diminished to ',h,' at point ',x
	if(x+h.eq.x) then
		print *, 'step size underflow.'
		print *,x,h
		stop
	endif
	goto 1
	end
c----------------------------------------------------------------------
c             Modified MIDpoint step
c
c  construct solution using an overlapping scheme:
c  y(xs+h)=y(x)+h*yp(xs)
c  y(xs+2*h)=y(x)+2*h*yp(xs+h)
c  y(xs+3*h)=y(xs+h)+2*h*yp(xs+2h)
c  y(xs+4*h)=y(xs+2*h)+2*h*ys(xs+3*h) and so on
c  stored value at x=xs+nh: ys=0.25*(y(x-h)+2*y(x)+y(x+h))
c----------------------------------------------------------------------
	subroutine bsint_mmid(y,nvar,xs,htot,nstep,isub,yout,ys,dys,sysmat)
	include 'bsint.h'
	include 'nvmax.h'
	external sysmat
	m4_function_type y(nvar),yout(nvar),swap,dy(nvmax)
	m4_function_type ym(nvmax),yn(nvmax),ys(nvmax,0:maxsub),dys(nvmax,0:maxsub)
	m4_function_type asys(nvmax,nvmax,0:48),rhs(nvmax,0:48)
	double precision xs,htot,h,x
	integer nvar,nstep,i,m,nelsm(imax),elsm(8,imax),isub,j
	data nelsm/3,2,4,4, 4, 8, 8/
	data elsm/   0,24,48, 0, 0, 0, 0, 0,
     1	            12,36, 0, 0, 0, 0, 0, 0,
     1	             8,16,32,40, 0, 0, 0, 0,
     1	             6,18,30,42, 0, 0, 0, 0,
     1	             4,20,28,44, 0, 0, 0, 0,
     1	             3, 9,15,21,27,33,39,45,
     1	             2,10,14,22,26,34,38,46/
	save nelsm,elsm,asys,rhs
c
c  compute new elements of sysmat needed for this subdivision
c  of the interval
c
c	write(6,'(a,7d11.4)') 'm=    0',xs,htot,(y(j),j=1,nvar)
	do j=1,nelsm(isub)
		x=xs+elsm(j,isub)*htot/48.d0
		call sysmat( x,asys(1,1,elsm(j,isub)),rhs(1,elsm(j,isub)) )
	enddo
c
	h=htot/nstep
	do i=1,nvar
		ym(i)=y(i)
	enddo
	do i=1,nvar
		dy(i)=rhs(i,0)
		do j=1,nvar
			dy(i)=dy(i)+asys(i,j,0)*ym(j)
		enddo
	enddo
	do i=1,nvar
		yn(i)=ym(i)+h*dy(i)
		ys(i,0)=ym(i)
		dys(i,0)=dy(i)
	enddo
c	write(6,'(a,7d11.4)') 'm=    0',xs,h,(dy(i),i=1,nvar)
c	write(6,'(a,7d11.4)') 'm=    0',xs,h,(yn(i),i=1,nvar)
	do m=1,nstep-1
		do i=1,nvar
			dy(i)=rhs(i,m*48/nstep)
			do j=1,nvar
				dy(i)=dy(i)+asys(i,j,m*48/nstep)*yn(j)
			enddo
		enddo
		do i=1,nvar
			swap=ym(i)+2.d0*h*dy(i)
			ym(i)=yn(i)
			yn(i)=swap			
			ys(i,m)=ym(i)
			dys(i,m)=dy(i)
		enddo
c		write(6,'(a,i5,7d11.4)') 'm=',m,xs+m*h,h,(dy(i),i=1,nvar)
c		write(6,'(a,i5,7d11.4)') 'm=',m,xs+m*h,h,(yn(i),i=1,nvar)
	enddo
	do i=1,nvar
		dy(i)=rhs(i,48)
		do j=1,nvar
			dy(i)=dy(i)+asys(i,j,48)*yn(j)
		enddo
	enddo
	do i=1,nvar
		yout(i)=0.5*(yn(i)+ym(i)+h*dy(i))
		ys(i,nstep)=yout(i)
		dys(i,nstep)=dy(i)
	enddo
c	write(6,'(a,i5,7d11.4)') 'm=',m,xs+m*h,h,(dy(i),i=1,nvar)
c	write(6,'(a,i5,7d11.4)') 'm=',m,xs+m*h,h,(yout(i),i=1,nvar)
c
	return
	end
c-----------------------------------------------------------------
c                     R Z E X T R
c
c        Rational function extrapolation 
c
c   Taken from:
c   Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c   1986, Cambridge University Press.
c-----------------------------------------------------------------
	subroutine bsint_rzextr(iest,xest,yest,yz,dy,nv)
	include 'nvmax.h'
	include 'bsint.h'
	integer j,nv,k,m1,iest
	m4_function_type yest(nv),yz(nv),dy(nv),d(nvmax,imax),v,c,b,b1,yy,ddy
	double precision x(imax),fx(imax),xest
	save
c
	x(iest)=xest
	if(iest.eq.1) then
c		write(6,'(a,2i4,f10.3)') 'rzextr: ',iest,iest,x(iest)
		do 11 j=1,nv
			yz(j)=yest(j)
			d(j,1)=yest(j)
			dy(j)=yest(j)
 11		continue
	else
		m1=iest
c		write(6,'(a,2i4,f10.3)') 'rzextr: ',m1,iest,x(iest)
		do 12 k=1,m1-1
			fx(k+1)=x(iest-k)/xest
c			print *,'fx: ',k,x(iest-k),fx(k+1)
 12		continue
c
c...................extrapolation part
c
		do 14 j=1,nv
			yy=yest(j)
			v=d(j,1)
			c=yy
			d(j,1)=yy
			do 13 k=2,m1
				b1=fx(k)*v
				b=b1-c
				if(b.ne.0.d0) then
					b=(c-v)/b
					ddy=c*b
					c=b1*b
				else
					ddy=v
				endif
				if(k.ne.m1) v=d(j,k)
				d(j,k)=ddy
				yy=yy+ddy
 13			continue
			dy(j)=ddy
			yz(j)=yy
 14		continue
	endif
	return
	end
c----------------------------------------------------------
c  extract solution at equidistant nodes
c
c  nvar:	number of components of function vector
c  x1:		begin of integration interval
c  x2:		end of integration interval
c  dx:		desired sampling
c  nmin:	minimum number of samples
c  nx:		number of created nodes
c  x:		x-values of nodes
c  y:		value of function vector at nodes
c----------------------------------------------------------
	subroutine bsint_sample_solution(nvar,x1,x2,dx,nmin,nx,x,y)
	include 'nvmax.h'
	include 'bsstore.h'
	double precision x1,x2,dx,x(nx),a,b,dxx
	m4_function_type y(nvmax,nx)
	integer nvar,nmin,nx,i,j,k
c
	dxx=dx
	nx=abs(nint((x2-x1)/dxx))
	dxx=(x2-x1)/nx
	if(nx.lt.nmin) then
		dxx=(x2-x1)/nmin
		nx=nmin
	endif
c
	nx=nx+1
	do i=1,nx
		x(i)=x1+(i-1)*dxx
		call dlocate(x(i),jstore,xstore,k)
		a=(xstore(k+1)-x(i))/(xstore(k+1)-xstore(k))
		b=1.d0-a
		do j=1,nvar
			y(j,i)=a*ystore(j,k)+b*ystore(j,k+1)
		enddo
	enddo
c
	return
	end
