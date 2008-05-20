c-----------------------------------------------------------------------
c  Do root count using stored integration values from the
c  last call to bsint_driver

	subroutine rootcount(nvar)
	include 'nvmax.h'
	include 'bsstore.h'
	integer nvar,j

	do j=2,jstore
		call counter(nvar,xstore(j),ystore(1,j),ystore(1,j-1),dystore(1,j))
	enddo
	return
	end
c-------------------------------------------------
c  evaluate kount
c-------------------------------------------------
	subroutine counter(nvar,x,y,yold,dydx)
	double precision y(nvar),yold(nvar),dydx(nvar),fz,fnp,x
	integer rt_docount,rt_eif,rt_count,nvar
	common/roots/rt_docount,rt_eif,rt_count

c	print *,x,y(nvar),yold(nvar)
	if(sign(1.d0,yold(nvar))*sign(1.d0,y(nvar)).lt.0.d0) then
		if(nvar.eq.5) then
			fz=sign(1.d0,y(4)+y(3))
			fnp=sign(1.d0,-dydx(5)-2.d0/x*y(5))
		else if(nvar.eq.2) then
			fz=sign(1.d0,y(1))
			fnp=sign(1.d0,-dydx(2)-1.d0/x*y(2))
		endif
		if(fz*fnp.gt.0.d0) rt_count=rt_count-1
		if(fz*fnp.lt.0.d0) rt_count=rt_count+1
c		print *,'Counter: Count ',rt_count,x,fz,fnp,y(nvar)
	endif
c
	return
	end
