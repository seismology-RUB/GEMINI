c-------------------------------------------------------
c  $Id: savemin.f,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c  copy contents of ysav after integration onto minor
c  and change to Woodhouse's definition
c-------------------------------------------------------
	subroutine savemin(js,sfl1,m)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	double precision m(6,nnd)
	double precision sfl1
	integer i,js,nsolid,n,nlay,ifliq
c
c  top layer if liquid
c
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	if(ifliq.eq.1) then
		nsolid=jwd
		do i=nnod,max(jwu,js),-1
			m(1,i)=ynod(1,i)
			m(2,i)=ynod(2,i)
		enddo
	else
		nsolid=nnod
	endif
c	
c  solid part
c
	do i=nsolid,js,-1
		m(1,i)=ynod(1,i)
		m(2,i)=sfl1*ynod(2,i)
		m(3,i)=sfl1*ynod(3,i)
		m(4,i)=sfl1*ynod(4,i)
		m(5,i)=sfl1*ynod(5,i)
		m(6,i)=-m(1,i)
c		print *,'Node: ',i,(m(n,i),n=1,6)
	enddo
c
c  below start radius
c
	do i=js-1,1,-1
	  do n=1,6
	    m(n,i)=0.d0
         enddo
	enddo
c
	return
	end
