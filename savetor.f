c-------------------------------------------------------
c  copy contents of ynod after integration onto yeig
c-------------------------------------------------------
	subroutine savetor(js,yeig)
	include 'eigdim.h'
	include 'nodesdim.h'
	include 'nodes.h'
	include 'nodesstore.h'
	double precision yeig(nvx,nnd)
	integer i,js,nsolid,n,nlay,ifliq
c
c  top layer if liquid
c
	call layer_getnlay(nlay)
	call layer_isfluid(nlay,ifliq)
	if(ifliq.eq.1) then
		nsolid=jwd
	else
		nsolid=nnod
	endif
c
	do i=nnod,nsolid+1,-1
		yeig(1,i)=0.d0
		yeig(2,i)=0.d0
	enddo
c	
c  solid part
c
	do i=nsolid,js,-1
		yeig(1,i)=ynod(1,i)
		yeig(2,i)=ynod(2,i)
	enddo
c
c  below start radius
c
	do i=js-1,1,-1
	  do n=1,2
	    yeig(n,i)=0.d0
         enddo
	enddo
c
	return
	end
