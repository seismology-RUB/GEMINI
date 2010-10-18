c---------------------------------------------------------------------
c                           SPHDERIVS
c
c  compute derivatives w.r.t. r from DSV and overwrite DSV with U,UP,V,VP	
c  treats also fluid case correctly
c---------------------------------------------------------------------------
	subroutine sphderivs(gf)
	include 'nodesdim.h'
	include 'nodes.h'
	include 'layer.h'
	include 'zero.h'
	include 'nvmax.h'
	double complex gf(nnd,6),yh(4)
	double complex zom,a(nvmax,nvmax),b(nvmax)
	double precision om,r,ro,rop
	integer i,j,top
	integer nr,nl
	common/omega/zom,om
c
c  loop over nodes
c
	do nr=jstu,nnod
		r=rnod(nr)
		call layer_getindex(r,nl,top)
		if(top.eq.1.and.dabs(rnod(nr-1)-rnod(nr)).lt.zero) nl=nl+1
		call layer_setnlactive(nl)
		if(layer_iflso(nl).eq.1) then
			call sysliq(r,a,b)
			call earthmodel_rholiq(r,ro,rop)
			yh(1)=0.d0
			do j=1,2
				yh(1)=yh(1)+a(1,j)*gf(nr,j)
			enddo
c
c  V = -R/(ro*om^2*r)
c  dV/dr = -1/(ro*om^2*r) R' + R/(ro^2*om^4*r^2)*(-(-1)*(rop*om^2*r+ro*om^2)
c        = +1/r*U + R*(rop*r+ro)/(ro^2*om^2*r^2)
c  R' = -om^2*ro*U
c
			gf(nr,3) = -gf(nr,2)/(ro*zom*zom*r)
			gf(nr,4) = 0.d0
			gf(nr,5) = yh(1)
			gf(nr,6) = gf(nr,1)/r+gf(nr,2)*(rop*r+ro)/(ro*zom*r)**2
		else
			call sysspher(r,a,b)
			do i=1,3,2
				yh(i)=0.d0
				do j=1,4
					yh(i)=yh(i)+a(i,j)*gf(nr,j)
				enddo
			enddo
			gf(nr,5) = yh(1)
			gf(nr,6) = yh(3)
		endif
	enddo
c
	return
	end
