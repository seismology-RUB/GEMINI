c--------------------------------------------------------
c  Locate
c
c  locate index of point xx with x(i) < xx <= x(i+1)
c  returns value in j
c  returns 0 or n if xx is off-scale
c-------------------------------------------------------
      subroutine Locate(xx,n,x,j)
	implicit integer(i-n)
	real x(n),xx
c
c  search interval
c
      jl=0
	ju=n+1
 10	if(ju-jl.gt.1) then
		jm=(ju+jl)/2
		if((x(n).gt.x(1)).eqv.(xx.gt.x(jm))) then
			jl=jm
		else
			ju=jm
		endif
		goto 10
	endif
	j=jl
c
      return
      end
