c--------------------------------------------------------
c  $Id: dlocate.f,v 1.1.1.1 2003/01/13 14:27:04 wolle Exp $
c  Locate (double precision)
c
c  Find index j such that xx is between x(j) and x(j+1).
c  If x is ascending:  x(j) <  xx <= x(j+1)
c  If x is descending: x(j) >= xx >  x(j+1)
c
c  returns value in j
c  returns 0 or n if xx is off-scale
c  Sitzt xx auf einer Stuetzstelle, waehlt dlocate immer
c  das Intervall links von der Stuetzstelle !
c-------------------------------------------------------
      subroutine dlocate(xx,n,x,j)
	implicit integer(i-n)
	double precision x(n),xx
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
