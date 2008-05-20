c----------------------------------------------------------------------------------
c  $Id: reciprocity.f,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c                           RECIPROCITY
c
c  Get Green function g(basis-solution,re,r) for fixed component
c  but variable basis solution from g(comp,r,re) for
c  fixed basis solution but variable component.
c
c  nb:            index of basis solution
c  r:             receiver node radius
c  re:            receiver radius
c  gh(1-4):       fixed basis-solution=nb, g(comp=1-4,r,re)
c  green(1-4):    comp(nb)=fixed, g(basis-solution=1-4,re,r)
c                 mit comp(1)=2, comp(2)=1, comp(3)=4, comp(4)=3
c                 see nb4comp in flgevas.f
c---------------------------------------------------------------------------------
	subroutine reciprocity(nb,r,re,elp1,gh,green)
	double complex gh(4),green(4)
	double precision r,re,ratio,elp1
	integer nb
c  
	ratio=(r/re)**2
c
c  here we get the second component of green(re,r) for basis solutions 1-4
c
	if(nb.eq.1) then
		green(1)= ratio*gh(2)
		green(2)=-ratio*gh(1)
		green(3)= ratio*elp1*gh(4)
		green(4)=-ratio*elp1*gh(3)
c
c  here we get the first component of green(re,r) for basis solutions 1-4
c
	else if(nb.eq.2) then
		green(1)=-ratio*gh(2) 
		green(2)= ratio*gh(1)
		green(3)=-ratio*elp1*gh(4)
		green(4)=+ratio*elp1*gh(3)
c
c  here we get the fourth component of green(re,r) for basis solutions 1-4
c
	else if(nb.eq.3) then
		green(1)=+ratio/elp1*gh(2)
		green(2)=-ratio/elp1*gh(1)
		green(3)=+ratio*gh(4)
		green(4)=-ratio*gh(3)
c
c  here we get the third component of gree(re,r) for basis solutions 1-4
c
	else if(nb.eq.4) then
		green(1)=-ratio/elp1*gh(2)
		green(2)=+ratio/elp1*gh(1)
		green(3)=-ratio*gh(4)
		green(4)=+ratio*gh(3)
	endif
c
	return
	end

