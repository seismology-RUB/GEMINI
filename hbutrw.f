c---------------------------------------------------------
c         H B U T R W
c
c  computes spectral response of a Butterworth filter with
c  real frequencies
c
c  ckind: either 'low ' or 'high'
c  nord:  order of filter
c  nf:    number of frequencies
c  wc:    corner frequency
c  wa:    first frequency
c  dw:    frequency spacing
c  h:     complex transfer function
c---------------------------------------------------------
c
      subroutine hbutrw(ckind,nord,nf,wc,wa,dw,h)
      complex h(nf),ci,ep(20)
	real wc,dw,pi,wa,w
	integer i,j,nord,nf
      character ckind*4
c
      pi=4.*atan(1.)
      ci=(0.,1.)
c
      do 20 j=1,nord
 20   ep(j)=cexp(ci*pi*(float(j)-.5)/float(nord))
c
c--------------------------------    Low pass
      if(ckind.eq.'low ') then
       do 10 i=1,nf
       w=(wa+(i-1)*dw)/wc
       h(i)=1.
       do 10 j=1,nord
 10    h(i)=-h(i)*ci/(w-ep(j))
c-------------------------------------  High pass
      else if(ckind.eq.'high') then
       do 30 i=1,nf
       w=(wa+(i-1)*dw)/wc
       h(i)=1.
       do 30 j=1,nord
 30    h(i)=h(i)*w/(w-ep(j))
      endif
c
      return
      end
