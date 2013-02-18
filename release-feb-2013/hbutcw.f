c--------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
c---------------------------------------------------------
c         H B U T C W
c
c  computes spectral response of a Butterworth filter with
c  complex frequencies
c
c  ckind: either 'low ' or 'high'
c  nord:  order of filter
c  nf:    number of frequencies
c  wc:    real corner frequency
c  wa:    complex first frequency
c  dw:    real frequency spacing
c  h:     complex transfer function
c---------------------------------------------------------
c
      subroutine hbutcw(ckind,nord,nf,wc,wa,dw,h)
      complex h(nf),ci,ep(20),wa,w
	real wc,dw,pi
	integer i,j,nord,nf
      character ckind*4
c
      pi=3.141592
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
