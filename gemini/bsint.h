c-----------------------------------------------------------------
c  $Id: bsint.h,v 1.2 2003/04/03 12:09:38 wolle Exp $
c
c  Class declarations for bsint
c-----------------------------------------------------------------
	integer maxstep,imax,nuse,maxsub,imin
	double precision tiny,hmin,shrink,grow,eps
	parameter(maxstep=20000,tiny=1.d-10,hmin=1.d-12,eps=1.d-6)
	parameter(shrink=0.95,grow=1.05,imax=7,nuse=5,maxsub=24,imin=2)
