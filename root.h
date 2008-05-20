c------------------------------------------------------------
c  $Id: root.h,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c   class for root finding procedures
c------------------------------------------------------------
	integer nvall
	parameter(nvall=512)
	integer kbr(nvall),nval
	double precision m5(nvall),elbr(nvall)
	common/rootfind/elbr,m5,kbr,nval
