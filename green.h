c-----------------------------------------------------------
c  $Id: green.h,v 1.3 2003/10/25 12:16:30 wolle Exp $
c
c  Header file for routines associated with a
c  Green kf-spectrum, dimensions in greendim.h
c----------------------------------------------------------
	integer green_nf1,green_nf2,green_nf,green_comp,green_jump
	integer green_nwn(nff),green_nwnmax
	real green_zs,green_ze,green_df,green_dwn,green_sigma,green_rearth
	character*1 green_sourcetype
c
	common/green/green_zs,green_ze,green_df,green_dwn,green_sigma,
     1	             green_nwn,green_nwnmax,green_nf1,green_nf2,green_nf,
     2	             green_rearth,green_jump,green_comp,green_sourcetype
