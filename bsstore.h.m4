c--------------------------------------------------------------------
c  $Id: bsstore.h,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c  Declarations for storage of solution vector
c
c  jstore:    array index after last element of solution vector stored
c  xstore:    x-values
c  ystore:    values of stored solution vector
c
c  requires previous include of nvmax.h 
c---------------------------------------------------------------------
c
	define(m4_function_type,`double complex')
c---------------------------------------------------------------------
	integer nnstore
	parameter(nnstore=40000)
	integer jstore
	double precision xstore(nnstore)
	m4_function_type ystore(nvmax,nnstore),dystore(nvmax,nnstore)
	common/bsstore/ystore,dystore,xstore,jstore
