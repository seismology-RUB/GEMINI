c---------------------------------------------------------------
c  Declarations for nodes-class
c
c  nnd:		max number of nodes
c  ynod:		solution vector stored at nodes
c
c  require previous include of nodes.h
c---------------------------------------------------------------
c  Type specification before compilation using m4
c
	define(m4_function_type,`double complex')
c---------------------------------------------------------------
	m4_function_type ynod(6,nnd)
	common/nodesstore/ynod

