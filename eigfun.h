c--------------------------------------------------------------
c  Class for eigenfunction values at neighbouring nodes
c  needed for interpolation and integration over
c  eigenfunctions
c-------------------------------------------------------------
	include 'eigdim.h'
	integer ef_node
	character ef_typ*1
	double precision ef_elp1
	double precision ef_y1(nvx),ef_y2(nvx)
	common/eigfun/ef_elp1,ef_y1,ef_y2,ef_node,ef_typ
