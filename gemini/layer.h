c--------------------------------------------------------------
c  $Id: layer.h,v 1.1.1.1 2003/01/13 14:27:03 wolle Exp $
c
c  Declarations for layer description of an earth model
c
c  nlayer:			max number of layers
c  nlay:			number of layers
c  noc:			layer index of outer core
c  nic:			layer index of inner core
c  rb:			radii of layer boundaries
c  iflso:			fluid-solid flag (1=fluid, 0=solid)
c  layer_nlactive:		 active layer index
c  nrna(nl):   node index of bottom node in layer
c  nrne(nl):   node index of top node in layer
c--------------------------------------------------------------
	include 'laydim.h'
	double precision layer_bali
	parameter(layer_bali=0.1)
	double precision layer_stepsize_wlfrac
	double precision layer_rb(0:nlayer)
	integer layer_iktop(0:nlayer),layer_ikbot(nlayer),layer_nlay,layer_iflso(nlayer)
	integer layer_nic,layer_noc,layer_nlactive
	integer layer_nrna(nlayer),layer_nrne(nlayer)
	common/layer/layer_rb,layer_stepsize_wlfrac,
     1	             layer_iktop,layer_ikbot,layer_nlay,layer_iflso,layer_nic,
     1	             layer_noc,layer_nlactive,layer_nrna,layer_nrne
