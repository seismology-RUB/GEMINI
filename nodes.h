c---------------------------------------------------------------
c  Declarations for nodes-class
c
c  nnd:		max number of nodes
c  rnod:	radii of nodes in km
c  nnod:	number of nodes
c  jsd,jsu:	indices of rs- and rs+
c  jed,jeu:	indices of re- and re+
c  jwd,jwu: 	indices of rw- and rw+, where rw is the water-solid
c		interface if present
c  jstu:	index of first node above starting radius
c
c  requires previous include of nodesdim.h
c---------------------------------------------------------------
	double precision rnod(nnd)
	integer nnod,jsd,jsu,jed,jeu,jwd,jwu,jstu
	common/rnodes/rnod,nnod,jsd,jsu,jed,jeu,jwd,jwu,jstu
