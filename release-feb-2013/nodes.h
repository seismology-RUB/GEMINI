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
