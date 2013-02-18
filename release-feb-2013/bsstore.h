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
c--------------------------------------------------------------------
c  Declarations for storage of solution vector
c
c  jstore:    array index after last element of solution vector stored
c  xstore:    x-values
c  ystore:    values of stored solution vector
c
c  requires previous include of nvmax.h 
c---------------------------------------------------------------------
c
	
c---------------------------------------------------------------------
	integer nnstore
	parameter(nnstore=40000)
	integer jstore
	double precision xstore(nnstore)
	double complex ystore(nvmax,nnstore),dystore(nvmax,nnstore)
	common/bsstore/ystore,dystore,xstore,jstore
