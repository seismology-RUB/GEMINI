#!/bin/sh
#--------------------------------------------------------------------------
#	Copyright 2013 Thomas Forbriger, Wolfgang Friederich
#
#	This file is part of Gemini II.
#
#	Gemini II is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 2 of the License, or
#	any later version.
#
#	Gemini II is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
# ============================================================================
#
for dir in . $@
do
for d in $dir/*.f
do
  fina=`basename $d .f`.o
  depna=`cat $d | grep '^	include' | cut -f 2 -d \' | sort | uniq`
  echo $fina: $depna
done
for d in $dir/*.f.m4
do
  fina=`basename $d .f.m4`.o
  depna=`cat $d | grep '^	include' | cut -f 2 -d \' | sort | uniq`
  echo $fina: $depna
done
for d in $dir/*.f90
do
  fina=`basename $d .f90`.o
  depna=`cat $d | grep '^	use' | cut -f 2 -d ' ' | cut -f 1 -d ',' | sort | uniq`
  for name in $depna
  do
	  echo $fina: $name.o
  done
done
done
# ----- END OF incdep ----- 
