#--------------------------------------------------------------------------
#	Copyright 2013 Wolfgang Friederich
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
#---------------------------------------------------------------
#  Makefile for benchmarking
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/tcsh
#-----------------------------------------------------------
#           BENCHMARK
workdir = ..
#-----------------------------------------------------------------
#            DIVERSE DIRECTORIES
#
models       = $(workdir)/models
info         = $(workdir)/info
sff          = $(workdir)/sff
temp         = $(workdir)/temp
ref          = $(workdir)/ref
#--------------------------------------------------------------------
#  Benchmark
#-------------------------------
dirs:
	if (! -e $(temp)) mkdir -p $(temp)
	if (! -e $(sff)) mkdir -p $(sff)
gfseam:
	gfdsvrkf -fmin .006 -fmax 1000. -tlen .128 -klim ".5 .5" -ze 0. -xlen .5 -s 0 -fnod 1000. -sig 5. -wmrg .1 \
		-msd "9. 1" -on nodes.seam.out -o $(temp)/seam $(models)/seam1.nm
gfshallow:
	gfdsvrkf -fmin .006 -fmax 100. -tlen 1.28 -klim "10. 10." -ze 0. -xlen 1. -s 0 -fnod 100. -sig 5. -wmrg .1 \
		-msd "10. 1" -on nodes.shallow.out -o $(temp)/shallow $(models)/flach.nm
gfregio:
	gfdsvrkf -fmin .006 -fmax .2 -tlen 512. -klim ".8 .8" -ze 0. -xlen 5000. -s 1 -fnod .2 -sig 5. -wmrg .1 \
		-msd "30000. 1" -on nodes.regio.out -o $(temp)/regio $(models)/ak135q.nm
gfseamseis:
	gfdsvrkfseis -lp "5 500." -hp "5 40." -e 40000 -gvw "2. .1" -dt .00025 -zs 8.5 \
		-o $(sff)/seam.gfdsvrkf.sff $(info)/seam.info $(temp)/seam.sf
gfshallowseis:
	gfdsvrkfseis -lp "5 50." -hp "2 10." -e 40000 -gvw "2. 1." -dt .0025 -zs 0.0 \
		-o $(sff)/shallow.gfdsvrkf.sff $(info)/shallow.info $(temp)/shallow.sf
gfregioseis:
	gfdsvrkfseis -lp "5 .1" -hp "2 .01" -e 40000 -gvw "2. 50." -dt 1. -zs 26999. \
		-o $(sff)/regio.gfdsvrkf.sff $(info)/regio.info $(temp)/regio.sm
diffseam:
	diffDataSynthetics -o misfit.seam $(sff)/seam.gfdsvrkf.sff $(ref)/seam.gfdsvrkf.sff
diffshallow:
	diffDataSynthetics -o misfit.shallow $(sff)/shallow.gfdsvrkf.sff $(ref)/shallow.gfdsvrkf.sff
diffregio:
	diffDataSynthetics -o misfit.regio $(sff)/regio.gfdsvrkf.sff $(ref)/regio.gfdsvrkf.sff
