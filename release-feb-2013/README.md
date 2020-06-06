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
Program package GEMINI:
=====================================

A program package to calculate Green functions
in a depth-dependent 1D spherical or plane earth model. Can be used from
continental scale calculations down to ultrasonic waves if epicentral
distances, penetration depth and max frequency scale appropriately.
Tested frequency range from 0.01 Hz to 500000 Hz.

Requirements:
----------------

- gfortran (at least 4.2) and gcc for compilation
- GNU make
- m4 macro processor

Compilation:
-------------
Set the variable F95 in your shell environment to the path to your Fortran compiler.
Use Makefile and Make.f90 for compiling the Green function codes:

* touch *.m4
* make gfdsvrkf
* make gfdsvrkf_mpi
* make -f Make.f90 gfdsvrkfseis

The m4 macro processor is used to automatically generate code needed for
gfdsvrkf from a common template. If anything goes wrong do a "make clean"
and repeat the compilation procedure.

Usage:
--------
The code is self-explaining. Just enter the name of the executable and you get a
description of arguments and options. For theoretical background look into the
documentation or into the paper by Friederich and Dalkolmo (GJI, 1995).

Info file:  all information about source and receiver is provided by an info file (see examples folder).

gfdsvrkf:      calculate frequency-wavenumber spectra for the displacement stress vector (DSV)
               with components U, R, V, S, W, T and optionally spatial derivatives either for
               for one source and many receiver depths or for one receiver and many source depths.

gfdsvrkf_mpi:  an embarassingly parallel version of gfdsvrkf which distributes calculations for
               different frequencies to available processors.

gfdsvrkfseis:  calculate synthetic seismograms (takes stations and components from an info file),
               current ourput formats: SFF and SSA (a special file with allows direct access).

Look into the examples folder where you find Makefiles for different applications of gfdsvrkf.

Seismogram format SFF:
----------------------
SFF (Stuttgart File Format) is a slight modification of GSE2.0
with CMPR6 compression. In the folder libstuff you find some Fortran modules which provide
routines for reading and wirting SFF files. As far as I know SFF files can be read in using Seismic Handler.
Look into stuff.f for a detailed description of the format.

Format SSA:
------------
SSA is very helpful when a big amount of seismograms are written because it permits direct access to any seismogram
upon reading. There is a Fortran module "include/ssaDataset.f90" which provides routines for reading SSA-files.
