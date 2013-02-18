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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/tcsh
#----------------------------------------------------------------
#  General definitions
#
bindir = ./bin
obsdir = ./objects
moduledir = ./mod
#
CFLAGS = -O3 
ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
	FFLAGS = -O3 -J$(moduledir) -Wunused-variable -Wuninitialized
endif
#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = 
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f files.
#  Output of f2c stays in directory where .f files are.
#  Substitute .f by .c in dependency and compile target.
#  Then remove .c files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $(obsdir)/$@
%.o: %.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
%.o: %.f
	$(F95) -c $(FFLAGS) -fimplicit-none -ffixed-line-length-132 $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Library paths
#
#-------------------------------------------------------------
.PHONY:
#--------------------------------------------------------------
STUFF = stuff.o gse20.o sunfortran.o sffFree.o sffSource.o sffInfo.o sffHeader.o sffDatablock.o sffDataSection.o
BESSEL_NETLIB = bessj0y0.o bessj1y1.o 
TIME_SERIES = timeSeries.o fourierTransform.o filterCoefficientsDecimate.o recursiveFilterCoefficients.o pythag.o
SA = streamAccess.o flexibleType.o simpleString.o kindDefinitions.o
SSA = ssaDataset.o streamAccess.o flexibleType.o simpleString.o
#----------------------------------------------------------------
#  create dependencies on .h-files
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make90.incdep:
	./makeDepFromUseInclude.py $(DEPDIRS) > make90.incdep
	if (! -e $(bindir)) mkdir -p $(bindir)
	if (! -e $(moduledir)) mkdir -p $(moduledir)
	if (! -e $(obsdir)) mkdir -p $(obsdir)
-include make90.incdep
#---------------------------------------------------------------
clean:
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
	-rm -f make90.incdep
#----------------------------------------------------------------
#
gfdsvrkfseis: %: %.o commandLine.o mathConstants.o sourceReceiver.o frequencyTime.o sortArray.o \
	pythag.o geo2epi.o hbutrw.o hbutcw.o locatePoint.o greenSection.o fourierTransform.o \
	primitiveTypeEncoding.o greenFrequencyWavenumber.o greenDisplacementSpectrum.o wavenumberIntegralMapping.o \
	dayOfYear.o instrumentResponse.o errorMessage.o realloc.o ssaDataset.o fileUnitHandler.o timeUtils.o \
	geminiEarthModel.o dataSegy.o binarySegy.o traceSegy.o dateTime.o seismicEvent.o seismicStation.o \
	$(SA) $(STUFF) $(BESSEL_NETLIB) $(TIME_SERIES)
	$(F95) -o $(bindir)/$@ $(obstring)
diffDataSynthetics: %: %.o commandLine.o dateTime.o timeUtils.o realloc.o $(STUFF) $(TIME_SERIES)
	$(F95) -o $(bindir)/$@ $(obstring)
