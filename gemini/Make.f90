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
bindir = ../bin
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
vpath %.f90 ../include ../plotlib ../libstuff
vpath %.c ../libstuff
vpath %.f ../libstuff ../wolib ../include
#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = ../plotlib ../libstuff ../include
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
x11 = -L/usr/X11R6/lib -lX11
pgplot = -lpgplot
#-------------------------------------------------------------
.PHONY:
#--------------------------------------------------------------
STUFF = stuff.o gse20.o sunfortran.o sffFree.o sffSource.o sffInfo.o sffHeader.o sffDatablock.o sffDataSection.o
BESSEL = bessj0.o bessy0.o bessj1.o bessy1.o bessj.o bessy.o 
TIME_SERIES = timeSeries.o fourierTransform.o filterCoefficientsDecimate.o recursiveFilterCoefficients.o pythag.o dfftpack.o
SA = streamAccess.o flexibleType.o simpleString.o kindDefinitions.o
SSA = ssaDataset.o streamAccess.o flexibleType.o simpleString.o
FLMODE = flmodeHeader.o flmodeBody.o flmodeContainer.o flmodeOvertones.o splineInterpol.o hankelFunction.o locatePoint.o linearInterpol.o
KERNEL = flmodeKernelBody.o flmodeKernelContainer.o flmodeKernelOvertones.o locatePoint.o nodeEarthmodel.o
GREEN = 
#----------------------------------------------------------------
#  create dependencies on .h-files
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	../bin/makeDepFromUseInclude.py $(DEPDIRS) > make.incdep
	if (! -e mod) mkdir -p mod
	if (! -e objects) mkdir -p objects
-include make.incdep
#---------------------------------------------------------------
clean:
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
	-rm -f make.incdep
#----------------------------------------------------------------
#
gfdsvrkfseis: %: %.o commandLine.o mathConstants.o sourceReceiver.o frequencyTime.o sortArray.o \
	pythag.o geo2epi.o hbutrw.o hbutcw.o locatePoint.o greenSection.o fourierTransform.o dfftpack.o \
	primitiveTypeEncoding.o greenFrequencyWavenumber.o greenDisplacementSpectrum.o wavenumberIntegralMapping.o \
	dayOfYear.o instrumentResponse.o errorMessage.o realloc.o ssaDataset.o $(SA) $(STUFF) $(BESSEL) $(TIME_SERIES)
	$(F95) -o $(bindir)/$@ $(obstring)
plotGreenFKSpectra: %: %.o commandLine.o mathConstants.o pgPlotWindow.o pgPlotImage.o zoomBindingsPgPlotSelect.o \
	greenFrequencyWavenumber.o wavenumberIntegralMapping.o locatePoint.o realloc.o errorMessage.o plotUtilities.o \
	fourierTransform.o dfftpack.o pgColorScale.o $(SA) $(BESSEL)
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot)
#
