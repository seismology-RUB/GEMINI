#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#
#  General definitions
#
bindir = .
obsdir = ./objects
moduledir = ./mod
#
CFLAGS = -O3 
ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
	FFLAGS = -O3 -J$(moduledir)
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
x11 = -L/usr/X11R6/lib -lX11
pgplot = -lpgplot /usr/lib/libpng12.so.0 -lz -lX11
#-------------------------------------------------------------
.PHONY:
#--------------------------------------------------------------
STUFF = stuff.o gse20.o sunfortran.o sffFree.o sffSource.o sffInfo.o sffHeader.o sffDatablock.o
BESSEL = bessj0.o bessy0.o bessj1.o bessy1.o bessj.o bessy.o 
#----------------------------------------------------------------
#  - create some h-files and f-files from their m4 templates
#    thus one .h file can be used in different programs with
#    varying values of dimensions, type definitions
#  - to change the default settings in the .m4 files say
#    touch *.h.m4
#    make VARNAME=VALUE target
#
%.f: %.f.m4
	m4 -Dm4_function_type=\`double\ precision\' -DMODES \
	   -Dm4_elcon_type=\`double\ precision\' $< > $@ 
%.h: %.h.m4
	m4 -Dm4_nnd=$(NND) -Dm4_nstatt=$(NSTATT) -Dm4_function_type=\`double\ precision\' \
	   -Dm4_elcon_type=\`double\ precision\' $< > $@
#-------------------------------------------------------------
#  create dependencies on .h-files
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep: *.f90
	incdep-wf.sh $(DEPDIRS) > make.incdep
-include make.incdep
#---------------------------------------------------------------
clean:
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
	-rm -f make.incdep
#----------------------------------------------------------------
#
modeseis: %: %.o commandLine.o mathConstants.o flmodeHeader.o flmodeBody.o sourceReceiver.o frequencyTime.o \
	timeSeries.o splineInterpol.o sortArray.o dlocate.o pythag.o geo2epi.o hbutrw.o sft.o ssft.o \
	$(STUFF) $(BESSEL) hankelFunction.o
	$(F95) -o $(bindir)/gemini/$@ $(obstring)
plotFlmodes: %: %.o commandLine.o mathConstants.o flmodeHeader.o flmodeBody.o splineInterpol.o pgPlotXY.o \
	      zoomBindingsPgPlotSelect.o dlocate.o hankelFunction.o $(BESSEL)
	$(F95) -o $(bindir)/gemini/$@ $(obstring) $(pgplot)
plotFldispersion: %: %.o commandLine.o mathConstants.o flmodeHeader.o flmodeBody.o splineInterpol.o pgPlotXY.o \
	      zoomBindingsPgPlotSelect.o dlocate.o hankelFunction.o $(BESSEL)
	$(F95) -o $(bindir)/gemini/$@ $(obstring) $(pgplot)
