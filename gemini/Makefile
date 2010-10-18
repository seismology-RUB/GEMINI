#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/tcsh
#---------------------------------------------------------------------
#  General definitions
#
F2C = f2c
AR = ar
F2CFLAGS = -f -u
CFLAGS = -O3
FFLAGS = -O3 -J$(moduledir) -Wunused-variable -Wuninitialized
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f files.
#  Output of f2c stays in directory where .f files are.
#  Substitute .f by .c in dependency and compile target.
#  Then remove .c files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
%.o: %.f
	$(F95) -c $(FFLAGS) -fimplicit-none -ffixed-line-length-132 $< -o $(obsdir)/$@
#%.o: %.f
#	$(F2C) $(F2CFLAGS) -d$(dir $<) $<
#	$(CC) -c $(CFLAGS) $(dir $<)$*.c -o $(obsdir)/$@
#	-rm -f $(dir $<)$*.c
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Paths
#
bindir = ../bin
obsdir = ./objects
moduledir = ./mod
#
#  Library paths
#
f2clib = -L/usr/lib -lf2c -lm
misc = -L$(WHOME)/lib -lwo$(ARCH)
stuff = -L$(WHOME)/lib -lf2cstuff$(ARCH)
mpi = -lmpi
#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
vpath %.f ../wolib ../libstuff
vpath %.c ../libstuff
vpath %.f90 ../include
#---------------------------------------------------------
#  additional folders to be searched for dependencies
#
DEPDIRS = ../include
#-----------------------------------------------------------
.PHONY: all
all: flgevask borndsvk stfinvtd seismogram bornspecforce bornspecgreen
#
# Dependent o-files
#
OBS_FLGEVAS = propag.o sysminors.o stavani.o sysliq.o earthmodel_procs.o \
	syssolidt.o nodes_procs.o startr.o flnm_procs.o bsint_procs.o zetl.o \
	sterms.o minorint.o greenint.o reciprocity.o zelcon.o \
	frechet.o layer_procs.o flgevas_init.o isoptbelpar.o sphancf.o \
	greenkfsph.o greenkftor.o systor.o torint.o frechettor.o hatkernel.o \
	stavsurf.o sphan2cf.o
OBS_GFDSV = propag.o sysminors.o stavani.o sysliq.o earthmodel_procs.o \
	syssolidt.o nodes_procs.o startr.o flnm_procs.o bsint_procs.o zetl.o \
	sterms.o minorint.o greenint.o reciprocity.o zelcon.o \
	layer_procs.o systor.o torint.o stavsurf.o sphan2cf.o sphancf.o
#-------------------------------------------------------------
#  create dependencies on .h-files
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	incdep-wf.sh $(DEPDIRS) > make.incdep
	if (! -e mod) mkdir -p mod
	if (! -e objects) mkdir -p objects
-include make.incdep
#---------------------------------------------------------------
clean:
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
	-rm -f make.incdep
#---------------------------------------------------------------
STUFF = stuff.o gse20.o sunfortran.o
BESSEL = bessj0.o bessy0.o bessj1.o bessy1.o bessj.o bessy.o
#----------------------------------------------------------------
#  - create some h-files and f-files from their m4 templates
#    thus one .h file can be used in different programs with
#    varying values of the dimensions
#  - to change the default settings in the .m4 files say
#    touch *.h.m4
#    make VARNAME=VALUE target
#
NSTATT=200
NKFMAX=150000
NND=500
NPP=1000
%.h: %.h.m4
	m4 -Dm4_nstatt=$(NSTATT) -Dm4_nnd=$(NND) -Dm4_nkfmax=$(NKFMAX) \
	   -Dm4_npp=$(NPP) \
	   -Dm4_elcon_type=\`double\ complex\' $< > $@
%.f: %.f.m4
	m4 -Dm4_function_type=\`double\ complex\' -DGREEN \
	   -Dm4_elcon_type=\`double\ complex\' $< > $@ 
#
gfdsvrkf.o: earthmodel.h flnm.h nodesdim.h nodes.h
gfdsvrkf_mpi.o: earthmodel.h flnm.h nodesdim.h nodes.h
#----------------------------------------------------------------
#       Targets
#
gfdsvrkf.o: gfdsvrkf.f90
	$(F95) -c $(FFLAGS) -ffixed-form -ffixed-line-length-132 $< -o $(obsdir)/$@
gfdsvrkf: %: %.o $(OBS_GFDSV) jump_solution_sph.o jump_solution_tor.o commandLine.o mathConstants.o \
	trim.o dlocate.o dsplco.o zsplco.o streamAccess.o primitiveTypeEncoding.o flexibleType.o simpleString.o \
	sphderivs.o torderivs.o sysspher.o systor.o sysliq.o
	$(F95) -o $(bindir)/$@ $(obstring)
gfdsvrkf_mpi.o: gfdsvrkf_mpi.f90
	mpif90 -c $(FFLAGS) -ffixed-form -ffixed-line-length-132 $< -o $(obsdir)/$@
mpiSupport.o: mpiSupport.f90
	mpif90 -c -fintrinsic-modules-path /usr/lib/openmpi/lib -J$(moduledir) $< -o $(obsdir)/$@
gfdsvrkf_mpi: %: %.o $(OBS_GFDSV) jump_solution_sph.o jump_solution_tor.o commandLine.o mathConstants.o mpiSupport.o \
	trim.o dlocate.o dsplco.o zsplco.o streamAccess.o primitiveTypeEncoding.o flexibleType.o simpleString.o \
	sphderivs.o torderivs.o  sysspher.o systor.o sysliq.o
	mpif90 -o $(bindir)/$@ $(obstring)
