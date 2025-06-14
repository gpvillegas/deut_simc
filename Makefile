## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## ARGONNE DEFAULT SETUP FLAGS:
#simcdir = .
#Csoft = /disk1/users/reinhold/Csoft/05Dec1996

## CEBAF DEFAULT SETUP FLAGS:
simcdir = .
#Csoft = /group/hallc/Csoft/Analyzer
#Csoft = /group/c-gep/jones/Linux_fc8/

## U Regina SETUP FLAGS 
#simcdir = .
#Csoft = /home/huberg/r2d2/simc/

## THE REST SHOULD BE OK WITHOUT MODIFICATION.

## This tells make not to delete these target files on error/interrupt (see man page)
.PRECIOUS: *.o sos/*.o hms/*.o hrsl/*.o hrsr/*.o shms/*.o calo/*.o mini_cern/*.o

RM        = rm -f 
SHELL     = /bin/sh
S	= $(simcdir)/sos/
H	= $(simcdir)/hms/
L	= $(simcdir)/hrsl/
R	= $(simcdir)/hrsr/
A	= $(simcdir)/shared/
SH	= $(simcdir)/shms/
T       = $(simcdir)/cteq5/
C       = $(simcdir)/calo/
CM      = $(simcdir)/cern_mini/
EL      = $(simcdir)/eloss/

#OBJ1	= target.o brem.o gauss1.o NtupleInit.o NtupleClose.o enerloss_new.o

OBJ1	= target.o brem.o gauss1.o NtupleInit.o NtupleClose.o
OBJ2	= radc.o init.o dbase.o physics_kaon.o physics_pion.o physics_delta.o physics_proton.o loren.o sf_lookup.o response_functions_mod.o
OBJ3    = semi_physics.o rho_physics.o rho_decay.o generate_rho.o trg_track.o semi_dilution.o LagetXsec1.o MSxsec.o MSxsec_f.o
OBJ4	= results_write.o event.o mt19937.o jacobians.o
OBJ5	= $(A)musc.o $(A)musc_ext.o $(A)project.o $(A)transp.o
OBJ6	= $(A)rotate_haxis.o $(A)rotate_vaxis.o $(A)locforunt.o
OBJ7	= $(H)mc_hms.o $(H)mc_hms_hut.o $(H)mc_hms_recon.o
OBJ8	= $(S)mc_sos.o $(S)mc_sos_hut.o $(S)mc_sos_recon.o
OBJ9	= $(R)mc_hrsr.o $(R)mc_hrsr_hut.o $(R)mc_hrsr_recon.o
OBJA	= $(L)mc_hrsl.o $(L)mc_hrsl_hut.o $(L)mc_hrsl_recon.o
OBJB	= $(SH)mc_shms.o $(SH)mc_shms_hut.o $(SH)mc_shms_recon.o
OBJC    = $(T)Ctq5Pdf.o
OBJD    = $(C)mc_calo.o $(C)mc_calo_recon.o
OBJE    = $(CM)lfit.o $(CM)lfitw.o $(CM)fint.o 
OBJEL   = $(EL)enerloss_new.o $(EL)eloss_mean.o $(EL)straggling.o $(EL)glandg.o $(EL)gvaviv.o $(EL)glande.o $(EL)glands.o $(EL)ranecu.o
my_objs	=  $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJA) $(OBJB) $(OBJC) $(OBJD) $(OBJE) $(OBJEL)

MYOS := $(subst -,,$(shell uname))
CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

#For use with gfortran compiler
# -fno-automatic - all program storage treated as static
ifeq ($(MYOS),Linux)
#  LIBROOT = CTP/O.Linux/Linux/lib
# JLab
#  CERN_ROOT = /apps/cernlib/i386_fc8/2005
# FIU
  CERN_ROOT = /cern/pro/
# 32 bit, standard Fedora distributuion
#  CERN_ROOT = /usr/lib/cernlib/2006
# 64 bit, standard Fedora distributuion
#  CERN_ROOT =  /usr/lib64/cernlib/2006 
# for debugging
# FFLAGSA=-g -O0 -W -ffixed-line-length-132 -fno-automatic -fdefault-real-8 -fdefault-double-8  -finit-local-zero -fbounds-check
#
# type mismatch troubleshooting: added flags -fallow-argument-mismatch -std=legacy --gpv
  FFLAGSA=-O -W -ffixed-line-length-132 -fno-automatic -fdefault-real-8 -fdefault-double-8  -finit-local-zero -fallow-argument-mismatch -std=legacy -fno-align-commons
#
  #FFLAGSA=-O -W -ffixed-line-length-132 -fno-automatic -fdefault-real-8 -fdefault-double-8 -finit-local-zero -fno-align-commons -ffixed-line-length-none -g -fbounds-check
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
#  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
  FC  := gfortran
  F77 := gfortran
endif

# Mac OSX, Leopard (checked on version 10.5.6)
# Tested using patched gfortran compiler from www-jlc.kek.jp/~fujiik/macosx/10.5.X/HEPonX
# I had trouble compiling geant stuff using the deafult gfortran that shipped with
# my version of Leopard - YMMV (your mileage may vary)
# Only change needed was to add  the -fno-range-check flag: there was some problem
# with the default integer size, specifically in mt19937.f.
# Note that the CTP libraries still end up in the O.Linux directory...
ifeq ($(MYOS),Darwin)
#  LIBROOT = CTP/O.Linux/Linux/lib
#  CERN_ROOT = /apps/cernlib/i386_fc8/2005
  # FFLAGSA=-O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  # the addition was necessary due to brem.f, rho_physics.f
# for debugging
#FFLAGSA=-g -O0 -W -ffixed-line-length-132 -fno-automatic -fdefault-real-8 -fdefault-double-8 -finit-local-zero -fbacktrace -fbounds-check
  FFLAGSA= -O -W -ffixed-line-length-132 -fno-automatic -fdefault-real-8 -fdefault-double-8 -finit-local-zero -fbacktrace 
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
#  OTHERLIBS = -L$(LIBROOT) -lctp \
#        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
  FC  := gfortran
  F77 := gfortran
endif

%.o: %.f90
	$(F77) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

$(S)/%.o: $(S)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(H)/%.o: $(H)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(L)/%.o: $(L)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(R)/%.o: $(R)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(SH)/%.o: $(SH)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(CM)/%.o: $(CM)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(EL)/%.o: $(EL)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

none: simc

all: simc, test_parameters, test_nml_default

simc: simc.o $(my_objs) Makefile 
	$(F77) $(OSF_SHARED) -o $@ $(FFLAGS) $(my_objs) simc.o $(OTHERLIBS)

test_parameters: test_parameters.o dbase_namelists.inc 
	$(F77) $(OSF_SHARED) -o $@ $(FFLAGS) test_parameters.o $(OTHERLIBS)

test_nml_default: test_nml_default.o dbase_namelists.inc  
	$(F77) $(OSF_SHARED) -o $@ $(FFLAGS) test_nml_default.o $(OTHERLIBS)

#CTP/O.Linux/Linux/lib/libctp.a: 
#	make -C CTP
# These routines have HP problems, and need to be compiled without optimization.

simc.o: simc.f
	$(F77) $(FFLAG1) simc.f

#init.o: init.f
#	$(F77) $(FFLAG1) init.f
#
#event.o: event.f
#	$(F77) $(FFLAG1) event.f
#
#dbase.o: dbase.f
#	$(F77) $(FFLAG1) dbase.f
#
#physics.o: physics.f
#	$(F77) $(FFLAG1) physics.f
#
#results_write.o: results_write.f
#	$(F77) $(FFLAG1) results_write.f


clean:
	$(RM) *.o $(H)*.o $(S)*.o $(L)*.o $(R)*.o $(SH)*.o $(A)*.o $(T)*.o $(C)*.o $(CM)*.o $(EL)*.o simc

