####################################################################
### INPUT ##########################################################
####################################################################
SHELL:=/bin/bash

export F_UFMTENDIAN=big

## name of the executable file
EXE = boundary_layer_0.7
## source code path
SRC = pwd

F90=mpiifort

####################################################################
# END INPUT ########################################################
####################################################################

WORK = $(PWD)
PREBUILD = $(WORK)/prebuild
OBJ = $(WORK)/build

FFTWDIR=$(FFTW_DIR)
LAPACKDIR=$(LAPACK_DIR)

#################################################### compiler flags

## general
FFTFLAGS = -I$(FFTWDIR)/include/

#################################################### linker flags

## general 
FFT_LFLAGS = -L$(FFTWDIR)/lib/ -lfftw3_mpi -lfftw3 -lm -lmpi
LAPACK_LFLAGS = -L$(LAPACKDIR) -llapack

LFLAGS += $(FFT_LFLAGS) 
LFLAGS += $(LAPACK_LFLAGS)

#################################################### extra flags

## !!! if_set beg !!!
## INTEL FLAGS
## test    -> -O0 -g -traceback -check all -check bounds
## profile -> -pg -O0 -fno-inline-functions
## run     -> -O3 -ipo
FEXDEB = -O0 -g -traceback -check all -check bounds
FEX = -O3 -ipo



#################################################### objects alpha
OBJECTS = $(OBJ)/mpi.o $(OBJ)/global.o $(OBJ)/Newton_solver.o $(OBJ)/interpolation.o \
	$(OBJ)/equations.o $(OBJ)/rescaled_inlet_bc.o $(OBJ)/boundary_conditions.o \
	$(OBJ)/pressure.o $(OBJ)/input_output.o $(OBJ)/fftz.o $(OBJ)/statistics.o \
	$(OBJ)/initialization.o $(OBJ)/subgrid.o $(OBJ)/wallmodel.o $(OBJ)/finalization.o \
	$(OBJ)/projection.o $(OBJ)/time_integration.o \
	$(OBJ)/monitor.o $(OBJ)/main.o

#################################################### compile 
$(OBJ)/%.o: $(PREBUILD)/%.f90
	@echo compiling $<
	cd $(OBJ);\
		$(F90) -c $(FFTFLAGS)  $(FEX) $< -o $@

.SECONDARY:
$(PREBUILD)/%.f90: $(WORK)/%.f90 $(WORK)/Makefile
	@echo 
	@echo -e "\033[7;35m WORK \033[0m  copying $<"
	@echo to $@
	@cp $< $@
.SECONDARY:
$(PREBUILD)/%.f90: $(SRC)/%.f90 $(WORK)/Makefile 
	@echo 
	@echo -e "\033[7;36m SRC \033[0m copying  $<"
	@echo to $@
	@cp $< $@

############################################################ dir setup

DIRSETUP:
	@if ! [ -d $(OBJ)      ]; then mkdir -p $(OBJ);     echo OBJECTS directory created; fi
	@if ! [ -d $(PREBUILD) ]; then mkdir -p $(PREBUILD);echo PRECOMP directory created; fi

############################################################ build
debug: DIRSETUP $(OBJECTS) 
	@echo " "
	@echo "LINKING... "
	@echo " "
	cd $(OBJ);\
	$(F90) -o $(WORK)/$@ $(OBJECTS) $(LFLAGS) $(FEXDEB)
	@echo " "
	@echo " "	
	echo $@ built, congratulations.


############################################################ build
$(EXE): DIRSETUP $(OBJECTS) 
	@echo " "
	@echo "LINKING... "
	@echo " "
	cd $(OBJ);\
	$(F90) -o $(WORK)/$@ $(OBJECTS) $(LFLAGS) $(FEX)
	@echo " "
	@echo " "	
	echo $@ built, congratulations.

########################################################## message
clean:
	@if [ -d $(OBJ)       ]; then rm -r $(OBJ);      echo OBJECTS directory removed ; fi
	@if [ -d $(PREBUILD)  ]; then rm -r $(PREBUILD); echo PRECOMP directory removed ; fi
	@if [ -f $(EXE)       ]; then rm $(EXE)        ;echo EXECUT. file       removed ; fi
	@if [ -f debug        ]; then rm debug         ;echo EXECUT. file       removed ; fi

all: $(EXE)
