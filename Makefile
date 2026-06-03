####################################################################
### GPU Build with NVIDIA HPC SDK + OpenACC                      ###
####################################################################

SHELL:=/bin/bash

export F_UFMTENDIAN=big

## name of the executable file
EXE = boundary_layer_gpu
## source code path
SRC = pwd

## NVIDIA HPC SDK paths (set NVHPC_ROOT or use default)
NVHPC_ROOT ?= /opt/nvidia/hpc_sdk/Linux_x86_64/24.3
F90 = $(NVHPC_ROOT)/comm_libs/mpi/bin/mpifort

####################################################################
# END INPUT ########################################################
####################################################################

WORK = $(PWD)
PREBUILD = $(WORK)/prebuild
OBJ = $(WORK)/build

FFTWDIR ?= $(FFTW_DIR)

#################################################### compiler flags

## OpenACC + GPU flags
ACCFLAGS = -acc -gpu=cc80 -Minfo=accel

## FFTW (still needed for CPU fallback of FFTW-MPI pressure solver)
FFTFLAGS = -I$(FFTWDIR)/include/

#################################################### linker flags

## FFTW (CPU, for pressure solver until cuFFT replacement is complete)
FFT_LFLAGS = -L$(FFTWDIR)/lib/ -lfftw3_mpi -lfftw3 -lm

## cuFFT
CUDA_LFLAGS = -cudalib=cufft

## LAPACK (from NVHPC)
LAPACK_LFLAGS = -llapack

LFLAGS += $(FFT_LFLAGS)
LFLAGS += $(LAPACK_LFLAGS)
LFLAGS += $(CUDA_LFLAGS)

#################################################### extra flags

FEXDEB = -O0 -g -Mbounds $(ACCFLAGS)
FEXOPT = -O2 $(ACCFLAGS)

debug:  FEX = $(FEXDEB)
$(EXE): FEX = $(FEXOPT)

#################################################### objects alpha
OBJECTS = $(OBJ)/mpi.o $(OBJ)/global.o $(OBJ)/cufft_solver.o $(OBJ)/Newton_solver.o $(OBJ)/interpolation.o \
	$(OBJ)/equations.o $(OBJ)/rescaled_inlet_bc.o $(OBJ)/boundary_conditions.o \
	$(OBJ)/pressure.o $(OBJ)/input_output.o $(OBJ)/fftz.o $(OBJ)/statistics.o \
	$(OBJ)/initialization.o $(OBJ)/subgrid.o $(OBJ)/wallmodel.o $(OBJ)/finalization.o \
	$(OBJ)/projection.o $(OBJ)/time_integration.o \
	$(OBJ)/params.o $(OBJ)/miscel.o $(OBJ)/monitor.o $(OBJ)/main.o

#################################################### compile
$(OBJ)/%.o: $(PREBUILD)/%.f90
	@echo compiling $<
	cd $(OBJ);\
		$(F90) -c $(FFTFLAGS) $(FEX) $< -o $@

# Files that trigger nvfortran 24.3 ICE at -O2: compile at -O0
$(OBJ)/statistics.o: $(PREBUILD)/statistics.f90
	@echo "compiling $< (O0 workaround for nvfortran ICE)"
	cd $(OBJ);\
		$(F90) -c $(FFTFLAGS) -O0 -acc -gpu=cc80 $< -o $@

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
.PHONY: clean cleandebug cleanobj
clean:
	@if [ -d $(OBJ)       ]; then rm -r $(OBJ);      echo OBJECTS directory removed ; fi
	@if [ -d $(PREBUILD)  ]; then rm -r $(PREBUILD); echo PRECOMP directory removed ; fi
cleandebug: clean
	@if [ -f debug        ]; then rm debug         ;echo EXECUT. file       removed ; fi
cleanobj: clean
	@if [ -f $(EXE)       ]; then rm $(EXE)        ;echo EXECUT. file       removed ; fi

all: $(EXE)

.DEFAULT_GOAL := all
