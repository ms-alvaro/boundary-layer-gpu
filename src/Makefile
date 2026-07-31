# Makefile for the CUDA Fortran boundary-layer solver.
#
#   make            optimized build  -> ../boundary_layer_cuda
#   make DEBUG=1    -O0 -g, bounds checks, per-kernel launch-error checks
#   make clean
#
# Requires the NVIDIA HPC SDK (nvfortran). No MPI, no FFTW, no LAPACK —
# the only library dependency is cuFFT.

NVHPC_ROOT ?= /opt/nvidia/hpc_sdk/Linux_x86_64/24.3
FC          = $(NVHPC_ROOT)/compilers/bin/nvfortran

EXE     = ../boundary_layer_cuda
GPU     = -gpu=cc80
FFLAGS  = -cuda $(GPU) -Mpreprocess
LDFLAGS = -cuda $(GPU) -cudalib=cufft

ifeq ($(DEBUG),1)
  FFLAGS += -O0 -g -Mbounds -DDEBUG_KERNELS
else
  FFLAGS += -O2
endif

# Dependency order (each module uses only earlier ones)
SRC = precision_mod.f90 \
      param_mod.f90     \
      grid_mod.f90      \
      field_mod.f90     \
      reductions.cuf    \
      ic_inflow_mod.f90 \
      rhs_kernels.cuf   \
      bc_kernels.cuf    \
      poisson_mod.cuf   \
      timestep_mod.cuf  \
      io_mod.f90        \
      main.f90

OBJ = $(addsuffix .o,$(basename $(SRC)))

$(EXE): $(OBJ)
	$(FC) -o $@ $(OBJ) $(LDFLAGS)
	@echo "$(EXE) built."

# nvfortran 24.3 ICEs (fort1 signal 11) on io_mod.f90 at -O2 — same compiler
# bug family as the reference's statistics.f90 -O0 workaround. io_mod is
# host-only, gated I/O; -O1 compiles fine and costs nothing measurable.
io_mod.o: io_mod.f90
	$(FC) $(filter-out -O2,$(FFLAGS)) -O1 -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cuf
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) *.mod $(EXE)

.PHONY: clean
