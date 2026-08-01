# Makefile for the CUDA Fortran boundary-layer solver.
#
#   make            optimized build  -> ../boundary_layer_cuda
#   make DEBUG=1    -O0 -g, bounds checks, per-kernel launch-error checks
#   make clean
#
# Requires the NVIDIA HPC SDK (nvfortran). No MPI, no FFTW, no LAPACK —
# the only library dependency is cuFFT.

NVHPC_ROOT ?= /opt/nvidia/hpc_sdk/Linux_x86_64/24.3
# mpifort (NVHPC OpenMPI wrapper around nvfortran): required since P5.1 for
# the CUDA-aware MPI layer (mpi_mod). Single-GPU runs work unchanged and can
# still be launched WITHOUT mpirun (OpenMPI singleton mode).
FC          = $(NVHPC_ROOT)/comm_libs/mpi/bin/mpifort
# mpifort needs the underlying nvfortran on PATH
export PATH := $(NVHPC_ROOT)/compilers/bin:$(PATH)

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
      mpi_mod.f90       \
      param_mod.f90     \
      grid_mod.f90      \
      field_mod.f90     \
      reductions.cuf    \
      ic_inflow_mod.f90 \
      lund_inflow_mod.f90 \
      les_mod.cuf       \
      legacy_mpi.f90    \
      legacy_global.f90 \
      legacy_interpolation.f90 \
      legacy_bc.f90     \
      legacy_subgrid.f90 \
      legacy_newton.f90 \
      legacy_wallmodel.f90 \
      wm_bridge.f90     \
      hit_inflow_mod.cuf \
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

# les_mod: FMA off — the dynamic-Smagorinsky Germano ratio amplifies
# FMA-vs-serial rounding differences at ill-conditioned (near-laminar)
# points to O(1) local nu_t differences vs the CPU reference; nofma
# restores bit-level agreement with the gfortran baseline.
les_mod.o: les_mod.cuf
	$(FC) $(FFLAGS) -gpu=cc80,nofma -c $< -o $@

%.o: %.cuf
	$(FC) $(FFLAGS) -c $< -o $@

# Module dependencies (each .o is the proxy for its .mod). Without these,
# editing a module rebuilt only its own object and left downstream objects
# referencing a stale module layout — which manifested as corrupted module
# variables at runtime (LES_model guard tripping on LES=0 input).
mpi_mod.o:        precision_mod.o
param_mod.o:      precision_mod.o mpi_mod.o
grid_mod.o:       precision_mod.o param_mod.o mpi_mod.o
field_mod.o:      precision_mod.o param_mod.o grid_mod.o
reductions.o:     precision_mod.o
ic_inflow_mod.o:  precision_mod.o param_mod.o grid_mod.o field_mod.o
lund_inflow_mod.o: precision_mod.o param_mod.o grid_mod.o field_mod.o \
                  mpi_mod.o ic_inflow_mod.o
les_mod.o:        precision_mod.o param_mod.o grid_mod.o field_mod.o
legacy_global.o:  legacy_mpi.o
legacy_interpolation.o: legacy_global.o
legacy_bc.o:      legacy_global.o legacy_mpi.o
legacy_subgrid.o: legacy_global.o legacy_interpolation.o legacy_bc.o legacy_mpi.o
legacy_newton.o:  legacy_global.o
legacy_wallmodel.o: legacy_global.o legacy_interpolation.o legacy_subgrid.o \
                  legacy_bc.o legacy_newton.o
wm_bridge.o:      precision_mod.o param_mod.o grid_mod.o field_mod.o \
                  mpi_mod.o les_mod.o legacy_global.o legacy_mpi.o \
                  legacy_wallmodel.o
hit_inflow_mod.o: precision_mod.o param_mod.o grid_mod.o ic_inflow_mod.o
rhs_kernels.o:    precision_mod.o param_mod.o grid_mod.o field_mod.o les_mod.o
bc_kernels.o:     precision_mod.o param_mod.o grid_mod.o field_mod.o \
                  ic_inflow_mod.o hit_inflow_mod.o reductions.o mpi_mod.o \
                  wm_bridge.o
poisson_mod.o:    precision_mod.o grid_mod.o field_mod.o reductions.o mpi_mod.o
timestep_mod.o:   precision_mod.o param_mod.o grid_mod.o field_mod.o \
                  reductions.o rhs_kernels.o bc_kernels.o poisson_mod.o mpi_mod.o \
                  lund_inflow_mod.o les_mod.o
io_mod.o:         precision_mod.o param_mod.o grid_mod.o field_mod.o \
                  timestep_mod.o poisson_mod.o mpi_mod.o lund_inflow_mod.o les_mod.o
main.o:           param_mod.o grid_mod.o field_mod.o ic_inflow_mod.o mpi_mod.o \
                  hit_inflow_mod.o poisson_mod.o timestep_mod.o \
                  bc_kernels.o reductions.o io_mod.o lund_inflow_mod.o les_mod.o

clean:
	rm -f $(OBJ) *.mod $(EXE)

.PHONY: clean
