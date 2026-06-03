# Boundary layer GPU solver

Incompressible boundary layer DNS solver accelerated with NVIDIA GPUs.
2nd-order finite differences on a staggered mesh, explicit RK2/RK3 time
integration, fractional-step pressure projection. The pressure Poisson
equation is solved via cosine transform in x (Neumann BC), Fourier
transform in z (periodic), and tridiagonal Thomas solve in y — all on
the GPU using cuFFT and OpenACC.

Written by Adrian Lozano-Duran (CPU/MPI version) and Alvaro Martinez
Sanchez (GPU port).

## Performance

**22x faster than 16 CPU cores** on a single NVIDIA A100.

| Metric | Value |
|---|---|
| Time per cell per step | 2.09 ns |
| Speedup vs 16 CPU cores | 22x |
| Speedup vs 1 CPU core | 359x |
| Test grid | 814 x 257 x 193 (40.4M cells) |
| GPU | A100 80 GB SXM |

## Prerequisites

- **NVIDIA HPC SDK** 24.3+ (`nvfortran`, `mpifort`)
- **NVIDIA GPU** with compute capability 8.0+ (A100, H100, etc.)
- **FFTW 3.3+** compiled with the NVHPC Fortran compiler (for the
  MPI-parallel FFTW fallback in the pressure solver)
- **CUDA toolkit** (included with HPC SDK)
- **LAPACK** (included with HPC SDK)

## Compilation

1. Set environment variables:

```bash
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3  # adjust to your installation
export FFTW_DIR=/path/to/fftw-3.3.10-nvhpc                # FFTW built with nvfortran
```

2. Compile:

```bash
make clean && make
```

This produces the `boundary_layer_gpu` executable.

For debugging:

```bash
make clean && make debug
```

> **Note:** `statistics.f90` is compiled at `-O0` due to an internal
> compiler error in nvfortran 24.3. This has no measurable performance
> impact since statistics are computed infrequently.

### Building FFTW with NVHPC

If your system does not have FFTW compiled with nvfortran:

```bash
export FC=$NVHPC_ROOT/compilers/bin/nvfortran
export CC=$NVHPC_ROOT/compilers/bin/nvc
./configure --prefix=$HOME/opt/fftw-3.3.10-nvhpc --enable-mpi FC=$FC CC=$CC
make -j && make install
```

### CPU-only build

A CPU Makefile (`Makefile.cpu`) is included for reference. It uses
`mpifort` (gfortran) with FFTW and LAPACK:

```bash
make -f Makefile.cpu FFTW_DIR=/path/to/fftw LAPACK_DIR=/path/to/lapack
```

## Running

```bash
mpirun -np 1 ./boundary_layer_gpu -i <input_file>.turbb
```

Currently single-GPU only (MPI decomposition is in z, but the GPU solver
runs on rank 0).

## Test cases

### Laminar Blasius (`laminar_test/`)

Laminar flat-plate boundary layer validation. Compares skin friction and
velocity profiles against the Blasius similarity solution.

```bash
cd laminar_test
python3 generate_blasius.py    # generate blasius_solution.dat
mpirun -np 1 ../boundary_layer_gpu -i laminar.turbb
python3 postprocess.py         # validation plots
```

### TS-mode transition (`transition_test/`)

Boundary layer transition triggered by Tollmien-Schlichting modes at the
inflow (inflow_flag=1).

```bash
cd transition_test
python3 generate_temporal_modes_local.py  # generate TS mode file
mpirun -np 1 ../boundary_layer_gpu -i transition.turbb
python3 plot_snapshot.py data/BL_transition.NNNNN
```

## Input parameters

See `input_parameters.turbb` for a documented template. Key parameters:

| Parameter | Description |
|---|---|
| `nxyz` | Grid points (nx, ny, nz) |
| `boxsize` | Domain size (Lx, Ly, Lz, alpha_stretch) |
| `CFL` | CFL number (positive) or fixed dt (negative) |
| `nu` | Kinematic viscosity |
| `RKscheme` | Time integration: 1=Euler, 2=RK2, 3=RK3 |
| `LES` | SGS model: 0=DNS, 1=Smagorinsky, 2=DSM |
| `WM` | Wall model: 0=none, 1-9=various |
| `inflow_flag` | Inflow BC: 1=Blasius, 3=Lund, 4=Blasius+noise, 6=Blasius+HIT |
| `nsteps` | Total time steps |
| `nsave` | Save snapshot every N steps |

## Code structure

| File | Description |
|---|---|
| `main.f90` | Main program, OpenACC data region |
| `global.f90` | Global variable declarations |
| `initialization.f90` | Grid setup, Thomas LU precomputation, cuFFT plans |
| `equations.f90` | RHS computation (inlined convective + viscous terms) |
| `time_integration.f90` | RK2/RK3 stepping with OpenACC loops |
| `projection.f90` | GPU pressure solver (cuFFT + Thomas) |
| `cufft_solver.f90` | cuFFT plan creation (`cufftPlanMany`) |
| `pressure.f90` | Pressure Poisson equation setup |
| `boundary_conditions.f90` | Inflow/outflow/wall BCs (GPU-resident) |
| `rescaled_inlet_bc.f90` | Lund's rescaling for turbulent inflow |
| `statistics.f90` | Time-averaged flow statistics |
| `monitor.f90` | Runtime diagnostics (Cf, max U, divergence) |
| `input_output.f90` | Snapshot read/write |
| `params.f90` | Input parameter parsing |
| `fftz.f90` | FFT wrappers |
| `interpolation.f90` | Grid interpolation utilities |
| `subgrid.f90` | LES subgrid-scale models |
| `wallmodel.f90` | Wall stress models |
| `miscel.f90` | Miscellaneous utilities |
| `Newton_solver.f90` | Newton solver (for wall models) |
| `finalization.f90` | Cleanup |
| `mpi.f90` | MPI initialization |

## License

Original code by Adrian Lozano-Duran and the Computational Turbulence
Group at MIT.
