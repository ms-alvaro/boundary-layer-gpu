# Boundary layer GPU solver

Incompressible boundary-layer DNS solver written natively in **CUDA
Fortran** for NVIDIA GPUs. 2nd-order finite differences on a staggered
mesh, explicit RK2/RK3 time integration, fractional-step pressure
projection. The pressure Poisson equation is solved via cosine transform
in x (Neumann BC), Fourier transform in z (periodic), and tridiagonal
Thomas solve in y — all on the GPU using cuFFT. The full timestep is
GPU-resident: the host is touched only for gated statistics, monitoring,
and snapshot I/O.

Original CPU/MPI solver by Adrian Lozano-Duran (Computational Turbulence
Group, MIT); OpenACC port and CUDA Fortran rewrite by Alvaro Martinez
Sanchez.

## Repository layout

```
src/          CUDA Fortran solver (the active code) -> boundary_layer_cuda
legacy/       OpenACC + MPI reference implementation -> boundary_layer_gpu
validation/   automated harness: physics pass/fail, field-level diffs, timing
docs/         design documents (docs/CUDA_PORT_DESIGN.md)
laminar_test/     Blasius validation case
transition_test/  TS-mode transition case
tests/            turbulent Lund-inflow case (needs external restart file)
```

## Performance (A100-PCIE-40GB, nvhpc 24.3)

| build | laminar 302x64x8 | transition grid 814x125x66 (production cadence) |
|---|---|---|
| legacy OpenACC     | 2.35 ms/step | 13.2 ms/step |
| **CUDA (src/)**    | **0.83 ms/step (2.8x)** | **10.0 ms/step (1.32x)** |

The CUDA build reproduces the OpenACC reference to max relative field
difference ~4e-12 after 2000 laminar steps (and 2e-14 over 5000
transition steps with active TS-mode inflow) and is bitwise run-to-run
deterministic (no atomics; fixed-order reductions). Optimizations so
far: fused RHS/pack kernels, save+advance fusion, fully device-resident
mass-conservation reduction (no per-step synchronization), half-length
real-DCT Poisson pipeline.

Legacy reference numbers: 22x vs 16 CPU cores, 2.09 ns/cell/step at
814x257x193 on an A100-80GB-SXM.

**Grid-size tip:** the Poisson solve FFTs have length `nx-2` (x) and
`nz-2` (z). Choose grids so these factor into small primes (2, 3, 5, 7).
Example: `nx = 812` (810 = 2*3^4*5) runs the same physics ~11% faster
per step than `nx = 814` (812 = 4*7*29, which forces cuFFT through a
slow large-prime path).

## Prerequisites

- NVIDIA HPC SDK 24.3+ (`nvfortran`) and a GPU with compute capability
  8.0+ (the Makefile targets `cc80`; adjust `GPU` for H100/H200).
- That's all for the CUDA solver (cuFFT ships with the SDK). The legacy
  build additionally needs MPI + FFTW built with nvfortran + LAPACK.

## Build and run

```bash
cd src && make          # -> ./boundary_layer_cuda at the repo root
mpirun -np 1 ./boundary_layer_cuda -i laminar_test/laminar.turbb
# multi-GPU (P5.2, experimental): correctness-validated (fields match
# single-rank to 9e-15); buys grid CAPACITY beyond one card's 40 GB —
# see docs/MULTI_GPU_DESIGN.md for honest performance notes. Use the
# UCX CUDA transports:
CUDA_VISIBLE_DEVICES=0,1 mpirun -np 2 --mca pml ucx \
  -x UCX_TLS=self,sm,cuda_copy,cuda_ipc -x UCX_MEMTYPE_CACHE=n \
  ./boundary_layer_cuda -i case.turbb
```

(Launch through `mpirun` from the NVHPC `comm_libs` — OpenMPI singleton
mode does not work in this environment.)

`make DEBUG=1` builds -O0 with bounds checks and per-kernel launch-error
checking.

Legacy reference build (used as the comparison baseline by the harness):

```bash
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH=$NVHPC_ROOT/compilers/bin:$NVHPC_ROOT/comm_libs/mpi/bin:$PATH
export FFTW_DIR=$HOME/opt/fftw-3.3.10-nvhpc
cd legacy && make       # -> legacy/boundary_layer_gpu
```

## Validation

Every solver change is gated on the harness (see `validation/README.md`):

```bash
python3 validation/validate.py laminar --exe ./boundary_layer_cuda --label my-change
python3 validation/validate.py compare --a openacc-ref --b my-change
python3 validation/validate.py bench   --exe ./boundary_layer_cuda --label my-change
```

The classic plot-based check is still available:
`cd laminar_test && python3 postprocess.py` after a run.

## Test cases

- **`laminar_test/`** — laminar flat-plate Blasius validation (Cf and
  velocity profiles vs the similarity solution). Runs in seconds.
- **`transition_test/`** — transition triggered by Tollmien-Schlichting
  modes at the inflow; `generate_temporal_modes_local.py` builds the
  mode file, `plot_ic_inflow.py` previews IC + inflow before launching.

## Input parameters

See `legacy/input_parameters.turbb` for a documented template. Key ones:

| Parameter | Description |
|---|---|
| `nxyz` | Grid points (nx, ny, nz) |
| `boxsize` | Domain size (Lx, Ly, Lz, alpha_stretch) |
| `CFL` | CFL number (positive) or fixed dt (negative value = -dt) |
| `nu` | Kinematic viscosity |
| `RKscheme` | Time integration: 2=RK2, 3=RK3 (Euler is legacy-only) |
| `inflow_flag` | Inflow BC (CUDA solver: 1 = Blasius + temporal modes) |
| `nsteps` / `nsave` / `nstats` / `nmonitor` | Run control / output cadence |

The CUDA solver covers the full DNS production path: RK2/RK3,
`inflow_flag=1` (Blasius + temporal modes) and `inflow_flag=6` (Blasius
+ HIT plane library with no-wrap guard), `top_flag=0` (Dirichlet) and
`top_flag=4` (zero-shear), `ygrid_file` custom wall-normal grids,
subvolume `boxout_*` output for causal-analysis campaigns, no-slip
wall, restart, single GPU. LES, wall models, and Lund recycling inflow
currently run only in `legacy/` (next porting phase) — the solver stops
with a clear message if one is requested.

## Code structure (src/)

| File | Description |
|---|---|
| `main.f90` | Driver: setup sequence + time loop |
| `param_mod.f90` | Input-file parsing, run parameters |
| `grid_mod.f90` | Grid generation, metrics, device copies |
| `field_mod.f90` | Device state arrays + host mirrors |
| `ic_inflow_mod.f90` | Blasius IC, inlet/top profiles, inflow-mode tables |
| `rhs_kernels.cuf` | Momentum RHS kernels (convective + viscous) |
| `bc_kernels.cuf` | Boundary-condition kernels + mass conservation |
| `poisson_mod.cuf` | cuFFT plans, DCT/FFT/Thomas Poisson solve, projection |
| `timestep_mod.cuf` | RK2/RK3 drivers, advance kernels, CFL dt |
| `reductions.cuf` | Deterministic (fixed-order, atomic-free) reductions |
| `io_mod.f90` | Snapshots, statistics, monitor, restart |
| `precision_mod.f90` | Working precision, small helpers |

Design rationale and porting conventions: `docs/CUDA_PORT_DESIGN.md`.

## License

Original code by Adrian Lozano-Duran and the Computational Turbulence
Group at MIT.
