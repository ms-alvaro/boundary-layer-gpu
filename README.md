# boundary-layer-gpu

**A GPU-native incompressible boundary-layer DNS solver in CUDA Fortran** —
single-GPU by default, multi-GPU capable, with a validation harness that
gates every change on *same physics, measured speedup*.

The solver computes spatially-developing flat-plate boundary layers:
laminar, transitional (Tollmien-Schlichting modes or freestream-turbulence
bypass transition via a precursor-HIT inflow library), and turbulent.
2nd-order staggered finite differences, explicit RK2/RK3, fractional-step
projection with an FFT-diagonalized Poisson solve (half-length real-DCT in
x, Fourier in z, tridiagonal Thomas in y) — everything resident on the GPU;
the host is touched only for gated statistics, monitoring, and I/O.

![TBL with HIT freestream-turbulence inflow](docs/img/hit_tbl_snapshot.png)
*Boundary layer under 3% freestream turbulence from a precursor-HIT inflow
library, running on two A100 z-slabs: near-wall Klebanoff streaks (top)
and freestream eddies over the growing layer (bottom).*

## Highlights

- **Native CUDA Fortran** (no OpenACC, no FFTW, no LAPACK — just nvfortran
  + cuFFT), tuned for A100 (`cc80`), fused stencil/transform kernels,
  half-length real-DCT Poisson pipeline.
- **Bitwise-deterministic**: no atomics; fixed-order reductions. Two runs
  of the same case produce identical bits — which is also how the code is
  validated: field-level comparisons at 1e-15 tolerances.
- **Validated against its reference implementation** at machine precision
  (max relative field difference ~1e-14 across laminar, TS-transition and
  HIT-inflow cases), and at the physics level over a full 500k-step
  transition run (identical transition-onset location at every checkpoint):

  ![Transition physics validation](docs/img/validation_transition_physics.png)

- **Fast**: 2.8x the OpenACC reference on the laminar case, 1.33x at
  production cadence on the transition grid, on A100-PCIE-40GB.
- **Multi-GPU (z-slab decomposition, CUDA-aware MPI)**: distributed
  transpose Poisson, halos, statistics, box output and restart — validated
  to ~1e-14 against single-GPU at 2 and 4 ranks. On PCIe-only nodes this
  buys grid *capacity* (270M-cell cases that cannot fit on one 40 GB card),
  not turnaround; see the honest analysis in
  [docs/MULTI_GPU_DESIGN.md](docs/MULTI_GPU_DESIGN.md).
- **Production features**: time-resolved HIT plane-library inflow with
  no-recycling guard, custom wall-normal grids, zero-shear top boundary,
  float32 subvolume ("box") output for data campaigns, exact chain
  restarts (bitwise on one GPU; snapshots interchange across rank counts).

## Repository layout

```
src/          the CUDA Fortran solver          -> boundary_layer_cuda
legacy/       OpenACC + MPI reference implementation (validation baseline)
cases/        runnable cases: laminar/, transition/, turbulent_lund/,
              profiles/ (mean-profile data used by inflow generators)
validation/   the harness: physics pass/fail, field-level diffs, timing
docs/         design documents and validation records
```

## Build and run

Requires the NVIDIA HPC SDK (nvfortran, 24.3+) and a cc80+ GPU.

```bash
cd src && make                    # -> ./boundary_layer_cuda at the repo root
cd ..
mpirun -np 1 ./boundary_layer_cuda -i cases/laminar/laminar.turbb
```

`make DEBUG=1` builds with bounds checks and per-kernel launch-error
checking. Multi-GPU (experimental; requires `nx-2` and `nz-2` divisible by
the rank count, and UCX CUDA transports):

```bash
CUDA_VISIBLE_DEVICES=0,1 mpirun -np 2 --mca pml ucx \
  -x UCX_TLS=self,sm,cuda_copy,cuda_ipc -x UCX_MEMTYPE_CACHE=n \
  ./boundary_layer_cuda -i case.turbb
```

The legacy OpenACC build (`cd legacy && make`) additionally needs MPI,
FFTW built with nvfortran, and LAPACK; it is kept as the harness baseline.

**Grid-size tip:** the Poisson FFT lengths are `nx-2` and `nz-2` — choose
them with small prime factors (2, 3, 5, 7). A large prime factor (e.g.
`812 = 4*7*29`) forces cuFFT onto a slow path and costs ~10% per step.

## Validation

Every change is gated by the harness (see
[validation/README.md](validation/README.md) for the method and the full
baseline/result tables):

```bash
python3 validation/validate.py laminar --exe ./boundary_layer_cuda --label my-change
python3 validation/validate.py compare --a openacc-ref --b my-change
python3 validation/validate.py bench   --exe ./boundary_layer_cuda --label my-change
```

- `laminar`: Blasius case with numeric checks (Cf vs the similarity
  solution and a golden value, divergence at machine zero, freestream
  cleanliness) — runs in seconds.
- `compare`: field-by-field diff between two builds/runs + speedup.
- `bench`: production-sized grid, ns/cell/step.

## Cases

| case | what it is |
|---|---|
| `cases/laminar` | Blasius flat plate; the correctness gate (~8 s on an A100) |
| `cases/transition` | K-type transition from Tollmien-Schlichting inflow modes |
| `cases/turbulent_lund` | turbulent BL with Lund recycling inflow (legacy solver only; needs an external restart field) |

Input format: `key = value` text files (`.turbb`); see the commented
examples in each case directory and the parameter table in
`legacy/input_parameters.turbb`. The CUDA solver covers the DNS path
(RK2/RK3, Blasius/TS-mode/HIT inflows, Dirichlet or zero-shear top,
no-slip wall, restart, box output); LES and wall models currently run only
in `legacy/` and the solver stops with a clear message if one is requested.

## Documentation

- [docs/CUDA_PORT_DESIGN.md](docs/CUDA_PORT_DESIGN.md) — the port's design
  contract: module layout, verbatim-numerics and determinism conventions.
- [docs/MULTI_GPU_DESIGN.md](docs/MULTI_GPU_DESIGN.md) — the multi-GPU
  architecture, implementation status, and measured scaling (including
  what does *not* scale on PCIe and why).
- [validation/README.md](validation/README.md) — validation methodology
  and the complete record of baselines and results.

## Performance summary (A100-PCIE-40GB, nvfortran 24.3)

| case | legacy OpenACC | CUDA solver |
|---|---|---|
| laminar 302x64x8 | 2.35 ms/step | **0.83 ms/step (2.8x)** |
| transition 814x125x66, production cadence | 13.2 ms/step | **10.0 ms/step (1.33x)** |
| capacity: 3074x341x258 (270M cells, ~65 GB) | does not fit | runs on 4 GPUs |

## Credits

Original CPU/MPI solver by **Adrian Lozano-Duran** (Computational
Turbulence Group, MIT). OpenACC port, CUDA Fortran rewrite, and multi-GPU
implementation by **Alvaro Martinez Sanchez**.

## License

Not yet licensed for redistribution — contact the authors.
