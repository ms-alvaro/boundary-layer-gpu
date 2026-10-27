# Boundary layer solver

An incompressible boundary-layer DNS and LES solver written in CUDA Fortran.
It runs on one GPU by default and on several GPUs with a spanwise slab
decomposition. A test harness checks every change for the same physics and
records the measured speed.

The solver computes spatially developing flat-plate boundary layers:
laminar, transitional (Tollmien-Schlichting modes or bypass transition under
free-stream turbulence from a precursor inflow library) and turbulent. The
discretization is second-order staggered finite differences with explicit
Euler, RK2 or RK3 stepping and a fractional-step projection. The Poisson
equation is diagonalized by a half-length real cosine transform in $x$ and a
Fourier transform in $z$, and solved by the Thomas algorithm in $y$. The
fields stay resident on the GPU; the host is used for gated statistics,
monitoring and I/O.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/img/tbl_hit_transition_dark.gif">
  <img alt="Bypass transition under freestream turbulence" src="docs/img/tbl_hit_transition.gif">
</picture>

*Bypass transition of a flat-plate boundary layer under 3% free-stream
turbulence (precursor inflow library), computed on two A100 z-slabs:
$u'/U_\infty$ in a wall-normal $x$-$y$ plane at equal aspect ratio, axes in
inlet-$\delta_{99}$ units. Free-stream turbulence above drives streaks and
turbulent spots in the layer below.*

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/img/tbl_hit_transition_les_dark.gif">
  <img alt="Same case as wall-resolved LES" src="docs/img/tbl_hit_transition_les.gif">
</picture>

*The same case computed as a wall-resolved LES (dynamic Smagorinsky,
`LES = 2`) on a 12x coarser grid with double the time step, shown at the
same instants as the DNS above. The streak growth and breakdown follow the
same scenario at about 8x lower compute cost: one GPU for 51 minutes against
two GPUs for 3.4 hours.*

## Repository layout

```
src/          the CUDA Fortran solver          -> boundary_layer_cuda
cases/        runnable cases: laminar/, transition/, turbulent_lund/,
              profiles/ (mean-profile data used by inflow generators)
tests/        try the code + the gating harness: physics pass/fail,
              field-level diffs, timing
docs/         design documents, parameter reference, validation records
```

The earlier CPU (MPI) and OpenACC implementations are kept in git history
on the `gpu-openacc` branch. `main` is CUDA-only.

## Build and run

Requires the NVIDIA HPC SDK (nvfortran, 24.3+) and a cc80+ GPU.

```bash
cd src && make                    # -> ./boundary_layer_cuda at the repo root
cd ..
mpirun -np 1 ./boundary_layer_cuda -i cases/laminar/laminar.turbb
```

`make DEBUG=1` builds with bounds checks and per-kernel launch-error
checking. Multi-GPU runs are experimental. They require `nx-2` and `nz-2`
divisible by the rank count and the UCX CUDA transports:

```bash
CUDA_VISIBLE_DEVICES=0,1 mpirun -np 2 --mca pml ucx \
  -x UCX_TLS=self,sm,cuda_copy,cuda_ipc -x UCX_MEMTYPE_CACHE=n \
  ./boundary_layer_cuda -i case.turbb
```

The Poisson FFT lengths are `nx-2` and `nz-2`. Choose them with small prime
factors (2, 3, 5, 7). A large prime factor, for example `812 = 4*7*29`,
sends cuFFT onto a slow path and costs about 10% per step.

## Tests

The fastest end-to-end check is the laminar Blasius case, about 8 s on an
A100. The same harness gates every solver change; see
[tests/README.md](tests/README.md) for the method and the full baseline and
result tables.

```bash
python3 tests/validate.py laminar --exe ./boundary_layer_cuda --label my-run
python3 tests/validate.py compare --a openacc-ref --b my-run
python3 tests/validate.py bench   --exe ./boundary_layer_cuda --label my-run
```

- `laminar`: Blasius case with numeric checks (Cf against the similarity
  solution and a golden value, divergence at machine zero, free-stream
  cleanliness). It runs in seconds.
- `compare`: field-by-field diff between two builds or runs, plus speedup.
- `bench`: production-sized grid, ns/cell/step.

## Cases

| case | what it is |
|---|---|
| `cases/laminar` | Blasius flat plate; the correctness gate (about 8 s on an A100) |
| `cases/transition` | K-type transition from Tollmien-Schlichting inflow modes |
| `cases/turbulent_lund` | turbulent boundary layer with Lund recycling inflow (`inflow_flag = 3/5`), ported and validated to machine precision against the CPU solver; see [cases/turbulent_lund/README.md](cases/turbulent_lund/README.md) |

Input files are `key = value` text files (`.turbb`); see the commented
examples in each case directory and the parameter table in
[docs/input_parameters.turbb](docs/input_parameters.turbb). String values
must follow `=` without a blank, as the shipped cases do.

## Physics and models

The solver integrates the incompressible Navier-Stokes equations with
$U_\infty = 1$ and the inlet distance from the leading edge as the unit of
length, so the inlet Reynolds number is $1/\texttt{nu}$. The inflow can be a
Blasius profile, Blasius plus temporal Tollmien-Schlichting modes, Blasius
read from a file, Blasius plus random noise, Lund recycling of a downstream
plane (full or fluctuations only) or Blasius plus time-resolved planes from
a precursor homogeneous-turbulence library (`inflow_flag` 1 to 6). The top
boundary is a Dirichlet free stream, a Coleman or Abe blowing-and-suction
lid, a Falkner-Skan free stream or a zero-shear lid (`top_flag` 0 to 4).
The outflow is convective with a global mass correction, the spanwise
direction is periodic, and the wall is no-slip unless a wall model is
active. The LES models are constant-coefficient Smagorinsky and dynamic
Smagorinsky with spanwise averaging or pointwise coefficients (`LES` 1 to
3), all on the GPU. The wall models are the 14 options of the reference
solver (`WM` 1 to 9 and 11 to 15: slip-length Robin conditions from
momentum balances, a dynamic slip model, a least-squares Robin fit, log-law
models and skin-friction correlations) together with a log-law wall-stress
estimate (`TauwModel = 1`); the reference wall-model code runs on the host
and is bridged to the GPU at every substep. Seven of the wall models are
validated against the CPU solver to machine precision; the others run
identically to the reference on the probes tried but have not been exercised
to a physical result. The Euler, RK2 and RK3 schemes, restart, statistics,
box output and strided 3-D snapshots are all available. Each option, its
input keys, the equations it imposes and its validation status are in
[docs/PHYSICS_AND_MODELS.md](docs/PHYSICS_AND_MODELS.md). Two memory errors
found in the reference wall models during the port are repaired (see the
commit history).

## Documentation

- [docs/PHYSICS_AND_MODELS.md](docs/PHYSICS_AND_MODELS.md): governing
  equations and numerics, every boundary condition, LES and wall-model
  option with its input keys, output formats, and the features that exist
  only on branches.
- [docs/CUDA_PORT_DESIGN.md](docs/CUDA_PORT_DESIGN.md): the port's design
  contract: module layout, verbatim-numerics and determinism conventions.
- [docs/MULTI_GPU_DESIGN.md](docs/MULTI_GPU_DESIGN.md): the multi-GPU
  architecture, implementation status and measured scaling, including what
  does not scale on PCIe and why.
- [tests/README.md](tests/README.md): test and validation method and the
  record of baselines and results.

## Performance summary (A100-PCIE-40GB, nvfortran 24.3)

Speedup of one A100 GPU over the original CPU solver (16-core MPI):

| case | CPU solver (MPI) | CUDA solver | speedup |
|---|---|---|---|
| transition 814x125x66, production cadence | 87 ms/step (16 cores) | 10.0 ms/step (1 GPU) | 8.7x |
| laminar 302x64x8 | 25.5 ms/step (6 cores*) | 0.83 ms/step (1 GPU) | 31x |
| capacity: 3074x341x258 (270M cells, ~65 GB) | 28.1 s/step (16 cores) | 0.7-0.8 s/step (4 GPUs) | ~35-40x |

*The z-slab decomposition requires `nz-2` divisible by the rank count, so
the laminar case (nz = 8) caps at 6 ranks. CPU numbers are median
per-monitor-interval timings of the pre-port MPI solver rebuilt from this
repository's history (gfortran -O2 + FFTW 3.3.10), on the same machine as
the GPU runs (dual-socket 28-core Xeon host of the A100s).*

## Credits

Original CPU/MPI solver by Adrian Lozano-Duran (Computational Turbulence
Group, MIT). OpenACC port, CUDA Fortran rewrite and multi-GPU implementation
by Alvaro Martinez Sanchez.

## License

Not yet licensed for redistribution. Contact the authors.
