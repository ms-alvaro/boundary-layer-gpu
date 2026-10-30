# CUDA Fortran port: design document

Status: Phase 2 (DNS timestep path). This document is the contract between
modules. The names, shapes and conventions defined here are binding.

This document is kept as the design record of Phase 2 and describes the
scope at that time. The features listed below as not ported (LES, wall
models, Lund recycling inflow, the Falkner-Skan and blowing-suction top
conditions, Euler time stepping and multi-rank execution) have since been
ported. The current state of the solver is described in
[PHYSICS_AND_MODELS.md](PHYSICS_AND_MODELS.md) and in the option tables of
the top-level README, and the multi-GPU implementation in
[MULTI_GPU_DESIGN.md](MULTI_GPU_DESIGN.md). The harness directory called
`validation/` below is now `tests/`.

## Goal

Replace the OpenACC solver with a native CUDA Fortran implementation that is

1. correct: the laminar harness passes, and the final fields match the
   OpenACC reference within max relative diff 1e-8 over 2000 steps
   (`validation/validate.py compare`),
2. organized: one module per concern, kernels documented with the
   discrete equations they implement, no dead code,
3. fast: never slower than the reference at any merge point. Phase 3
   then optimizes (fusion, real-DCT, streams) on top of this port.

## Scope of Phase 2 (and explicit non-goals)

Ported: DNS RHS, RK2/RK3, GPU boundary conditions (`inflow_flag=1` including
temporal modes, convective outflow, periodic z, Dirichlet top, no-slip wall,
global mass conservation), cuFFT Poisson + Thomas, projection, monitor,
statistics, snapshot/stats/restart I/O, and the CFL dt as a device
reduction, which improves on the reference, where it runs on the host.

Not ported in Phase 2: LES (`subgrid.f90`), wall models (`wallmodel.f90`),
Lund recycling inflow (flags 3/5), Falkner-Skan and blowing-suction top BCs,
Euler time stepping (GPU-incoherent in the reference, dropped deliberately),
multi-rank (`nprocs > 1` aborts with a clear message), and the true-pressure
solve (`pressure.f90`, dead in the reference). The input parser still
accepts these flags, but the code stops with an informative error if one is
requested.

## Source layout

```
src/
  precision_mod.f90    dp = real64; pi
  param_mod.f90        input-file parsing, runtime parameters   [from params.f90]
  grid_mod.f90         grid generation + metrics, host & device [from initialization.f90]
  field_mod.f90        all device state arrays + host mirrors   [from global.f90]
  ic_inflow_mod.f90    Blasius profile, IC, inflow-mode tables  [from initialization.f90]
  rhs_kernels.cuf      RHS kernels for u, v, w (DNS)            [from equations.f90]
  bc_kernels.cuf       boundary-condition kernels               [from boundary_conditions.f90:127-262]
  poisson_mod.cuf      cuFFT plans, DCT pack/unpack, Thomas,
                       projection                               [from projection.f90, cufft_solver.f90,
                                                                 initialization.f90 (operators)]
  timestep_mod.cuf     RK2/RK3 driver, advance kernels, dt      [from time_integration.f90]
  reductions.cuf       deterministic device reductions (max/sum)
  io_mod.f90           snapshots, stats, restart, monitor       [from input_output.f90,
                                                                 statistics.f90, monitor.f90]
  main.f90             driver
  Makefile             produces ../boundary_layer_cuda
```

`.cuf` files contain `attributes(global)` kernels and launch code. `.f90`
files are host-only. Everything compiles with nvfortran (`-cuda -gpu=cc80`).

## Conventions (binding)

- `use precision_mod, only: dp` everywhere, and all reals are `real(dp)`.
  Literal constants are written in the `1.0_dp` style (or `0.5_dp*...`),
  never `d0` in new expressions. See "verbatim numerics" below.
- Lower-case modern Fortran (`do ... end do`), `implicit none` in every
  module, `private` default with explicit `public` lists.
- Device arrays live in `field_mod` (state) and `grid_mod` (metrics) as
  module-level `device, allocatable` variables with the same names and
  shapes as the reference `global.f90`, suffixed `_d`:
  `U_d(nx,nyg,nzg)`, `V_d(nxg,ny,nzg)`, `W_d(nxg,nyg,nz)`,
  `P_d(nxg,nyg,nzg)`, `Uo_d, Vo_d, Wo_d`, `rhs_p_d`, stage arrays
  `Fu1_d..Fw3_d`, cuFFT buffers `plane_d, rhs_hat_d`, Thomas factors
  `thomas_dl_fact_d, thomas_d_pivot_d, thomas_du_d`, `dct_twiddle_d`, etc.
  Host mirrors (no suffix) exist only where I/O or setup needs them
  (`U, V, W, P` in field_mod, grids in grid_mod).
- Scalar run parameters (nx, nyg, dt, nu, ...) are host module variables in
  `param_mod`/`grid_mod`. Kernels receive them by value as arguments
  (`integer, value :: nx, nyg, nzg`, `real(dp), value :: dt, nu, ...`).
  Array arguments to kernels are declared with explicit shape from those
  value dims. Kernels access big arrays via use-association from
  `field_mod`/`grid_mod` or as arguments. This is the writer's choice, but
  it must be consistent within a file, and arguments are preferred for
  testability.
- Thread mapping for 3D kernels over (i,j,k) with i fastest:
  `tBlock = dim3(64,4,1)`, `grid = dim3(ceiling((ihi-ilo+1)/64.),
  ceiling((jhi-jlo+1)/4.), khi-klo+1)`, and inside the kernel
  `i = ilo + (blockIdx%x-1)*blockDim%x + threadIdx%x - 1` etc., with a
  bounds guard `if (i > ihi .or. j > jhi .or. k > khi) return`.
  Provide this via a small helper in `precision_mod` or per file. Do not
  invent per-kernel exotic configurations in Phase 2 (that is Phase 3 work).
- Verbatim numerics: kernel loop bodies are transcribed from the
  reference loop bodies unchanged, with the same expressions and the same
  operation order. The `d0` literals are converted to `_dp`, but the
  arithmetic structure is not "improved". Interpolation weights, metric
  factors and RK tableau values are copied exactly. This is what makes the
  1e-8 comparison achievable. Optimization comes later, gated by the
  harness.
- Deterministic reductions: no atomics. `reductions.cuf` provides
  grid-level max/sum. Kernel 1 writes one partial per block (fixed-order
  sequential loop within each thread's stride range, then a fixed-order
  shared-memory tree). Kernel 2 (single block) combines the partials in
  index order. Same input gives the same result, every run. Used by: CFL
  dt (max), mass-conservation flux (sum), divergence check (max).
- Error handling: `istat = cudaGetLastError()` after each kernel launch in
  debug builds (`make DEBUG=1` defines `-DDEBUG_KERNELS`). The host aborts
  with the kernel name on failure.
- Monitor output format strings are copied verbatim from `monitor.f90` and
  `time_integration.f90`, because the harness parses them with regular
  expressions: `Mean Cf:     :`, `Maximum divergence          :`,
  `Elapsed time (s)            :`, `Maximum U    :` etc., and the
  `PROFILE n= ...` line with fields `copy= RHS= adv= BC= proj=`.
- I/O binary formats are byte-identical to the reference
  (`input_output.f90`). A snapshot is t, nu (8B each), magic -73 (int4),
  istep, nx, x(nx), ny, y(ny), nz, z(nz), nxm, xm, nym, ym, nzm, zm, then
  U, V, W each as n1, n2, n3 (int4 x3) + n1*n2*(n3-1) doubles (stream
  access, little endian, matching what the reference writes on this
  machine). Stats and restart formats are likewise copied from the
  reference code.

## Algorithm reference (what each module implements)

The numerical method is unchanged from the reference:

- Staggered 2nd-order FD, u(nx,nyg,nzg) at x-faces, v at y-faces, w at
  z-faces, p at centers. i is the fastest (contiguous) index. The interior
  is 2..n-1, and planes 1 and n are ghosts or boundaries. The last x-plane
  duplicates nx-1.
- RK substep (s = 1..2 for RK2): save Uo, Vo, Wo once per step. Per
  substep: RHS (divergence-form convection + viscous Laplacian + dPdx),
  then advance `U = Uo + dt*sum_q rk_coef(s,q)*Fu_q` (interior), then BCs,
  projection, BCs.
- Projection: rhs_p = div(u*) on centers, not divided by dt. The dt is
  folded into the eigenvalue scaling exactly as in the reference. The
  Poisson solve is a DCT-II in x (symmetric-extension double-length complex
  trick), a DFT in z, and a tridiagonal solve in y with precomputed Thomas
  factors. Then u = u* - grad(phi). P is saved from rhs_p on substep 1
  with the reference's `rhs_p/(dt*rk2_coef(1,1))` convention.
- BCs (GPU path semantics, boundary_conditions.f90:127-262): inflow
  plane from Blasius + temporal Fourier modes, convective outflow using
  Uo and rk_t(s)*dt, z-periodicity, Dirichlet top from U_top/V_top/W_top,
  no-slip wall (U, W antisymmetric ghosts, V=0 at the wall face), and
  global mass conservation (net-flux reduction, then a Delta_U correction
  on the outflow plane).
- cuFFT: three Z2Z plans made once (2D batched forward, strided 1D z
  inverse, 1D x inverse), executed on device pointers. These are plain
  calls on `device` arrays, so no host_data is needed in CUDA Fortran.

Consult the reference file listed per module in the layout table. When this
document and the reference disagree on numerics, the reference wins, and
the discrepancy is reported.

## Build & validate

```bash
cd src && make            # -> ../boundary_layer_cuda
cd .. && python3 validation/validate.py laminar --exe ./boundary_layer_cuda --label cuda-v1
python3 validation/validate.py compare --a openacc-ref --b cuda-v1
python3 validation/validate.py bench   --exe ./boundary_layer_cuda --label cuda-v1
```

Merge gate: all laminar checks pass, compare tol 1e-8, bench not slower
than openacc-ref (2.31 ns/cell/step).
