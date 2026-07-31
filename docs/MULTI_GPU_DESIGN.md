# Multi-GPU design (Phase 5 proposal — not yet implemented)

Status: design document. The CUDA solver (src/) is single-GPU; this is the
plan for distributing it. Nothing here is built yet.

## Why

1. **Memory**: production-scale grids outgrow one A100-40GB. The solver
   holds ~30 full-size double arrays; a 3074x341x258 grid (270M cells,
   the fst_tbl_ctr scale) needs ~65 GB of state — impossible on one
   40 GB card, comfortable on 4 (~16 GB/rank).
2. **Turnaround**: 5 A100s on tifa sit mostly idle while one runs.

## Hardware reality (tifa)

5x A100-PCIE-40GB, **no NVLink** — peer-to-peer goes over PCIe Gen4
(~20-25 GB/s effective per direction, shared lanes). Communication is
therefore the scarce resource; the design minimizes global data motion
and overlaps what remains.

## Decomposition: z-slabs

Split the spanwise (periodic) direction across ranks, as the legacy MPI
code did (`mod(nz-2, nprocs) = 0` constraint — note nz = 258 means 4
ranks, not 5).

Why z and not x/y:
- z is periodic: halo logic is uniform, no special inflow/outflow ranks.
- y hosts the Thomas solve (serial recurrence) — splitting y would
  serialize across ranks.
- x works too (and balances better for long domains), but the transform
  bookkeeping matches the existing single-GPU code best with x kept
  local (the half-length DCT stays entirely on-rank).

Per-rank state: U(nx, nyg, nzl+2) etc., nzl = (nz-2)/P local planes plus
one halo plane per side.

## What communicates, and how much

Per RK2 substep (2 substeps/step):

| operation | pattern | volume/rank (transition grid, P=4) |
|---|---|---|
| RHS halos (U,V,W, 1 plane/side) | neighbor sendrecv | ~0.8 MB x 6 -> ~5 MB, <1 ms |
| mass-conservation flux | 4-scalar allreduce | negligible (latency only) |
| Poisson forward transpose | all-to-all | ~1/P of the 103 MB mode array |
| Poisson inverse transpose | all-to-all | same |

The Poisson solve is the expensive part. Pipeline per solve:
1. local divergence + Makhoul pack + **x-FFT** (x complete on-rank),
2. **all-to-all transpose**: z-slabs -> kx-slabs (each rank ends with all
   z and all y for a contiguous subset of x-modes),
3. local **z-FFT + Thomas in y** (both directions complete on-rank —
   this is why y and z must not be split simultaneously),
4. inverse z-FFT, **all-to-all back**, local x-IFFT + unpack.

Estimated overhead at P=4 on PCIe: ~3-4 ms/step on the transition grid
vs ~10 ms/step single-GPU compute -> expected strong scaling ~1.8x at 2
ranks, ~3-3.5x at 4 (PCIe-bound; on NVLink/NVSwitch hardware this design
scales near-linearly). Weak scaling (bigger grids per rank) fares
better since compute grows faster than the transpose.

## Communication layer

Start with **CUDA-aware OpenMPI** (ships with NVHPC; `mpirun` workflow
identical to legacy; device pointers passed straight to MPI_Alltoall /
MPI_Sendrecv). Move the transposes to **NCCL** later only if profiling
shows MPI overhead matters. cuFFTMp (NVSHMEM) is the turnkey
alternative for step 2-4 but imposes its own data layout and does not
know about our custom half-length DCT in x; the hand-rolled slab
transpose keeps the validated single-GPU kernels unchanged.

## Determinism is preserved

Per-rank reductions stay fixed-order (reductions.cuf unchanged); the
inter-rank combination gathers per-rank partials and combines them in
rank order on rank 0 (or identically on all ranks) — same input, same
result, any launch, any P... at fixed P. Note results at P=4 will differ
from P=1 in the last bits (different reduction split), exactly like the
legacy MPI code.

## Feature interactions

- **HIT inflow (flag 6)**: every rank holds inlet columns (z-decomp), so
  each rank loads/interpolates only its z-window of the plane library —
  the buffered reader already indexes planes independently.
- **Box output**: each rank writes its z-window; a post-hoc cat, or the
  writer gathers to rank 0 (boxes are small).
- **Statistics/monitor**: z-averages become allreduces (as in legacy).
- **Restart I/O**: rank 0 reads/writes with per-rank z-plane scatter
  (the legacy reader's structure, minus its int32 offset bugs — the
  fixed Int64 POS= addressing is already in both trees).

## Implementation phases (each harness-gated)

1. **P5.1 — infrastructure**: MPI init, per-rank grids/allocation, halo
   exchange module (device-direct), distributed RHS/advance/BCs; Poisson
   still on rank 0 via gather/scatter. Validate on 2 ranks vs single-GPU
   (fields must match to reduction-split roundoff, ~1e-12).
2. **P5.2 — distributed Poisson**: the transpose pipeline above.
   Remove the gather. Bench 2/4 ranks.
3. **P5.3 — distributed features**: HIT inflow z-windows, stats
   allreduces, box output, restart scatter/gather.
4. **P5.4 — scaling study**: strong + weak scaling table on tifa (2, 4
   ranks), added to validation/README.md; production-grid shakedown
   (the 270M-cell config on 4 ranks).

## Non-goals

Multi-node (tifa is one node); y-decomposition; overlapping
communication with the Thomas solve (revisit only if the transpose
dominates at P=4).
