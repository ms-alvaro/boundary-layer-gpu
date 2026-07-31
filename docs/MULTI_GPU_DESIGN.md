# Multi-GPU design (Phase 5 proposal — not yet implemented)

Status update 3: **P5.3 + P5.4 implemented and validated** (2026-07-31):
- **Statistics**: distributed via partial z-sums + array allreduces (the
  overlapping-slab interiors tile the global z exactly, so the reduced
  statistics are numerically IDENTICAL to single-rank — verified 0.00
  difference on the .stats.txt values).
- **Box output**: per-rank z-window extraction + float32 plane gather to
  rank 0; identical headers/grids, payloads within one float32 ulp.
- **Restart**: fully parallel windowed reads — every rank POS-seeks its own
  z-slab in the global snapshot (no MPI in the reader); the writer's
  gathered global snapshots make restart files interchangeable between any
  rank counts. A post-load z-ghost sync (the legacy 'force periodicity'
  fixup) makes single-rank restart continuation BITWISE-exact vs the
  straight-through run; multi-rank matches at reduction roundoff (7e-15).
- **Scaling (P5.4)**, transition grid 6.7M cells, PCIe A100s:
  P=1: 10.0, P=2: 15.0, P=4: 18.7 ms/step — strong scaling is NEGATIVE at
  this size on PCIe (as flagged below: transposes + syncs dominate small
  grids without NVLink). Multi-GPU on tifa is a CAPACITY feature: the
  270M-cell fst_tbl_ctr-scale configuration (~65 GB of state, impossible
  on one 40 GB card) runs on 4 ranks — see validation/README.md.
- Bug found by the 2-rank HIT validation: the plane-library z-period check
  used the LOCAL nz (half the period at P=2); fixed to nz_global.

Status update 2: **P5.2 implemented and validated** (2026-07-31): the
distributed transpose Poisson is in (no gather; every rank x-FFTs its own
z-slab, all-to-alls the REAL DCT coefficients — half the bytes — to x-mode
slabs, runs z-FFT + its Thomas-factor slice locally, and transposes back).
2-rank transition fields identical to single-rank (8.7e-15). The P5.1
gathered path remains as a fallback via BL_POISSON_GATHER=1.

Measured 2-rank timing on the transition grid (6.7M cells): 15 ms/step vs
10 single-rank — i.e. **no speedup at P = 2 on PCIe at this size**, and
UCX CUDA-IPC transports are REQUIRED to get even that (26.8 ms/step with
the default ob1 transports; launch with
`--mca pml ucx -x UCX_TLS=self,sm,cuda_copy,cuda_ipc -x UCX_MEMTYPE_CACHE=n`).
The honest position: on tifa's PCIe-only A100s, multi-GPU buys CAPACITY
(grids beyond one card's 40 GB, e.g. the 270M-cell fst_tbl_ctr config on
4 ranks) rather than turnaround at test sizes; the compute does scale
(RHS halves per rank), so larger per-rank grids and NVLink-class hardware
both improve the ratio. Next levers if speed at small P matters: NCCL
all-to-all, packed halo messages, comm/compute overlap.

Status: **P5.1 implemented and validated** (2026-07-31): z-slab
decomposition, device-direct halo/periodic exchanges, distributed
RHS/advance/BCs and mass conservation, Poisson via rank-0 gather at global
size, gathered snapshot writer. 2-rank transition fields match single-rank
at 9e-15 after 2000 steps. As expected for the gather-Poisson interim,
2 ranks are ~2x SLOWER than 1 (the 4 gather/scatters per step dominate at
test-grid size) — P5.2 (distributed transposes) is the speed step. P5.1
restrictions: restart, statistics, and boxout are single-rank only.

Implementation lesson recorded: CUDA-aware MPI is not stream-aware — a
cudaDeviceSynchronize is required before handing device pointers to MPI
(kernels still writing the buffer otherwise race the transfer; symptom was
a timing-dependent blow-up that the z-uniform laminar case could not see).

The remaining phases below are unchanged.

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
