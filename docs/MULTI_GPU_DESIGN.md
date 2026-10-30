# Multi-GPU design

This document holds the Phase 5 proposal for multi-GPU execution and the
status of its implementation. The status comes first. The design sections
that follow are the original proposal and are kept as written.

## Status

P5.1 was implemented and validated on 2026-07-31: z-slab decomposition,
device-direct halo and periodic exchanges, distributed RHS, advance, BCs
and mass conservation, a Poisson solve via rank-0 gather at global size,
and a gathered snapshot writer. The 2-rank transition fields match the
single-rank fields at 9e-15 after 2000 steps. As expected for the
gather-Poisson interim, 2 ranks were about 2x slower than 1, because the 4
gathers and scatters per step dominate at test-grid size. P5.2
(distributed transposes) was the speed step. Under P5.1, restart,
statistics and box output were single-rank only.

One implementation lesson is recorded here. CUDA-aware MPI is not
stream-aware, so a cudaDeviceSynchronize is required before handing device
pointers to MPI. Otherwise kernels still writing the buffer race the
transfer. The symptom was a timing-dependent blow-up that the z-uniform
laminar case could not see.

P5.2 was implemented and validated on 2026-07-31. The distributed
transpose Poisson is in, with no gather: every rank x-FFTs its own z-slab,
all-to-alls the real DCT coefficients (half the bytes) to x-mode slabs,
runs the z-FFT and its slice of the Thomas factors locally, and transposes
back. The 2-rank transition fields are identical to single-rank (8.7e-15).
The P5.1 gathered path remains as a fallback via `BL_POISSON_GATHER=1`.

The measured 2-rank timing on the transition grid (6.7M cells) is 15
ms/step against 10 single-rank, so there is no speedup at P = 2 on PCIe at
this size. The UCX CUDA-IPC transports are required to get even that
(26.8 ms/step with the default ob1 transports). Launch with
`--mca pml ucx -x UCX_TLS=self,sm,cuda_copy,cuda_ipc -x UCX_MEMTYPE_CACHE=n`.
On tifa's PCIe-only A100s, multi-GPU buys capacity (grids beyond one
card's 40 GB, for example the 270M-cell fst_tbl_ctr configuration on 4
ranks) rather than turnaround at test sizes. The compute does scale (the
RHS halves per rank), so larger per-rank grids and NVLink-class hardware
both improve the ratio. The next levers if speed at small P matters are
NCCL all-to-all, packed halo messages, and overlap of communication with
compute.

P5.3 and P5.4 were implemented and validated on 2026-07-31:

- Statistics are distributed via partial z-sums and array allreduces. The
  overlapping-slab interiors tile the global z exactly, so the reduced
  statistics are numerically identical to single-rank. The difference on
  the `.stats.txt` values is 0.00.
- Box output: per-rank z-window extraction and a float32 plane gather to
  rank 0. The headers and grids are identical and the payloads are within
  one float32 ulp.
- Restart: fully parallel windowed reads. Every rank POS-seeks its own
  z-slab in the global snapshot, with no MPI in the reader. The writer's
  gathered global snapshots make restart files interchangeable between any
  rank counts. A post-load z-ghost sync (the legacy "force periodicity"
  fixup) makes single-rank restart continuation bitwise-exact against the
  straight-through run. Multi-rank matches at reduction roundoff (7e-15).
- Scaling (P5.4), transition grid of 6.7M cells, PCIe A100s: P=1: 10.0,
  P=2: 15.0, P=4: 18.7 ms/step. Strong scaling is negative at this size on
  PCIe, as flagged in the design below: transposes and syncs dominate
  small grids without NVLink. Multi-GPU on tifa is a capacity feature.
  The 270M-cell fst_tbl_ctr-scale configuration (about 65 GB of state,
  impossible on one 40 GB card) runs on 4 ranks. See
  [../tests/README.md](../tests/README.md).
- A bug found by the 2-rank HIT validation: the plane-library z-period
  check used the local nz (half the period at P=2). It now uses
  nz_global.

The sections below are the original proposal.

## Why

1. Memory: production-scale grids outgrow one A100-40GB. The solver
   holds about 30 full-size double arrays. A 3074x341x258 grid (270M
   cells, the fst_tbl_ctr scale) needs about 65 GB of state, which is
   impossible on one 40 GB card and comfortable on 4 (about 16 GB/rank).
2. Turnaround: 5 A100s on tifa sit mostly idle while one runs.

## Hardware reality (tifa)

5x A100-PCIE-40GB with no NVLink. Peer-to-peer goes over PCIe Gen4
(about 20-25 GB/s effective per direction, shared lanes). Communication is
therefore the scarce resource. The design minimizes global data motion
and overlaps what remains.

## Decomposition: z-slabs

Split the spanwise (periodic) direction across ranks, as the legacy MPI
code did. The constraint is `mod(nz-2, nprocs) = 0`, so nz = 258 means 4
ranks, not 5.

Why z and not x or y:
- z is periodic: the halo logic is uniform, with no special inflow or
  outflow ranks.
- y hosts the Thomas solve (a serial recurrence). Splitting y would
  serialize across ranks.
- x works too, and balances better for long domains, but the transform
  bookkeeping matches the existing single-GPU code best with x kept
  local (the half-length DCT stays entirely on-rank).

Per-rank state: U(nx, nyg, nzl+2) etc., with nzl = (nz-2)/P local planes
plus one halo plane per side.

## What communicates, and how much

Per RK2 substep (2 substeps/step):

| operation | pattern | volume/rank (transition grid, P=4) |
|---|---|---|
| RHS halos (U,V,W, 1 plane/side) | neighbor sendrecv | ~0.8 MB x 6 -> ~5 MB, <1 ms |
| mass-conservation flux | 4-scalar allreduce | negligible (latency only) |
| Poisson forward transpose | all-to-all | ~1/P of the 103 MB mode array |
| Poisson inverse transpose | all-to-all | same |

The Poisson solve is the expensive part. The pipeline per solve is:
1. local divergence + Makhoul pack + x-FFT (x complete on-rank),
2. all-to-all transpose from z-slabs to kx-slabs (each rank ends with all
   z and all y for a contiguous subset of x-modes),
3. local z-FFT + Thomas in y (both directions complete on-rank, which is
   why y and z must not be split simultaneously),
4. inverse z-FFT, all-to-all back, local x-IFFT + unpack.

The estimated overhead at P=4 on PCIe is about 3-4 ms/step on the
transition grid against about 10 ms/step of single-GPU compute. The
expected strong scaling is therefore about 1.8x at 2 ranks and 3-3.5x at
4. This is PCIe-bound. On NVLink or NVSwitch hardware this design scales
near-linearly. Weak scaling (bigger grids per rank) fares better, since
compute grows faster than the transpose.

## Communication layer

Start with CUDA-aware OpenMPI, which ships with NVHPC. The `mpirun`
workflow is identical to legacy, and device pointers are passed straight
to MPI_Alltoall and MPI_Sendrecv. Move the transposes to NCCL later only if
profiling shows that MPI overhead matters. cuFFTMp (NVSHMEM) is the
turnkey alternative for steps 2-4, but it imposes its own data layout and
does not know about our custom half-length DCT in x. The hand-rolled slab
transpose keeps the validated single-GPU kernels unchanged.

## Determinism is preserved

Per-rank reductions stay fixed-order (reductions.cuf unchanged). The
inter-rank combination gathers per-rank partials and combines them in
rank order on rank 0 (or identically on all ranks), so the same input
gives the same result on any launch at fixed P. Results at P=4 differ
from P=1 in the last bits, because the reduction split differs, exactly
like the legacy MPI code.

## Feature interactions

- HIT inflow (flag 6): every rank holds inlet columns (z-decomposition),
  so each rank loads and interpolates only its z-window of the plane
  library. The buffered reader already indexes planes independently.
- Box output: each rank writes its z-window. Either a post-hoc cat, or the
  writer gathers to rank 0 (boxes are small).
- Statistics and monitor: z-averages become allreduces (as in legacy).
- Restart I/O: rank 0 reads and writes with a per-rank z-plane scatter.
  This is the legacy reader's structure, minus its int32 offset bugs. The
  fixed Int64 POS= addressing is already in both trees.

## Implementation phases (each harness-gated)

1. P5.1, infrastructure: MPI init, per-rank grids and allocation, a halo
   exchange module (device-direct), distributed RHS, advance and BCs.
   Poisson still on rank 0 via gather and scatter. Validate on 2 ranks
   against single-GPU (fields must match to reduction-split roundoff,
   about 1e-12).
2. P5.2, distributed Poisson: the transpose pipeline above. Remove the
   gather. Bench 2 and 4 ranks.
3. P5.3, distributed features: HIT inflow z-windows, statistics
   allreduces, box output, restart scatter and gather.
4. P5.4, scaling study: strong and weak scaling table on tifa (2, 4
   ranks), added to [../tests/README.md](../tests/README.md), and a
   production-grid shakedown (the 270M-cell config on 4 ranks).

## Non-goals

Multi-node (tifa is one node), y-decomposition, and overlapping
communication with the Thomas solve (revisit only if the transpose
dominates at P=4).
