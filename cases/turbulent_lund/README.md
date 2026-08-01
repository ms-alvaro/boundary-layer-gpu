# Turbulent boundary layer — Lund–Wu–Squires recycling inflow

A zero-pressure-gradient turbulent boundary layer sustained by Lund-type
recycling: the solution at a downstream plane (`Lund_ix`) is rescaled to
inlet units every step and re-injected as the inflow condition
(`inflow_flag = 3` full rescaling, `= 5` fluctuation-only rescaling on
top of a prescribed mean profile; see
[docs/input_parameters.turbb](../../docs/input_parameters.turbb)).

## Files

- `turbulent_lund.turbb` — reference input for this case (grid, Lund plane,
  rescaling parameters). Start from this file.
- `Mean_profile_Retheta1100_ny192_Ly0.60.dat` — turbulent mean inlet
  profile (Re_theta = 1100) used by `inflow_flag = 5`. **Big-endian**
  binary (Intel-era, written under `F_UFMTENDIAN=big`); the CUDA solver
  reads it with `convert='big_endian'` automatically.
- `pyfiles/` — post-processing helpers (`BL.py`, `check.py`).

## Status in the CUDA solver: PORTED and validated

The Lund recycling inflow (`inflow_flag = 3/5`) and the blowing/suction
lid (`top_flag = 1/2`, Coleman 2018 / Abe 2017) are implemented in the
CUDA solver (`src/lund_inflow_mod.f90` + `bc_kernels.cuf`). The
rescaling math runs on the host exactly as the reference wrote it (the
recycling planes are tiny); only the two y-z planes move between device
and host each step.

Validated against the pre-port CPU MPI solver, both codes restarted from
the same noisy snapshot and advanced 200 steps (302x64x8 grid; max
absolute field differences):

| configuration | max U diff | max V diff |
|---|---|---|
| `inflow_flag = 3`, `top_flag = 0` | 3.6e-14 | 7.6e-15 |
| `inflow_flag = 5`, `top_flag = 0` | 2.5e-13 | 1.9e-14 |
| `inflow_flag = 3`, `top_flag = 1` | 4.4e-15 | 9.5e-16 |
| `inflow_flag = 5`, `top_flag = 1` | 4.1e-14 | 5.5e-15 |

Porting notes (also in the source comments):

- The reference's `zero_wz_top` re-imposition writes `V(i,nyg,:)` — one
  row past V's y extent. The out-of-bounds write lands on the wall row
  of a neighboring z-plane and is immediately overwritten by the wall
  BC, so its *effective* semantics (V keeps the blowing/suction profile
  over the whole lid; only the U and W top ghosts are re-imposed up to
  `Lund_ix`) are what the port reproduces — without the memory
  corruption.
- The `inflow_flag = 5` spanwise index swap for W
  (`mod(k+16,nz)+1`) uses the local `nz`, so the reference's result
  depends on the MPI rank count; the port matches the single-rank
  semantics.
- For `inflow_flag = 3` restarts, the time-averaged rescaling means are
  carried in a `<snapshot>.mean.rescaling` companion file (written
  native-endian by both this port and the gfortran/nvfortran CPU
  builds); if absent, the averages restart from the initial field's
  spanwise means, as in the reference.

## Running this case

The reference input expects a developed turbulent restart field
(`filein`) — a cold start from the mean profile alone would need a very
long laminar-to-turbulent transient. Two practical routes:

1. restart from an existing production snapshot of compatible grid, or
2. bootstrap: run the transition case to a turbulent state, interpolate
   onto this grid, and hand the result to `filein`.

The original study ran this case with `LES = 1` (constant-coefficient
Smagorinsky). LES is **not** ported — the CUDA solver covers the DNS
path only, so the reference input here sets `LES = 0`; on this grid
(1156x192x34) that makes it a coarse DNS, adequate for exercising the
recycling machinery, not for publication-grade statistics.
