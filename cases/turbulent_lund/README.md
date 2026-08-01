# Turbulent boundary layer — Lund–Wu–Squires recycling inflow

A zero-pressure-gradient turbulent boundary layer sustained by Lund-type
recycling: the solution at a downstream plane (`Lund_ix`) is rescaled to
inlet units every step and re-injected as the inflow condition
(`inflow_flag = 3` full rescaling, `= 5` fluctuation-only rescaling; see
[docs/input_parameters.turbb](../../docs/input_parameters.turbb)).

## Files

- `turbulent_lund.turbb` — reference input for this case (grid, Lund plane,
  rescaling parameters). Start from this file.
- `Mean_profile_Retheta1100_ny192_Ly0.60.dat` — target mean profile
  (Re_theta = 1100) used to seed the initial condition.
- `pyfiles/` — post-processing helpers (`BL.py`, `check.py`).

## Status in the CUDA solver

The Lund recycling inflow (`inflow_flag = 3/5`) is **not yet ported** to the
CUDA solver — requesting it stops the run with a clear message. Additionally
the case needs an external turbulent restart field to avoid a very long
laminar-to-turbulent development transient.

Port checklist (contributions welcome):
1. Device-side plane extraction at `Lund_ix` (exists for box output — reuse).
2. Inner/outer rescaling of mean + fluctuations (Lund, Wu & Squires 1998),
   inlet-unit conversion via `Lund_deltai`, `Lund_T` running average.
3. Injection through the existing inflow ghost-plane path
   (`hit_inflow_mod.cuf` shows the pattern used by `inflow_flag = 6`).
4. Validate against the pre-port MPI implementation (`gpu-openacc` branch)
   on this case: mean profile, Cf, and spectra at the recycling plane.
