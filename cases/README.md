# Cases

Runnable configurations, smallest first. All use the `.turbb` input format
(`key = value` text; unknown keys ignored; see `legacy/input_parameters.turbb`
for the annotated parameter list).

## laminar/ — Blasius flat plate (the correctness gate)

302x64x8, 2000 steps, ~8 s on one A100. The inlet feeds the Blasius
similarity solution; the run is checked against it (Cf, profiles).

```bash
cd cases/laminar
python3 generate_blasius.py                     # similarity solution table
mpirun -np 1 ../../boundary_layer_cuda -i laminar.turbb
python3 postprocess.py                          # figures vs Blasius
```

The automated version of this check is `validation/validate.py laminar`.

## transition/ — K-type transition from TS inflow modes

814x125x66, 500k steps (~10 flow-throughs, ~1.5 h on one A100).
Tollmien-Schlichting eigenmodes (from the Orr-Sommerfeld solver in
`generate_temporal_modes_local.py`) enter at the inlet, grow, and break
down to turbulence. `plot_ic_inflow.py <config>` previews the IC and
inflow before you commit the compute; `plot_snapshot.py` renders results.

## turbulent_lund/ — turbulent BL with Lund recycling inflow

**Legacy solver only** (Lund inflow and LES are not in the CUDA port), and
the config restarts from an external flow field that is not shipped in the
repo. Kept as the reference configuration for the recycling-inflow path.

## profiles/

Mean boundary-layer profiles (`Re_theta.*.prof`) consumed by the inflow
generators (`generate_inflow*.py`) when building mean-profile-based inflow
files.

## HIT freestream-turbulence inflow (inflow_flag = 6)

The bypass-transition configuration (Blasius base flow + time-resolved
turbulence planes from a precursor HIT simulation) needs a plane library
(format v2) produced by the accompanying preprocessing pipeline, which is
not part of this repository. The solver-side reader, its consistency
checks, and the no-recycling guard are in `src/hit_inflow_mod.cuf`; input
keys: `hit_file`, `N_buffer_hit`, `ygrid_file`.
