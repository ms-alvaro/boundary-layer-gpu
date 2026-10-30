# Cases

Runnable configurations, smallest first. All use the `.turbb` input format,
which is `key = value` text with unknown keys ignored. The annotated
parameter list is [docs/input_parameters.turbb](../docs/input_parameters.turbb).

## laminar/: Blasius flat plate (the correctness gate)

302x64x8, 2000 steps, about 8 s on one A100. The inlet feeds the Blasius
similarity solution and the run is checked against it (Cf, profiles).

```bash
cd cases/laminar
python3 generate_blasius.py                     # similarity solution table
mpirun -np 1 ../../boundary_layer_cuda -i laminar.turbb
python3 postprocess.py                          # figures vs Blasius
```

The automated version of this check is `tests/validate.py laminar`.

## transition/: K-type transition from TS inflow modes

814x125x66, 500k steps (about 10 flow-throughs, about 1.5 h on one A100).
Tollmien-Schlichting eigenmodes from the Orr-Sommerfeld solver in
`generate_temporal_modes_local.py` enter at the inlet, grow, and break
down to turbulence. `plot_ic_inflow.py <config>` plots the initial
condition and the inflow before the run is launched. `plot_snapshot.py`
renders results.

## turbulent_lund/: turbulent boundary layer with Lund recycling inflow

The Lund recycling inflow (`inflow_flag = 3/5`) and the constant-coefficient
Smagorinsky model (`LES = 1`) of this case run in the CUDA solver and were
validated against the CPU solver to machine precision. The shipped input
restarts from an external flow field that is not in the repository, so the
case is kept as the reference configuration for the recycling-inflow path.
See [turbulent_lund/README.md](turbulent_lund/README.md) for the file
list, the validation table and the porting notes.

## profiles/

Mean boundary-layer profiles (`Re_theta.*.prof`) read by the inflow
generators (`generate_inflow*.py`) when they build mean-profile inflow
files.

## HIT free-stream turbulence inflow (inflow_flag = 6)

The bypass-transition configuration, a Blasius base flow plus time-resolved
turbulence planes from a precursor HIT simulation, needs a plane library
(format v2) produced by a preprocessing pipeline that is not part of this
repository. The solver-side reader, its consistency checks and the
no-recycling guard are in `src/hit_inflow_mod.cuf`. The input keys are
`hit_file`, `N_buffer_hit` and `ygrid_file`.
