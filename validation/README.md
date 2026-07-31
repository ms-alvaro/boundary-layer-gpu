# Validation and benchmarking harness

Every solver change is gated on this harness: **same physics, measured
speedup**. It runs the solver, checks the results numerically (no eyeballing),
and records timing in a machine-readable form.

## Quick start

```bash
# from the repo root
python3 validation/validate.py laminar --exe ./boundary_layer_gpu --label openacc-ref
python3 validation/validate.py bench   --exe ./boundary_layer_gpu --label openacc-ref
python3 validation/validate.py compare --a openacc-ref --b cuda-v1
```

## What each subcommand does

- **laminar** — runs `laminar_test/laminar.turbb` (302x64x8 Blasius BL,
  2000 steps, ~8 s on an A100) and checks:
  - mean Cf vs the golden value from the validated OpenACC baseline
  - max divergence < 1e-12 (projection works)
  - freestream overshoot |max(U)-1| < 1e-4
  - pointwise Cf(x) vs Blasius 0.664/sqrt(Re_x) within 3% (buffers excluded)

  It also plots the inlet U,V profiles against the analytical Blasius
  solution (`ic_inflow.png` — inflow BC check) and writes `summary.json`.

- **bench** — runs a production-sized grid (814x125x66, 6.7M cells, the
  transition-case grid with plain Blasius inflow) for `--nsteps` steps and
  reports ns/cell/step. The laminar grid is launch-latency dominated, so
  speedups quoted for real runs must come from this case.

- **compare** — field-by-field diff (U,V,W of the final laminar snapshot)
  between two labeled runs plus a timing table and speedup. Default
  tolerance: max relative diff < 1e-8. The solver is deterministic
  (verified: two identical runs are bitwise equal), so any diff is caused
  by the code change being tested.

## Layout

```
validation/
  validate.py        the harness (single file, no exotic deps: numpy/scipy/mpl)
  cases/bench_dns.turbb
  runs/<label>/...   run dirs + summary.json (disposable, gitignored)
  runs/compare_*.json
```

## Baseline record (A100-PCIE-40GB, nvhpc 24.3, CUDA driver on tifa)

| build | case | ms/step | ns/cell/step |
|---|---|---|---|
| openacc-ref (commit 0537a6c) | laminar 302x64x8 | 2.35 | 15.2 |
| openacc-ref (commit 0537a6c) | bench 814x125x66 | 15.53 | 2.31 |

(Laminar-grid timing is kernel-launch-latency dominated; the bench number is
the one that predicts production performance.)
