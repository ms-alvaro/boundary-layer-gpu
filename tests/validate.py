#!/usr/bin/env python3
"""
Validation and benchmarking harness for the boundary-layer solver.

Every solver change (OpenACC reference -> CUDA port -> optimizations) is gated
on this harness: same physics, measured speedup.

Subcommands
-----------
  laminar  Run the Blasius laminar case with a given executable, check the
           physics (Cf vs Blasius, divergence, freestream overshoot), plot the
           inflow/inlet profiles, and write a machine-readable summary.json.
  bench    Run a production-sized DNS grid for a fixed number of steps and
           report ns/cell/step (timing only, no physics checks).
  compare  Compare two laminar runs field-by-field (final snapshot) and report
           timing side by side. Used to verify a new build against a reference.

Examples
--------
  python3 tests/validate.py laminar --exe ./boundary_layer_gpu --label openacc-ref
  python3 tests/validate.py bench   --exe ./boundary_layer_gpu --label openacc-ref
  python3 tests/validate.py compare --a openacc-ref --b cuda-v1 --tol 1e-8

Runs live in tests/runs/<label>/{laminar,bench}/ and are disposable;
summary.json in each run dir is the durable record.
"""
import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(HERE, 'runs')

# Golden values from the validated OpenACC baseline (A100, nvhpc 24.3,
# cases/laminar/laminar.turbb, 2000 steps). Physics tolerances are loose enough
# to accept floating-point reordering from kernel rewrites, tight enough to
# catch real bugs.
GOLDEN = {
    'mean_cf': 9.7410751327e-4,   # monitor "Mean Cf" at step 2000
    'mean_cf_rtol': 1e-4,
    'max_div': 1e-12,             # machine-zero divergence after projection
    'max_u_overshoot': 1e-4,      # |max(U) - 1| in the freestream
    'cf_blasius_rtol': 0.03,      # pointwise Cf vs 0.664/sqrt(Re_x), buffers excluded
}


# ----------------------------------------------------------------------------
# File readers (formats defined by input_output.f90 / statistics.f90)
# ----------------------------------------------------------------------------
def read_stats(fname):
    """Read a .stats.txt file: header, Cf(x), x, y, Umean(x,y), Vmean(x,y)."""
    with open(fname) as f:
        h = f.readline().strip().lstrip('%').split()
        t, nu = float(h[0]), float(h[1])
        nx, ny = int(h[2]), int(h[3])
        cf = [float(v) for v in f.readline().split()]
        x = [float(v) for v in f.readline().split()]
        y = [float(v) for v in f.readline().split()]
        umean = [[float(v) for v in f.readline().split()] for _ in range(nx)]
        vmean = [[float(v) for v in f.readline().split()] for _ in range(nx)]
    return {'t': t, 'nu': nu, 'nx': nx, 'ny': ny,
            'Cf': cf, 'x': x, 'y': y, 'Umean': umean, 'Vmean': vmean}


def read_snapshot(fname):
    """Read a binary snapshot: grids + U, V, W fields (see input_output.f90)."""
    import numpy as np
    with open(fname, 'rb') as f:
        t, nu = struct.unpack('dd', f.read(16))
        marker = struct.unpack('i', f.read(4))[0]
        if marker == -73:
            _istep, nx = struct.unpack('ii', f.read(8))
        else:
            nx = marker
        x = np.frombuffer(f.read(8 * nx), dtype=np.float64)
        ny = struct.unpack('i', f.read(4))[0]
        y = np.frombuffer(f.read(8 * ny), dtype=np.float64)
        nz = struct.unpack('i', f.read(4))[0]
        z = np.frombuffer(f.read(8 * nz), dtype=np.float64)
        grids = {}
        for name in ('xm', 'ym', 'zm'):
            n = struct.unpack('i', f.read(4))[0]
            grids[name] = np.frombuffer(f.read(8 * n), dtype=np.float64)

        def field():
            n1, n2, n3 = struct.unpack('iii', f.read(12))
            n3 -= 1  # code writes planes 1..n3-1
            a = np.frombuffer(f.read(8 * n1 * n2 * n3), dtype=np.float64)
            return a.reshape((n3, n2, n1)).transpose(2, 1, 0)

        U, V, W = field(), field(), field()
    return {'t': t, 'nu': nu, 'x': x, 'y': y, 'z': z, **grids,
            'U': U, 'V': V, 'W': W}


def parse_monitor_log(logfile):
    """Extract per-interval timings and final monitor values from stdout."""
    text = open(logfile, errors='replace').read()
    grab = lambda pat: [float(v) for v in re.findall(pat, text)]
    out = {
        'elapsed_intervals': grab(r'Elapsed time \(s\)\s*:\s*([0-9.Ee+-]+)'),
        'mean_cf': grab(r'Mean Cf:\s*:\s*([0-9.Ee+-]+)'),
        'max_div': grab(r'Maximum divergence\s*:\s*([0-9.Ee+-]+)'),
        'max_u': grab(r'Maximum U\s*:\s*([0-9.Ee+-]+)'),
    }
    prof = re.findall(r'PROFILE n=\s*\d+ \(substep1 avg ms\): copy=\s*([0-9.]+)'
                      r' RHS=\s*([0-9.]+) adv=\s*([0-9.]+) BC=\s*([0-9.]+)'
                      r' proj=\s*([0-9.]+)', text)
    if prof:
        last = prof[-1]
        out['profile_substep1_ms'] = dict(
            zip(('copy', 'RHS', 'adv', 'BC', 'proj'), map(float, last)))
    return out


# ----------------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------------
def pick_free_gpu():
    """Return the index of a fully free GPU (house rule: 0 MiB-ish, 0%)."""
    q = subprocess.run(
        ['nvidia-smi', '--query-gpu=index,memory.used,utilization.gpu',
         '--format=csv,noheader,nounits'], capture_output=True, text=True)
    for line in q.stdout.strip().splitlines():
        idx, mem, util = [int(v) for v in line.split(',')]
        if mem < 20 and util == 0:
            return idx
    sys.exit('ERROR: no fully free GPU available (house rule: wait for one).')


def run_case(exe, case_files, turbb, rundir, nsteps=None, gpu=None):
    """Set up rundir with the case inputs and run the solver in it."""
    exe = os.path.abspath(exe)
    os.makedirs(rundir, exist_ok=True)
    for f in case_files:
        shutil.copy(f, rundir)
    if nsteps is not None:  # override step count for benchmarking
        src = os.path.join(rundir, os.path.basename(turbb))
        txt = re.sub(r'nsteps\s*=\s*\d+', f'nsteps   = {nsteps}',
                     open(src).read())
        open(src, 'w').write(txt)
    if os.path.isdir(os.path.join(rundir, 'data')):
        shutil.rmtree(os.path.join(rundir, 'data'))

    gpu = pick_free_gpu() if gpu is None else gpu
    env = dict(os.environ)
    nvhpc = env.get('NVHPC_ROOT', '/opt/nvidia/hpc_sdk/Linux_x86_64/24.3')
    env['PATH'] = (f'{nvhpc}/compilers/bin:{nvhpc}/comm_libs/mpi/bin:'
                   + env['PATH'])
    env['CUDA_VISIBLE_DEVICES'] = str(gpu)

    logfile = os.path.join(rundir, 'run.log')
    print(f'  running {os.path.basename(exe)} on GPU {gpu} in {rundir}')
    with open(logfile, 'w') as log:
        p = subprocess.run(
            ['mpirun', '-np', '1', exe, '-i', os.path.basename(turbb)],
            cwd=rundir, env=env, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        sys.exit(f'ERROR: solver exited with code {p.returncode}, '
                 f'see {logfile}')
    return logfile


# ----------------------------------------------------------------------------
# Physics checks (laminar)
# ----------------------------------------------------------------------------
def check_laminar(rundir):
    """Return (checks dict, all_passed) for a completed laminar run."""
    import numpy as np
    mon = parse_monitor_log(os.path.join(rundir, 'run.log'))
    stats = read_stats(os.path.join(rundir, 'data',
                                    'BL_laminar.00002000.stats.txt'))
    checks, results = {}, []

    def check(name, value, ok, detail=''):
        checks[name] = {'value': value, 'pass': bool(ok), 'detail': detail}
        results.append(ok)
        print(f'  [{"PASS" if ok else "FAIL"}] {name}: {value:.6e} {detail}')

    cf = mon['mean_cf'][-1]
    check('mean_cf_vs_golden', cf,
          abs(cf - GOLDEN['mean_cf']) / GOLDEN['mean_cf'] < GOLDEN['mean_cf_rtol'],
          f'(golden {GOLDEN["mean_cf"]:.6e}, rtol {GOLDEN["mean_cf_rtol"]})')

    div = max(mon['max_div'])
    check('max_divergence', div, div < GOLDEN['max_div'],
          f'(< {GOLDEN["max_div"]})')

    over = max(abs(u - 1.0) for u in mon['max_u'])
    check('freestream_overshoot', over, over < GOLDEN['max_u_overshoot'],
          f'(< {GOLDEN["max_u_overshoot"]})')

    # Pointwise Cf vs Blasius 0.664/sqrt(Re_x), inlet/outlet buffers excluded
    x = np.array(stats['x'])
    cf_sim = np.array(stats['Cf'])
    s = slice(3, -5)
    cf_bl = 0.664 / np.sqrt(x[s] / stats['nu'])
    err = float(np.max(np.abs(cf_sim[s] - cf_bl) / cf_bl))
    check('cf_vs_blasius_maxrelerr', err, err < GOLDEN['cf_blasius_rtol'],
          f'(< {GOLDEN["cf_blasius_rtol"]})')

    return checks, all(results), mon, stats


def plot_ic_inflow(rundir, stats):
    """Standing rule: plot the inflow BC (and evolved inlet state) for every
    BL run. Inlet U,V profiles vs the analytical Blasius solution."""
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    sol = solve_ivp(lambda e, s: [s[1], s[2], -s[0] * s[2]],
                    (0, 20), [0, 0, 0.46960],
                    t_eval=np.linspace(0, 20, 2000), rtol=1e-12, atol=1e-12)
    eta_bl, f_bl, df_bl = sol.t, sol.y[0], sol.y[1]

    nu = stats['nu']
    x0 = stats['x'][0]
    y = np.array(stats['y'])
    eta = y * np.sqrt(1.0 / (2 * nu * x0))
    U_in = np.array(stats['Umean'][0])
    V_in = np.array(stats['Vmean'][0])
    # Blasius V = sqrt(nu U / 2x) * (eta f' - f)
    V_bl = np.sqrt(nu / (2 * x0)) * (eta_bl * df_bl - f_bl)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].plot(U_in, eta, 'bo', ms=3, label='inlet (sim)')
    axes[0].plot(df_bl, eta_bl, 'k-', lw=1, label='Blasius')
    axes[0].set(xlabel='U', ylabel=r'$\eta$', ylim=(0, 8),
                title=f'Inflow BC: U at x={x0:.2f}')
    axes[1].plot(V_in, eta, 'ro', ms=3, label='inlet (sim)')
    axes[1].plot(V_bl, eta_bl, 'k-', lw=1, label='Blasius')
    axes[1].set(xlabel='V', ylabel=r'$\eta$', ylim=(0, 8),
                title='Inflow BC: V')
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(rundir, 'ic_inflow.png')
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f'  inflow figure: {out}')


# ----------------------------------------------------------------------------
# Timing
# ----------------------------------------------------------------------------
def timing_summary(mon, ncells, steps_per_interval):
    """ns/cell/step from the per-interval monitor clock (first interval
    dropped: it absorbs warm-up)."""
    iv = mon['elapsed_intervals']
    if len(iv) > 2:
        iv = iv[1:]
    per_step = sorted(iv)[len(iv) // 2] / steps_per_interval  # median
    return {'median_s_per_step': per_step,
            'ns_per_cell_per_step': per_step / ncells * 1e9,
            'intervals_s': iv}


# ----------------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------------
def cmd_laminar(args):
    case = os.path.join(REPO, 'cases', 'laminar')
    rundir = os.path.join(RUNS, args.label, 'laminar')
    files = [os.path.join(case, f)
             for f in ('laminar.turbb', 'blasius_solution.dat')]
    run_case(args.exe, files, 'laminar.turbb', rundir, gpu=args.gpu)

    checks, ok, mon, stats = check_laminar(rundir)
    ncells = 302 * 64 * 8
    timing = timing_summary(mon, ncells, steps_per_interval=100)
    print(f'  timing: {timing["median_s_per_step"]*1e3:.3f} ms/step '
          f'({timing["ns_per_cell_per_step"]:.2f} ns/cell/step)')
    plot_ic_inflow(rundir, stats)

    summary = {'label': args.label, 'exe': os.path.abspath(args.exe),
               'case': 'laminar', 'passed': ok, 'checks': checks,
               'timing': timing,
               'profile_substep1_ms': mon.get('profile_substep1_ms')}
    with open(os.path.join(rundir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(f'  {"ALL CHECKS PASSED" if ok else "*** CHECKS FAILED ***"} '
          f'-> {rundir}/summary.json')
    sys.exit(0 if ok else 1)


def cmd_bench(args):
    case = os.path.join(HERE, 'cases')
    rundir = os.path.join(RUNS, args.label, 'bench')
    files = [os.path.join(case, 'bench_dns.turbb'),
             os.path.join(REPO, 'cases', 'laminar', 'blasius_solution.dat')]
    run_case(args.exe, files, 'bench_dns.turbb', rundir,
             nsteps=args.nsteps, gpu=args.gpu)
    mon = parse_monitor_log(os.path.join(rundir, 'run.log'))
    ncells = 814 * 125 * 66
    timing = timing_summary(mon, ncells, steps_per_interval=50)
    summary = {'label': args.label, 'exe': os.path.abspath(args.exe),
               'case': 'bench_dns (814x125x66)', 'timing': timing,
               'profile_substep1_ms': mon.get('profile_substep1_ms'),
               'max_div_final': mon['max_div'][-1] if mon['max_div'] else None}
    with open(os.path.join(rundir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(f'  bench: {timing["median_s_per_step"]*1e3:.2f} ms/step, '
          f'{timing["ns_per_cell_per_step"]:.3f} ns/cell/step '
          f'-> {rundir}/summary.json')


def cmd_compare(args):
    import numpy as np
    a_dir = os.path.join(RUNS, args.a, 'laminar')
    b_dir = os.path.join(RUNS, args.b, 'laminar')
    snap = 'data/BL_laminar.00002000'
    A = read_snapshot(os.path.join(a_dir, snap))
    B = read_snapshot(os.path.join(b_dir, snap))

    print(f'  field-level diffs ({args.a} vs {args.b}), final snapshot:')
    ok = True
    diffs = {}
    for name in ('U', 'V', 'W'):
        d = np.abs(A[name] - B[name])
        scale = max(np.abs(A[name]).max(), 1e-300)
        rel = d.max() / scale
        diffs[name] = {'max_abs': float(d.max()), 'max_rel': float(rel),
                       'rms': float(np.sqrt((d ** 2).mean()))}
        passed = bool(rel < args.tol)
        ok = ok and passed
        print(f'  [{"PASS" if passed else "FAIL"}] {name}: '
              f'max|diff|={d.max():.3e}  rel={rel:.3e}  (tol {args.tol})')

    for label, d in ((args.a, a_dir), (args.b, b_dir)):
        s = json.load(open(os.path.join(d, 'summary.json')))
        t = s['timing']
        print(f'  {label:>20}: {t["median_s_per_step"]*1e3:.3f} ms/step '
              f'({t["ns_per_cell_per_step"]:.2f} ns/cell/step)')
    sa = json.load(open(os.path.join(a_dir, 'summary.json')))
    sb = json.load(open(os.path.join(b_dir, 'summary.json')))
    speedup = (sa['timing']['median_s_per_step']
               / sb['timing']['median_s_per_step'])
    print(f'  speedup ({args.b} vs {args.a}): {speedup:.2f}x')

    out = {'a': args.a, 'b': args.b, 'tol': args.tol, 'passed': ok,
           'field_diffs': diffs, 'speedup_b_over_a': speedup}
    outfile = os.path.join(RUNS, f'compare_{args.a}_vs_{args.b}.json')
    json.dump(out, open(outfile, 'w'), indent=2)
    print(f'  {"FIELDS MATCH" if ok else "*** FIELDS DIFFER ***"} -> {outfile}')
    sys.exit(0 if ok else 1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    pl = sub.add_parser('laminar', help='run + validate the Blasius case')
    pl.add_argument('--exe', required=True)
    pl.add_argument('--label', required=True,
                    help='run name, e.g. openacc-ref, cuda-v1')
    pl.add_argument('--gpu', type=int, default=None)
    pl.set_defaults(fn=cmd_laminar)

    pb = sub.add_parser('bench', help='production-sized timing run')
    pb.add_argument('--exe', required=True)
    pb.add_argument('--label', required=True)
    pb.add_argument('--nsteps', type=int, default=500)
    pb.add_argument('--gpu', type=int, default=None)
    pb.set_defaults(fn=cmd_bench)

    pc = sub.add_parser('compare', help='diff two laminar runs + speedup')
    pc.add_argument('--a', required=True, help='reference run label')
    pc.add_argument('--b', required=True, help='candidate run label')
    pc.add_argument('--tol', type=float, default=1e-8,
                    help='max relative field difference to pass')
    pc.set_defaults(fn=cmd_compare)

    args = p.parse_args()
    args.fn(args)


if __name__ == '__main__':
    main()
