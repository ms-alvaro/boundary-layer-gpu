"""
GPU transition monitoring: plots inflow, BL parameters, and snapshots.
Run periodically to check simulation progress.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import struct
import os
import glob

figdir = 'figures'
os.makedirs(figdir, exist_ok=True)

# ============================================================
# Readers
# ============================================================
def read_stats(fname):
    with open(fname, 'r') as f:
        header = f.readline().strip().lstrip('%').split()
        t, nu = float(header[0]), float(header[1])
        nx, ny, nz = int(header[2]), int(header[3]), int(header[4])
        istep = int(header[5])
        Cf = np.array(f.readline().split(), dtype=float)
        x = np.array(f.readline().split(), dtype=float)
        y = np.array(f.readline().split(), dtype=float)
        Umean = np.zeros((nx, ny))
        for i in range(nx):
            Umean[i, :] = np.array(f.readline().split(), dtype=float)
        Vmean = np.zeros((nx, ny))
        for i in range(nx):
            Vmean[i, :] = np.array(f.readline().split(), dtype=float)
    return {'t': t, 'nu': nu, 'nx': nx, 'ny': ny, 'nz': nz, 'istep': istep,
            'Cf': Cf, 'x': x, 'y': y, 'Umean': Umean, 'Vmean': Vmean}

def read_snapshot(fname):
    with open(fname, 'rb') as f:
        t, nu = struct.unpack('dd', f.read(16))
        marker = struct.unpack('i', f.read(4))[0]
        if marker == -73:
            istep = struct.unpack('i', f.read(4))[0]
            nx = struct.unpack('i', f.read(4))[0]
        else:
            istep = 0
            nx = marker
        x = np.array(struct.unpack(f'{nx}d', f.read(8*nx)))
        ny = struct.unpack('i', f.read(4))[0]
        y = np.array(struct.unpack(f'{ny}d', f.read(8*ny)))
        nz = struct.unpack('i', f.read(4))[0]
        z = np.array(struct.unpack(f'{nz}d', f.read(8*nz)))
        nxm = struct.unpack('i', f.read(4))[0]; xm = f.read(8*nxm)
        nym = struct.unpack('i', f.read(4))[0]; ym = f.read(8*nym)
        nzm = struct.unpack('i', f.read(4))[0]; zm = f.read(8*nzm)
        nxg, nyg, nzg = nxm+2, nym+2, nzm+2
        # U
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        U = np.frombuffer(f.read(8*n1*n2*(n3-1)), dtype=np.float64).copy().reshape((n3-1,n2,n1)).transpose(2,1,0)
        # V
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        V = np.frombuffer(f.read(8*n1*n2*(n3-1)), dtype=np.float64).copy().reshape((n3-1,n2,n1)).transpose(2,1,0)
        # W
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        W = np.frombuffer(f.read(8*n1*n2*(n3-1)), dtype=np.float64).copy().reshape((n3-1,n2,n1)).transpose(2,1,0)
    xm_arr = np.array(struct.unpack(f'{nxm}d', xm))
    ym_arr = np.array(struct.unpack(f'{nym}d', ym))
    xg = np.zeros(nxm+2); xg[1:-1] = xm_arr
    xg[0] = xm_arr[0] - 2*(xm_arr[0]-x[0]); xg[-1] = xm_arr[-1] + 2*(x[-1]-xm_arr[-1])
    yg = np.zeros(nym+2); yg[1:-1] = ym_arr
    yg[0] = ym_arr[0] - 2*(ym_arr[0]-y[0]); yg[-1] = ym_arr[-1] + 2*(y[-1]-ym_arr[-1])
    return {'t': t, 'nu': nu, 'istep': istep, 'x': x, 'y': y, 'z': z,
            'U': U, 'V': V, 'W': W, 'xg': xg, 'yg': yg, 'nx': nx, 'ny': ny, 'nz': nz}

def blasius_solution(eta_max=12, n_points=500):
    def ode(eta, y): return [y[1], y[2], -y[0]*y[2]]
    sol = solve_ivp(ode, (0, eta_max), [0, 0, 0.46960],
                    t_eval=np.linspace(0, eta_max, n_points), rtol=1e-12, atol=1e-12)
    return sol.t, sol.y[1]

# ============================================================
# 1. Inflow profile
# ============================================================
def plot_inflow():
    stats_files = sorted(glob.glob('data/*.stats.txt'))
    if not stats_files: return
    s = read_stats(stats_files[0])
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[0].plot(s['Umean'][0, :], s['y'], 'b-', linewidth=2)
    ax[0].set_xlabel('U / U_inf'); ax[0].set_ylabel('y')
    ax[0].set_title('Inlet U profile')
    ax[0].grid(True, alpha=0.3)
    ax[1].plot(s['Vmean'][0, :], s['y'], 'r-', linewidth=2)
    ax[1].set_xlabel('V'); ax[1].set_ylabel('y')
    ax[1].set_title('Inlet V profile')
    ax[1].grid(True, alpha=0.3)
    fig.suptitle(f'Inflow profiles (step {s["istep"]})')
    fig.tight_layout()
    fig.savefig(f'{figdir}/inflow_profile.png', dpi=150)
    plt.close()
    print('Saved inflow_profile.png')

# ============================================================
# 2. Cf evolution
# ============================================================
def plot_cf():
    stats_files = sorted(glob.glob('data/*.stats.txt'))
    if not stats_files: return
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    eta_bl, df_bl = blasius_solution()

    for sf in stats_files:
        s = read_stats(sf)
        Rex = s['x'] / s['nu']
        step = os.path.basename(sf).split('.')[1]
        ax.plot(Rex[3:-5], s['Cf'][3:-5], '-', alpha=0.5, linewidth=0.8, label=f'step {step}')

    s0 = read_stats(stats_files[0])
    Rex = s0['x'] / s0['nu']
    Cf_bl = 0.664 / np.sqrt(Rex)
    ax.plot(Rex[3:-5], Cf_bl[3:-5], 'k--', linewidth=2, label='Blasius')

    ax.set_xlabel(r'$Re_x$', fontsize=12)
    ax.set_ylabel(r'$C_f$', fontsize=12)
    ax.legend(fontsize=7, ncol=3, loc='upper right')
    ax.set_title('Skin friction evolution')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{figdir}/Cf_evolution.png', dpi=150)
    plt.close()
    print('Saved Cf_evolution.png')

# ============================================================
# 3. Velocity profiles at several x-stations
# ============================================================
def plot_velocity_profiles():
    stats_files = sorted(glob.glob('data/*.stats.txt'))
    if not stats_files: return
    s = read_stats(stats_files[-1])  # latest
    eta_bl, df_bl = blasius_solution()

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    stations = np.linspace(5, min(s['nx']-10, 280), 5, dtype=int)
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(stations)))
    for idx, ix in enumerate(stations):
        Rex = s['x'][ix] / s['nu']
        eta = s['y'] * np.sqrt(1.0 / (2*s['nu']*s['x'][ix]))
        ax.plot(s['Umean'][ix, :], eta, 'o', color=colors[idx], markersize=3,
                label=f'$Re_x={Rex:.0f}$')
    ax.plot(df_bl, eta_bl, 'k-', linewidth=1.5, label='Blasius')
    ax.set_xlabel(r'$U / U_\infty$'); ax.set_ylabel(r'$\eta$')
    ax.set_ylim(0, 10); ax.set_xlim(-0.05, 1.3)
    ax.legend(fontsize=7)
    ax.set_title(f'Velocity profiles (step {s["istep"]})')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{figdir}/velocity_profiles.png', dpi=150)
    plt.close()
    print('Saved velocity_profiles.png')

# ============================================================
# 4. Snapshots (U, V at mid-z plane) for latest snapshot
# ============================================================
def plot_snapshots():
    snap_files = sorted([f for f in glob.glob('data/BL_transition.[0-9]*')
                         if not f.endswith('.stats.txt') and not f.endswith('.restart')])
    if not snap_files: return

    # Plot latest snapshot
    fname = snap_files[-1]
    try:
        snap = read_snapshot(fname)
    except Exception as e:
        print(f'Warning: could not read {fname}: {e}')
        return

    kz_mid = snap['U'].shape[2] // 2
    U_plane = snap['U'][:, :, kz_mid]
    V_plane = snap['V'][:, :, min(kz_mid, snap['V'].shape[2]-1)]
    W_plane = snap['W'][:, :, min(kz_mid, snap['W'].shape[2]-1)]

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    ix_end = -3
    x_plot = snap['x'][:ix_end]
    y_plot = snap['y']

    for idx, (field, name, vlims) in enumerate([
        (U_plane[:ix_end,:], f'U (step {snap["istep"]}, t={snap["t"]:.4f})', (0, 1.2)),
        (V_plane[:ix_end,:], 'V', None),
        (W_plane[:ix_end,:], 'W', None)]):

        if vlims is None:
            vv = max(abs(np.nanmin(field)), abs(np.nanmax(field)), 0.01)
            vlims = (-vv, vv)

        im = axes[idx].imshow(field.T, origin='lower', aspect='auto',
                              extent=[x_plot[0], x_plot[-1], y_plot[0], y_plot[-1]],
                              cmap='RdBu_r', vmin=vlims[0], vmax=vlims[1])
        axes[idx].set_ylabel('y')
        axes[idx].set_title(name)
        plt.colorbar(im, ax=axes[idx])

    fig.tight_layout()
    fig.savefig(f'{figdir}/snapshot_latest.png', dpi=150)
    plt.close()
    print(f'Saved snapshot_latest.png (step {snap["istep"]})')

    # Also plot z-plane (x-z slice at y near BL edge)
    jy = min(snap['ny']//3, snap['U'].shape[1]-1)
    U_xz = snap['U'][:ix_end, jy, :]
    nz_u = U_xz.shape[1]
    z_edges = np.linspace(snap['z'][0], snap['z'][-1], nz_u+1) if snap['nz'] > nz_u else snap['z']

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    im = ax.imshow(U_xz.T, origin='lower', aspect='auto',
                   extent=[snap['x'][0], snap['x'][ix_end], snap['z'][0], snap['z'][-1]],
                   cmap='RdBu_r', vmin=0, vmax=1.2)
    ax.set_xlabel('x'); ax.set_ylabel('z')
    ax.set_title(f'U at y={snap["y"][jy]:.4f} (step {snap["istep"]})')
    plt.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(f'{figdir}/snapshot_xz.png', dpi=150)
    plt.close()
    print(f'Saved snapshot_xz.png (y={snap["y"][jy]:.4f})')

# ============================================================
# 5. Monitor log (divergence, mean U)
# ============================================================
def plot_monitor():
    logfile = 'run_gpu_cufft2.log'
    if not os.path.exists(logfile): return

    steps, divs, meanU, elapsed = [], [], [], []
    import subprocess
    lines = subprocess.run(['strings', logfile], capture_output=True, text=True).stdout.split('\n')
    for i, line in enumerate(lines):
        if 'step number' in line:
            try:
                step = int(line.split(':')[1])
                steps.append(step)
            except: pass
        elif 'Maximum divergence' in line and steps:
            try:
                d = float(line.split(':')[1])
                divs.append(d)
            except: pass
        elif 'Elapsed time' in line and steps:
            try:
                e = float(line.split(':')[1])
                elapsed.append(e)
            except: pass
        elif 'Mean U' in line and 'Mean U :' not in line and steps:
            try:
                u = float(line.split(':')[1])
                meanU.append(u)
            except: pass

    n = min(len(steps), len(divs), len(elapsed))
    if n < 2: return
    steps, divs, elapsed = np.array(steps[:n]), np.array(divs[:n]), np.array(elapsed[:n])

    fig, axes = plt.subplots(2, 1, figsize=(10, 6))

    axes[0].semilogy(steps[1:], np.maximum(divs[1:], 1e-16), 'b.-')
    axes[0].set_ylabel('Max divergence')
    axes[0].set_title('Simulation monitor')
    axes[0].grid(True, alpha=0.3)

    # Per-step time
    dt_steps = np.diff(steps)
    dt_elapsed = np.diff(elapsed)
    ms_per_step = dt_elapsed / dt_steps * 1000
    axes[1].plot(steps[1:], ms_per_step, 'r.-')
    axes[1].set_xlabel('Step')
    axes[1].set_ylabel('ms/step')
    axes[1].set_title(f'Performance (avg: {np.mean(ms_per_step):.2f} ms/step)')
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(f'{figdir}/monitor.png', dpi=150)
    plt.close()
    print('Saved monitor.png')

# ============================================================
if __name__ == '__main__':
    plot_inflow()
    plot_cf()
    plot_velocity_profiles()
    plot_snapshots()
    plot_monitor()
    print(f'\nAll figures saved in {figdir}/')
