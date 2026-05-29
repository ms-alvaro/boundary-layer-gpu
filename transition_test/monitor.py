"""
Live monitoring + snapshot plotting for BL transition simulation.
Parses output.log, reads binary snapshots and stats files.

Usage: python monitor.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, re, glob, struct
from scipy.integrate import solve_ivp


# ============================================================
# Parse output log
# ============================================================
def parse_output_log(logfile='output.log'):
    steps, cfs, maxW, maxU, maxV, divs = [], [], [], [], [], []
    if not os.path.exists(logfile):
        return {'step': np.array([]), 'Cf': np.array([]), 'maxW': np.array([]),
                'maxU': np.array([]), 'maxV': np.array([]), 'div': np.array([])}
    with open(logfile, 'rb') as f:
        text = f.read().decode('ascii', errors='ignore')
    lines = text.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if 'step number' in line:
            m = re.search(r'step number\s*:\s*(\d+)', line)
            if m:
                step = int(m.group(1))
                cf_val = wmax = umax = vmax = div_val = None
                for j in range(max(0, i-10), min(len(lines), i+10)):
                    l = lines[j].strip()
                    m2 = re.search(r'Mean Cf:\s*:\s*([\d.eE+-]+)', l)
                    if m2: cf_val = float(m2.group(1))
                    m2 = re.search(r'Maximum W\s*:\s*([\d.eE+-]+)', l)
                    if m2: wmax = float(m2.group(1))
                    m2 = re.search(r'Maximum U\s*:\s*([\d.eE+-]+)', l)
                    if m2: umax = float(m2.group(1))
                    m2 = re.search(r'Maximum V\s*:\s*([\d.eE+-]+)', l)
                    if m2: vmax = float(m2.group(1))
                    m2 = re.search(r'Maximum divergence\s*:\s*([\d.eE+-]+)', l)
                    if m2: div_val = float(m2.group(1))
                steps.append(step)
                cfs.append(cf_val); maxW.append(wmax)
                maxU.append(umax); maxV.append(vmax); divs.append(div_val)
        i += 1
    return {'step': np.array(steps), 'Cf': np.array(cfs, dtype=float),
            'maxW': np.array(maxW, dtype=float), 'maxU': np.array(maxU, dtype=float),
            'maxV': np.array(maxV, dtype=float), 'div': np.array(divs, dtype=float)}


# ============================================================
# Read binary snapshot
# ============================================================
def read_snapshot(fname):
    with open(fname, 'rb') as f:
        t, nu = struct.unpack('dd', f.read(16))
        marker = struct.unpack('i', f.read(4))[0]
        if marker == -73:
            istep = struct.unpack('i', f.read(4))[0]
            nx = struct.unpack('i', f.read(4))[0]
        else:
            nx = marker
        x = np.array(struct.unpack(f'{nx}d', f.read(8*nx)))
        ny = struct.unpack('i', f.read(4))[0]
        y = np.array(struct.unpack(f'{ny}d', f.read(8*ny)))
        nz = struct.unpack('i', f.read(4))[0]
        z = np.array(struct.unpack(f'{nz}d', f.read(8*nz)))
        nxm = struct.unpack('i', f.read(4))[0]; xm = np.array(struct.unpack(f'{nxm}d', f.read(8*nxm)))
        nym = struct.unpack('i', f.read(4))[0]; ym = np.array(struct.unpack(f'{nym}d', f.read(8*nym)))
        nzm = struct.unpack('i', f.read(4))[0]; zm = np.array(struct.unpack(f'{nzm}d', f.read(8*nzm)))
        nxg, nyg, nzg = nxm+2, nym+2, nzm+2
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        U = np.array(struct.unpack(f'{n1*n2*n3}d', f.read(8*n1*n2*n3))).reshape((n3,n2,n1)).transpose(2,1,0)
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        V = np.array(struct.unpack(f'{n1*n2*n3}d', f.read(8*n1*n2*n3))).reshape((n3,n2,n1)).transpose(2,1,0)
        n1,n2,n3 = struct.unpack('iii', f.read(12))
        W = np.array(struct.unpack(f'{n1*n2*n3}d', f.read(8*n1*n2*n3))).reshape((n3,n2,n1)).transpose(2,1,0)
    yg = np.zeros(nyg); yg[1:-1] = ym
    yg[0] = ym[0]-2*(ym[0]-y[0]); yg[-1] = ym[-1]+2*(y[-1]-ym[-1])
    xg = np.zeros(nxg); xg[1:-1] = xm
    xg[0] = xm[0]-2*(xm[0]-x[0]); xg[-1] = xm[-1]+2*(x[-1]-xm[-1])
    return {'t': t, 'nu': nu, 'x': x, 'y': y, 'z': z, 'xg': xg, 'yg': yg,
            'xm': xm, 'ym': ym, 'zm': zm, 'U': U, 'V': V, 'W': W}


# ============================================================
# Read stats file
# ============================================================
def read_stats(fname):
    with open(fname, 'r') as f:
        header = f.readline().strip().lstrip('%').split()
        t, nu = float(header[0]), float(header[1])
        nx, ny = int(header[2]), int(header[3])
        Cf = np.array(f.readline().split(), dtype=float)
        x  = np.array(f.readline().split(), dtype=float)
        y  = np.array(f.readline().split(), dtype=float)
        Umean = np.zeros((nx, ny))
        for i in range(nx): Umean[i,:] = np.array(f.readline().split(), dtype=float)
        Vmean = np.zeros((nx, ny))
        for i in range(nx): Vmean[i,:] = np.array(f.readline().split(), dtype=float)
    return {'t': t, 'nu': nu, 'nx': nx, 'ny': ny, 'Cf': Cf, 'x': x, 'y': y,
            'Umean': Umean, 'Vmean': Vmean}


# ============================================================
# Blasius reference
# ============================================================
def blasius():
    sol = solve_ivp(lambda e,y: [y[1],y[2],-y[0]*y[2]], (0,20), [0,0,0.4696],
                    t_eval=np.linspace(0,20,1000), rtol=1e-12, atol=1e-12)
    return sol.t, sol.y[0], sol.y[1]


# ============================================================
# Plot progress from output.log
# ============================================================
def plot_progress(data, figdir='figures'):
    os.makedirs(figdir, exist_ok=True)
    dt = 2e-5
    time = data['step'] * dt
    n = len(data['step'])
    if n < 2:
        print(f'Only {n} data points, skipping progress plot.')
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    valid = ~np.isnan(data['Cf'])
    if valid.any():
        axes[0,0].plot(time[valid], data['Cf'][valid], 'b.-', markersize=2)
        axes[0,0].axhline(y=data['Cf'][valid][0], color='r', ls='--', alpha=0.5,
                          label=f'Initial = {data["Cf"][valid][0]:.4e}')
        axes[0,0].legend(fontsize=8)
    axes[0,0].set_xlabel('Time'); axes[0,0].set_ylabel('Mean Cf')
    axes[0,0].set_title('Skin friction'); axes[0,0].grid(True, alpha=0.3)

    valid = ~np.isnan(data['maxW'])
    if valid.any():
        axes[0,1].semilogy(time[valid], data['maxW'][valid], 'g.-', markersize=2)
    axes[0,1].set_xlabel('Time'); axes[0,1].set_ylabel('Max |W|')
    axes[0,1].set_title('Spanwise velocity'); axes[0,1].grid(True, alpha=0.3)

    valid = ~np.isnan(data['maxU'])
    if valid.any():
        axes[1,0].plot(time[valid], data['maxU'][valid], 'r.-', markersize=2)
    axes[1,0].set_xlabel('Time'); axes[1,0].set_ylabel('Max U')
    axes[1,0].set_title('Max streamwise velocity'); axes[1,0].grid(True, alpha=0.3)

    valid = ~np.isnan(data['div'])
    if valid.any():
        axes[1,1].semilogy(time[valid], np.abs(data['div'][valid])+1e-16, 'k.-', markersize=2)
    axes[1,1].set_xlabel('Time'); axes[1,1].set_ylabel('Max |div(u)|')
    axes[1,1].set_title('Divergence'); axes[1,1].grid(True, alpha=0.3)

    nsteps_total = 2500000
    fig.suptitle(f'Transition monitor — step {data["step"][-1]:,} / {nsteps_total:,}  '
                 f'(t = {time[-1]:.3f})', fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{figdir}/progress.png', dpi=150)
    plt.close()
    print(f'Saved {figdir}/progress.png')


# ============================================================
# Plot snapshot planes (x-y at mid-z, x-z near wall)
# ============================================================
def plot_snapshot_planes(snap, figdir='figures'):
    os.makedirs(figdir, exist_ok=True)
    t = snap['t']
    step_str = f't{t:.3f}'

    nzg = snap['U'].shape[2]
    kz_mid = nzg // 2

    # --- x-y plane at mid-z ---
    fig, axes = plt.subplots(3, 1, figsize=(14, 9))
    ix_end = -5

    U_xy = snap['U'][:ix_end, :, kz_mid]
    im = axes[0].pcolormesh(snap['x'][:ix_end], snap['yg'], U_xy.T,
                            cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
    axes[0].set_ylabel('y'); axes[0].set_title(f'U  (x-y plane, mid-z, t={t:.3f})')
    plt.colorbar(im, ax=axes[0])

    V_xy = snap['V'][:ix_end, :, kz_mid]
    vmax_v = max(abs(V_xy.min()), abs(V_xy.max())) * 0.8
    im = axes[1].pcolormesh(snap['xg'][:ix_end], snap['y'], V_xy.T,
                            cmap='RdBu_r', shading='auto', vmin=-vmax_v, vmax=vmax_v)
    axes[1].set_ylabel('y'); axes[1].set_title('V')
    plt.colorbar(im, ax=axes[1])

    W_xy = snap['W'][:ix_end, :, nzg//2]
    vmax_w = max(abs(W_xy.min()), abs(W_xy.max())) * 0.8
    if vmax_w < 1e-10: vmax_w = 0.01
    im = axes[2].pcolormesh(snap['xg'][:ix_end], snap['yg'], W_xy.T,
                            cmap='RdBu_r', shading='auto', vmin=-vmax_w, vmax=vmax_w)
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y'); axes[2].set_title('W')
    plt.colorbar(im, ax=axes[2])

    fig.tight_layout()
    fig.savefig(f'{figdir}/snapshot_xy_{step_str}.png', dpi=150)
    plt.close()

    # --- x-z plane near wall (j=3, inside BL) ---
    jj = 3
    fig, axes = plt.subplots(3, 1, figsize=(14, 6))

    U_xz = snap['U'][:ix_end, jj, :]
    im = axes[0].pcolormesh(snap['x'][:ix_end], snap['z'], U_xz.T,
                            cmap='RdBu_r', shading='auto')
    axes[0].set_ylabel('z'); axes[0].set_title(f'U  (x-z plane, y={snap["yg"][jj]:.4f}, t={t:.3f})')
    plt.colorbar(im, ax=axes[0])

    W_xz = snap['W'][:ix_end, jj, :]
    vmax_w = max(abs(W_xz.min()), abs(W_xz.max())) * 0.8
    if vmax_w < 1e-10: vmax_w = 0.01
    im = axes[1].pcolormesh(snap['xg'][:ix_end], snap['yg'][:snap['W'].shape[1]], W_xz.T,
                            cmap='RdBu_r', shading='auto', vmin=-vmax_w, vmax=vmax_w)
    axes[1].set_ylabel('z'); axes[1].set_title('W')
    plt.colorbar(im, ax=axes[1])

    # U fluctuation (subtract spanwise mean)
    U_mean_z = snap['U'][:ix_end, jj, :].mean(axis=1, keepdims=True)
    U_fluct = snap['U'][:ix_end, jj, :] - U_mean_z
    vmax_uf = max(abs(U_fluct.min()), abs(U_fluct.max())) * 0.8
    if vmax_uf < 1e-10: vmax_uf = 0.01
    im = axes[2].pcolormesh(snap['x'][:ix_end], snap['z'], U_fluct.T,
                            cmap='RdBu_r', shading='auto', vmin=-vmax_uf, vmax=vmax_uf)
    axes[2].set_xlabel('x'); axes[2].set_ylabel('z'); axes[2].set_title("U' (fluctuation)")
    plt.colorbar(im, ax=axes[2])

    fig.tight_layout()
    fig.savefig(f'{figdir}/snapshot_xz_{step_str}.png', dpi=150)
    plt.close()

    # --- y-z plane at several x-stations ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    nx_snap = snap['U'].shape[0]
    ix_stations = [5, nx_snap//4, nx_snap//2, 3*nx_snap//4, nx_snap-10, nx_snap-5]

    for col, ix in enumerate(ix_stations[:3]):
        U_yz = snap['U'][ix, :, :]
        im = axes[0, col].pcolormesh(snap['z'], snap['yg'], U_yz,
                                     cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
        axes[0, col].set_title(f'U at x={snap["x"][min(ix,len(snap["x"])-1)]:.2f}')
        axes[0, col].set_xlabel('z'); axes[0, col].set_ylabel('y')
        plt.colorbar(im, ax=axes[0, col])

    for col, ix in enumerate(ix_stations[3:]):
        U_yz = snap['U'][min(ix, nx_snap-1), :, :]
        im = axes[1, col].pcolormesh(snap['z'], snap['yg'], U_yz,
                                     cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
        x_val = snap['x'][min(ix, len(snap['x'])-1)]
        axes[1, col].set_title(f'U at x={x_val:.2f}')
        axes[1, col].set_xlabel('z'); axes[1, col].set_ylabel('y')
        plt.colorbar(im, ax=axes[1, col])

    fig.suptitle(f'y-z cross-sections of U (t={t:.3f})', fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{figdir}/snapshot_yz_{step_str}.png', dpi=150)
    plt.close()

    print(f'Saved snapshot planes for t={t:.3f}')


# ============================================================
# Plot stats (Cf vs x, U profiles)
# ============================================================
def plot_stats(figdir='figures'):
    os.makedirs(figdir, exist_ok=True)
    eta_bl, f_bl, df_bl = blasius()
    stats_files = sorted(glob.glob('data/*.stats.txt'))
    if not stats_files:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(stats_files)))

    for idx, sf in enumerate(stats_files):
        s = read_stats(sf)
        Rex = s['x'] / s['nu']
        step = sf.split('.')[-3]

        axes[0].plot(Rex[3:-5], s['Cf'][3:-5], '-', color=colors[idx],
                     linewidth=1.5, label=f'step {step} (t={s["t"]:.2f})')

        ix_mid = s['nx'] // 2
        eta = s['y'] * np.sqrt(1 / (2 * s['nu'] * s['x'][ix_mid]))
        axes[1].plot(s['Umean'][ix_mid, :], eta, 'o-', color=colors[idx],
                     markersize=2, linewidth=1, label=f'step {step}')

    Cf_bl = 0.664 / np.sqrt(Rex)
    axes[0].plot(Rex[3:-5], Cf_bl[3:-5], 'k--', linewidth=1.5, label='Blasius')
    axes[0].set_xlabel(r'$Re_x$'); axes[0].set_ylabel(r'$C_f$')
    axes[0].set_title('Skin friction vs x')
    axes[0].legend(fontsize=6, ncol=2); axes[0].grid(True, alpha=0.3)

    axes[1].plot(df_bl, eta_bl, 'k-', linewidth=2, label='Blasius')
    axes[1].set_xlabel(r'$U/U_\infty$'); axes[1].set_ylabel(r'$\eta$')
    axes[1].set_ylim(0, 8); axes[1].set_title('U profile at mid-domain')
    axes[1].legend(fontsize=6, ncol=2); axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(f'{figdir}/stats.png', dpi=150)
    plt.close()
    print(f'Saved {figdir}/stats.png')


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    data = parse_output_log()
    os.makedirs('figures', exist_ok=True)

    if len(data['step']) > 1:
        plot_progress(data)
    elif len(data['step']) == 1:
        print(f'Step {data["step"][0]}: Cf={data["Cf"][0]:.4e}, maxW={data["maxW"][0]:.4e}')

    # Plot stats
    plot_stats()

    # Plot snapshot planes for each available binary snapshot
    snap_files = sorted(glob.glob('data/BL_transition.[0-9]*'))
    snap_files = [f for f in snap_files if not f.endswith('.txt')]
    for sf in snap_files:
        try:
            snap = read_snapshot(sf)
            plot_snapshot_planes(snap)
        except Exception as e:
            print(f'Error reading {sf}: {e}')

    if len(data['step']) > 0:
        print(f'\nLatest: step {data["step"][-1]:,}, Cf={data["Cf"][-1]:.4e}, '
              f'maxW={data["maxW"][-1]:.4e}, maxU={data["maxU"][-1]:.4e}')
