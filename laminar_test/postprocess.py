"""
Post-processing script for laminar Blasius test case.
Generates figures comparing simulation results with analytical Blasius solution.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import struct
import os

# ============================================================
# 1. Read statistics file
# ============================================================
def read_stats(fname):
    """Read the .stats.txt file produced by the BL code."""
    with open(fname, 'r') as f:
        # Header: t, nu, nx, ny, nz, istep
        header = f.readline().strip().lstrip('%').split()
        t = float(header[0])
        nu = float(header[1])
        nx = int(header[2])
        ny = int(header[3])
        nz = int(header[4])
        istep = int(header[5])

        # Cf (nx values)
        Cf = np.array(f.readline().split(), dtype=float)
        # x (nx values)
        x = np.array(f.readline().split(), dtype=float)
        # y (ny values)
        y = np.array(f.readline().split(), dtype=float)

        # Umean (nx x ny)
        Umean = np.zeros((nx, ny))
        for i in range(nx):
            Umean[i, :] = np.array(f.readline().split(), dtype=float)
        # Vmean
        Vmean = np.zeros((nx, ny))
        for i in range(nx):
            Vmean[i, :] = np.array(f.readline().split(), dtype=float)

    return {'t': t, 'nu': nu, 'nx': nx, 'ny': ny, 'nz': nz,
            'Cf': Cf, 'x': x, 'y': y, 'Umean': Umean, 'Vmean': Vmean}


# ============================================================
# 2. Read binary snapshot
# ============================================================
def read_snapshot(fname):
    """Read binary flow snapshot."""
    with open(fname, 'rb') as f:
        t, nu = struct.unpack('dd', f.read(16))

        # Check for magic number
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

        nxm = struct.unpack('i', f.read(4))[0]
        xm = np.array(struct.unpack(f'{nxm}d', f.read(8*nxm)))

        nym = struct.unpack('i', f.read(4))[0]
        ym = np.array(struct.unpack(f'{nym}d', f.read(8*nym)))

        nzm = struct.unpack('i', f.read(4))[0]
        zm = np.array(struct.unpack(f'{nzm}d', f.read(8*nzm)))

        nxg = nxm + 2
        nyg = nym + 2
        nzg = nzm + 2

        def read_field(f):
            """Read a 3D field: header (3 ints) + data (n3-1 planes).
            The code writes n3-1 z-planes (excluding the last ghost cell)."""
            n1, n2, n3 = struct.unpack('iii', f.read(12))
            n3_actual = n3 - 1  # code writes (:,:,1:nz-1) or (:,:,1:nzg-1)
            total = n1 * n2 * n3_actual
            arr = np.frombuffer(f.read(8 * total), dtype=np.float64).copy()
            return arr.reshape((n3_actual, n2, n1)).transpose(2, 1, 0)

        U = read_field(f)
        V = read_field(f)
        W = read_field(f)

    return {'t': t, 'nu': nu, 'x': x, 'y': y, 'z': z,
            'xm': xm, 'ym': ym, 'zm': zm,
            'U': U, 'V': V, 'W': W}


# ============================================================
# 3. Analytical Blasius solution
# ============================================================
def blasius_solution(eta_max=20, n_points=1000):
    """Solve the Blasius ODE: f''' + ff'' = 0."""
    def ode(eta, y):
        return [y[1], y[2], -y[0] * y[2]]

    y0 = [0.0, 0.0, 0.46960]
    sol = solve_ivp(ode, (0, eta_max), y0,
                    t_eval=np.linspace(0, eta_max, n_points),
                    rtol=1e-12, atol=1e-12)
    return sol.t, sol.y[0], sol.y[1]  # eta, f, f'=U/U_inf


# ============================================================
# 4. Generate figures
# ============================================================
def main():
    # Read final stats
    stats = read_stats('data/BL_laminar.00002000.stats.txt')
    nu = stats['nu']
    x = stats['x']
    y = stats['y']
    Cf = stats['Cf']
    Umean = stats['Umean']
    Vmean = stats['Vmean']

    # Analytical Blasius
    eta_bl, f_bl, df_bl = blasius_solution()

    U_inf = 1.0
    figdir = 'figures'
    os.makedirs(figdir, exist_ok=True)

    # ---- Figure 1: Velocity profiles at several x stations ----
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    stations = [10, 30, 50, 70, 90]  # grid indices
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(stations)))

    for idx, ix in enumerate(stations):
        Rex = U_inf * x[ix] / nu
        eta = y * np.sqrt(U_inf / (2 * nu * x[ix]))
        ax.plot(Umean[ix, :] / U_inf, eta, 'o', color=colors[idx],
                markersize=3, label=f'$Re_x = {Rex:.0f}$')

    ax.plot(df_bl, eta_bl, 'k-', linewidth=1.5, label='Blasius')
    ax.set_xlabel(r'$U / U_\infty$', fontsize=12)
    ax.set_ylabel(r'$\eta = y \sqrt{U_\infty / (2 \nu x)}$', fontsize=12)
    ax.set_ylim(0, 8)
    ax.set_xlim(-0.05, 1.1)
    ax.legend(fontsize=8, loc='lower right')
    ax.set_title('Velocity profiles vs Blasius solution')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{figdir}/velocity_profiles.png', dpi=150)
    plt.close()
    print('Saved velocity_profiles.png')

    # ---- Figure 2: Skin friction coefficient Cf vs Rex ----
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    Rex = U_inf * x / nu
    Cf_blasius = 0.664 / np.sqrt(Rex)

    # Exclude first 3 and last 5 points (inlet/outflow buffer zones)
    s = slice(3, -5)
    ax.plot(Rex[s], Cf[s], 'bo', markersize=3, label='Simulation')
    ax.plot(Rex[s], Cf_blasius[s], 'r-', linewidth=1.5, label=r'Blasius: $C_f = 0.664 / \sqrt{Re_x}$')
    ax.set_xlabel(r'$Re_x$', fontsize=12)
    ax.set_ylabel(r'$C_f$', fontsize=12)
    ax.legend(fontsize=10)
    ax.set_title('Skin friction coefficient')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{figdir}/skin_friction.png', dpi=150)
    plt.close()
    print('Saved skin_friction.png')

    # ---- Figure 3: Snapshot of U velocity field (x-y plane, mid-z) ----
    try:
        snap = read_snapshot('data/BL_laminar.00002000')
        nzg = snap['U'].shape[2]
        kz_mid = nzg // 2

        # U at mid-z plane (U is at x-faces, centers in y and z)
        U_plane = snap['U'][:, :, kz_mid]

        # Create coordinate arrays for U (at x-faces, y-centers)
        xg = np.zeros(snap['xm'].size + 2)
        xg[1:-1] = snap['xm']
        xg[0] = snap['xm'][0] - 2*(snap['xm'][0] - snap['x'][0])
        xg[-1] = snap['xm'][-1] + 2*(snap['x'][-1] - snap['xm'][-1])

        yg = np.zeros(snap['ym'].size + 2)
        yg[1:-1] = snap['ym']
        yg[0] = snap['ym'][0] - 2*(snap['ym'][0] - snap['y'][0])
        yg[-1] = snap['ym'][-1] + 2*(snap['y'][-1] - snap['ym'][-1])

        # Exclude last 5 x-stations (outflow buffer)
        ix_end = -5

        fig, axes = plt.subplots(2, 1, figsize=(10, 5))

        # U
        im = axes[0].pcolormesh(snap['x'][:ix_end], yg,
                                U_plane[:ix_end, :].T,
                                cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
        axes[0].set_ylabel('y')
        axes[0].set_title(f'U velocity (t = {snap["t"]:.4f})')
        plt.colorbar(im, ax=axes[0], label='U')

        # V at mid-z
        V_plane = snap['V'][:, :, kz_mid]
        im = axes[1].pcolormesh(xg[:ix_end], snap['y'],
                                V_plane[:ix_end, :].T,
                                cmap='RdBu_r', shading='auto')
        axes[1].set_xlabel('x')
        axes[1].set_ylabel('y')
        axes[1].set_title('V velocity')
        plt.colorbar(im, ax=axes[1], label='V')

        fig.tight_layout()
        fig.savefig(f'{figdir}/snapshot_UVW.png', dpi=150)
        plt.close()
        print('Saved snapshot_UVW.png')

    except Exception as e:
        print(f'Warning: could not read snapshot: {e}')

    # ---- Figure 4: Convergence of Cf over time ----
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    stats_files = sorted([f for f in os.listdir('data') if f.endswith('.stats.txt')])
    for sf in stats_files:
        s = read_stats(f'data/{sf}')
        Rex_s = U_inf * s['x'] / s['nu']
        step = sf.split('.')[1]
        ax.plot(Rex_s[5:-5], s['Cf'][5:-5], '-', alpha=0.6, label=f'step {step}')

    ax.plot(Rex[5:-5], Cf_blasius[5:-5], 'k--', linewidth=2, label='Blasius analytical')
    ax.set_xlabel(r'$Re_x$', fontsize=12)
    ax.set_ylabel(r'$C_f$', fontsize=12)
    ax.legend(fontsize=7, ncol=2)
    ax.set_title('Cf convergence over time steps')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{figdir}/Cf_convergence.png', dpi=150)
    plt.close()
    print('Saved Cf_convergence.png')

    print('\nAll figures saved in figures/')


if __name__ == '__main__':
    main()
