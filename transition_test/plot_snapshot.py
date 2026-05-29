"""
Read and visualize a binary snapshot from the Fortran incompressible NS solver.

Binary format (stream/unformatted, little-endian):
  - t (float64), nu (float64)
  - magic=-73 (int32), nstep (int32)
  - For each of x, y, z, xm, ym, zm:  n (int32), data (n*float64)
  - For U, V, W, nu_t, P:  shape (3*int32), data (prod(shape)*float64)

Staggered grid:
  U(nx, nyg, nzg)   -- at x-faces, y-centers, z-centers
  V(nxg, ny, nzg)   -- at x-centers, y-faces, z-centers
  W(nxg, nyg, nz)   -- at x-centers, y-centers, z-faces
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import struct
import os

# ============================================================
# 1. Read binary snapshot
# ============================================================
import sys
snap_file = sys.argv[1] if len(sys.argv) > 1 else '/home/alvaroms/tbl/transition/data/BL_transition.00050000'
fig_dir = '/home/alvaroms/tbl/transition/figures/'
os.makedirs(fig_dir, exist_ok=True)

nu_param = 2e-6

def read_int(f, n=1):
    data = f.read(4 * n)
    return struct.unpack(f'{n}i', data) if n > 1 else struct.unpack('i', data)[0]

def read_double(f, n):
    data = f.read(8 * n)
    return np.array(struct.unpack(f'{n}d', data))

with open(snap_file, 'rb') as f:
    # Metadata
    t, nu = struct.unpack('dd', f.read(16))
    print(f't = {t:.6f},  nu = {nu:.2e}')

    magic = read_int(f)
    assert magic == -73, f'Expected magic -73, got {magic}'
    nstep = read_int(f)
    print(f'nstep = {nstep}')

    # Grids: x, y, z
    nx = read_int(f);  x = read_double(f, nx)
    ny = read_int(f);  y = read_double(f, ny)
    nz = read_int(f);  z = read_double(f, nz)
    print(f'nx={nx}, ny={ny}, nz={nz}')
    print(f'x: [{x[0]:.4f}, {x[-1]:.4f}]')
    print(f'y: [{y[0]:.6f}, {y[-1]:.6f}]')
    print(f'z: [{z[0]:.6f}, {z[-1]:.6f}]')

    # Staggered grids: xm, ym, zm
    nxm = read_int(f); xm = read_double(f, nxm)
    nym = read_int(f); ym = read_double(f, nym)
    nzm = read_int(f); zm = read_double(f, nzm)

    # Compute cell-center grids (for U positioning)
    # yg: cell centers in y (nyg = nym + 2 = ny + 1 points, with ghosts)
    nyg = nym + 2
    yg = np.zeros(nyg)
    for j in range(nym):
        yg[j+1] = 0.5 * (y[j] + y[j+1])
    yg[0] = yg[1] - 2.0 * (yg[1] - y[0])
    yg[-1] = yg[-2] + 2.0 * (y[-1] - yg[-2])

    nxg = nxm + 2
    xg = np.zeros(nxg)
    for i in range(nxm):
        xg[i+1] = 0.5 * (x[i] + x[i+1])
    xg[0] = xg[1] - 2.0 * (xg[1] - x[0])
    xg[-1] = xg[-2] + 2.0 * (x[-1] - xg[-2])

    nzg = nzm + 2
    zg = np.zeros(nzg)
    for k in range(nzm):
        zg[k+1] = 0.5 * (z[k] + z[k+1])
    zg[0] = zg[1] - 2.0 * (zg[1] - z[0])
    zg[-1] = zg[-2] + 2.0 * (z[-1] - zg[-2])

    # U field
    u_shape = read_int(f, 3)
    assert u_shape == (nx, nyg, nzg), f'U shape mismatch: {u_shape} vs ({nx},{nyg},{nzg})'
    U = np.frombuffer(f.read(nx * nyg * nzg * 8), dtype=np.float64).reshape((nx, nyg, nzg), order='F')
    print(f'U: shape={U.shape}, range=[{U.min():.4f}, {U.max():.4f}]')

    # V field
    v_shape = read_int(f, 3)
    assert v_shape == (nxg, ny, nzg), f'V shape mismatch: {v_shape} vs ({nxg},{ny},{nzg})'
    V = np.frombuffer(f.read(nxg * ny * nzg * 8), dtype=np.float64).reshape((nxg, ny, nzg), order='F')
    print(f'V: shape={V.shape}, range=[{V.min():.4f}, {V.max():.4f}]')

    # W field
    w_shape = read_int(f, 3)
    assert w_shape == (nxg, nyg, nz), f'W shape mismatch: {w_shape} vs ({nxg},{nyg},{nz})'
    W = np.frombuffer(f.read(nxg * nyg * nz * 8), dtype=np.float64).reshape((nxg, nyg, nz), order='F')
    print(f'W: shape={W.shape}, range=[{W.min():.4f}, {W.max():.4f}]')

print('Data read successfully.\n')

# ============================================================
# Helper: approximate BL thickness delta99 at each x
# ============================================================
U_inf = 1.0

def compute_delta99(U_field, x_arr, yg_arr):
    """Compute delta99 at each x location (spanwise-averaged)."""
    # Average over z (skip ghost cells)
    U_xz_avg = U_field[:, :, 1:-1].mean(axis=2)  # (nx, nyg)
    delta99 = np.zeros(len(x_arr))
    for i in range(len(x_arr)):
        for j in range(len(yg_arr) - 1, -1, -1):
            if U_xz_avg[i, j] < 0.99 * U_inf:
                # Linear interpolation
                if j < len(yg_arr) - 1:
                    frac = (0.99 * U_inf - U_xz_avg[i, j]) / (U_xz_avg[i, j+1] - U_xz_avg[i, j] + 1e-30)
                    delta99[i] = yg_arr[j] + frac * (yg_arr[j+1] - yg_arr[j])
                else:
                    delta99[i] = yg_arr[j]
                break
    return delta99

delta99 = compute_delta99(U, x, yg)
print(f'delta99 range: [{delta99.min():.6f}, {delta99.max():.6f}]')

# ============================================================
# 2. x-z planes of U at several y-heights
# ============================================================
print('Plotting x-z planes of U ...')

# Pick y-heights: near wall, lower BL, mid-BL, upper BL
# Use absolute y values since delta99 ~ Ly
y_targets = [yg[2], 0.001, 0.005, 0.015]
y_labels = ['near wall (1st cell)', 'y=0.001', 'y=0.005', 'y=0.015']

fig, axes = plt.subplots(len(y_targets), 1, figsize=(14, 2.8 * len(y_targets)), constrained_layout=True)
for idx, (yt, yl) in enumerate(zip(y_targets, y_labels)):
    jy = np.argmin(np.abs(yg - yt))
    ax = axes[idx]
    # U[i, j, k] with x along i, z along k (skip ghost cells in z)
    Uslice = U[:, jy, 1:-1].T  # (nz_interior, nx)
    im = ax.pcolormesh(x, zg[1:-1], Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
    ax.set_title(f'U at y = {yg[jy]:.5f} ({yl})')
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, label='U')
fig.suptitle(f'x-z planes of U  (t={t:.4f}, step={nstep})', fontsize=14)
fig.savefig(os.path.join(fig_dir, 'U_xz_planes.png'), dpi=200)
plt.close(fig)

# ============================================================
# 3. x-y plane of U at mid-z
# ============================================================
print('Plotting x-y plane of U ...')

kz_mid = nzg // 2
fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
Uslice = U[:, :, kz_mid].T  # (nyg, nx)
im = ax.pcolormesh(x, yg, Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
ax.plot(x, delta99, 'k--', lw=1, label='delta99')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title(f'U in x-y plane at z = {zg[kz_mid]:.5f} (mid-span)')
ax.set_ylim([0, min(y[-1], 3 * delta99.max())])
ax.legend()
plt.colorbar(im, ax=ax, label='U')
fig.savefig(os.path.join(fig_dir, 'U_xy_plane.png'), dpi=200)
plt.close(fig)

# ============================================================
# 4. y-z planes of U at several x-locations
# ============================================================
print('Plotting y-z planes of U ...')

x_targets = [x[0] + 0.1*(x[-1]-x[0]),
             x[0] + 0.3*(x[-1]-x[0]),
             x[0] + 0.5*(x[-1]-x[0]),
             x[0] + 0.7*(x[-1]-x[0]),
             x[0] + 0.9*(x[-1]-x[0])]

fig, axes = plt.subplots(1, len(x_targets), figsize=(4 * len(x_targets), 5), constrained_layout=True)
for idx, xt in enumerate(x_targets):
    ix = np.argmin(np.abs(x - xt))
    ax = axes[idx]
    Uslice = U[ix, :, 1:-1]  # (nyg, nz_interior)
    im = ax.pcolormesh(zg[1:-1], yg, Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
    ax.set_title(f'x = {x[ix]:.3f}')
    ax.set_xlabel('z')
    ax.set_ylabel('y')
    ax.set_ylim([0, min(y[-1], 3 * delta99[ix])])
    plt.colorbar(im, ax=ax, label='U', shrink=0.8)
fig.suptitle(f'y-z planes of U  (t={t:.4f}, step={nstep})', fontsize=14)
fig.savefig(os.path.join(fig_dir, 'U_yz_planes.png'), dpi=200)
plt.close(fig)

# ============================================================
# 5. Cf vs x
# ============================================================
print('Computing Cf ...')

# Wall shear stress: tau_w = nu * dU/dy |_{y=0}
# U is at y-centers (yg). Wall is at y[0]=0.
# First interior point is yg[1] (yg[0] is ghost below wall).
# For no-slip wall: U(wall)=0, so dU/dy ~ U(yg[1]) / yg[1]
# Average over z (skip ghost cells)
U_wall = U[:, 1, 1:-1].mean(axis=1)  # (nx,)
tau_w = nu_param * U_wall / yg[1]
Cf = 2.0 * tau_w / U_inf**2

# Blasius reference: Cf = 0.664 / sqrt(Rex)
Rex = U_inf * x / nu_param
Cf_blasius = 0.664 / np.sqrt(Rex + 1e-30)

# Turbulent reference: Cf ~ 0.0592 / Rex^(1/5)  (1/5-power law)
Cf_turb = 0.0592 / (Rex + 1e-30)**0.2

# Skip first few points near leading edge to avoid singularity in plot
i_start = 5

fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
ax.plot(x[i_start:], Cf[i_start:], 'b-', lw=1.5, label='DNS')
ax.plot(x[i_start:], Cf_blasius[i_start:], 'k--', lw=1, label='Blasius (laminar)')
ax.plot(x[i_start:], Cf_turb[i_start:], 'r--', lw=1, label=r'$0.0592\,Re_x^{-1/5}$ (turbulent)')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$C_f$')
ax.set_title(f'Skin-friction coefficient (t={t:.4f}, step={nstep})')
ax.set_ylim([0, max(Cf[i_start:].max() * 1.3, 0.005)])
ax.legend()
ax.grid(True, alpha=0.3)

# Secondary x-axis: Rex
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim()[0] * U_inf / nu_param, ax.get_xlim()[1] * U_inf / nu_param)
ax2.set_xlabel(r'$Re_x$')

fig.savefig(os.path.join(fig_dir, 'Cf_vs_x.png'), dpi=200)
plt.close(fig)

# ============================================================
# Bonus: delta99 vs x
# ============================================================
print('Plotting delta99 ...')
delta99_blasius = 4.91 * x / np.sqrt(Rex + 1e-30)

fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
ax.plot(x, delta99, 'b-', lw=1.5, label='DNS delta99')
ax.plot(x, delta99_blasius, 'k--', lw=1, label='Blasius delta99')
ax.set_xlabel('x')
ax.set_ylabel('delta99')
ax.set_title('Boundary-layer thickness')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(fig_dir, 'delta99_vs_x.png'), dpi=200)
plt.close(fig)

print(f'\nAll figures saved to {fig_dir}')
