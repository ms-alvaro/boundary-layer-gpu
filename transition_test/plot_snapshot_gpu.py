"""
Read and visualize a binary snapshot from the Fortran incompressible NS solver.
Scales all wall-normal plots by local delta99.

Usage: python3 plot_snapshot.py <snapshot_file> [--fig-dir <dir>]

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
import sys
import argparse

# ============================================================
# Parse arguments
# ============================================================
parser = argparse.ArgumentParser(description='Plot snapshot from BL solver')
parser.add_argument('snap_file', help='Path to binary snapshot file')
parser.add_argument('--fig-dir', default=None,
                    help='Output directory for figures (default: figures/snapshot_NNNNN/)')
args = parser.parse_args()

snap_file = args.snap_file

# ============================================================
# 1. Read binary snapshot
# ============================================================
def read_int(f, n=1):
    data = f.read(4 * n)
    return struct.unpack(f'{n}i', data) if n > 1 else struct.unpack('i', data)[0]

def read_double(f, n):
    data = f.read(8 * n)
    return np.array(struct.unpack(f'{n}d', data))

with open(snap_file, 'rb') as f:
    t, nu = struct.unpack('dd', f.read(16))
    print(f't = {t:.6f},  nu = {nu:.2e}')

    magic = read_int(f)
    assert magic == -73, f'Expected magic -73, got {magic}'
    nstep = read_int(f)
    print(f'nstep = {nstep}')

    nx = read_int(f);  x = read_double(f, nx)
    ny = read_int(f);  y = read_double(f, ny)
    nz = read_int(f);  z = read_double(f, nz)
    print(f'nx={nx}, ny={ny}, nz={nz}')
    print(f'x: [{x[0]:.4f}, {x[-1]:.4f}]')
    print(f'y: [{y[0]:.6f}, {y[-1]:.6f}]')
    print(f'z: [{z[0]:.6f}, {z[-1]:.6f}]')

    nxm = read_int(f); xm = read_double(f, nxm)
    nym = read_int(f); ym = read_double(f, nym)
    nzm = read_int(f); zm = read_double(f, nzm)

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

    u_shape = read_int(f, 3)
    n3_u = u_shape[2]; u_actual = nx * nyg * n3_u  # handle nprocs=1 (n3-1 planes), f'U shape mismatch: {u_shape} vs ({nx},{nyg},{nzg})'
    U = np.frombuffer(f.read(nx * nyg * (n3_u-1) * 8), dtype=np.float64).reshape((nx, nyg, n3_u-1), order='F')
    print(f'U: shape={U.shape}, range=[{U.min():.4f}, {U.max():.4f}]')

    v_shape = read_int(f, 3)
    assert v_shape == (nxg, ny, nzg), f'V shape mismatch: {v_shape} vs ({nxg},{ny},{nzg})'
    V = np.frombuffer(f.read(nxg * ny * nzg * 8), dtype=np.float64).reshape((nxg, ny, nzg), order='F')
    print(f'V: shape={V.shape}, range=[{V.min():.4f}, {V.max():.4f}]')

    w_shape = read_int(f, 3)
    assert w_shape == (nxg, nyg, nz), f'W shape mismatch: {w_shape} vs ({nxg},{nyg},{nz})'
    W = np.frombuffer(f.read(nxg * nyg * nz * 8), dtype=np.float64).reshape((nxg, nyg, nz), order='F')
    print(f'W: shape={W.shape}, range=[{W.min():.4f}, {W.max():.4f}]')

print('Data read successfully.\n')

# ============================================================
# Output directory
# ============================================================
if args.fig_dir:
    fig_dir = args.fig_dir
else:
    fig_dir = f'figures/snapshot_{nstep:07d}'
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Derived quantities
# ============================================================
U_inf = 1.0
Lx = x[-1] - x[0]
Ly = y[-1]
Lz = z[-1]

# delta99 at each x (spanwise-averaged)
def compute_delta99(U_field, x_arr, yg_arr):
    U_xz_avg = U_field[:, :, 1:-1].mean(axis=2)
    delta99 = np.zeros(len(x_arr))
    for i in range(len(x_arr)):
        for j in range(len(yg_arr) - 1, -1, -1):
            if U_xz_avg[i, j] < 0.99 * U_inf:
                if j < len(yg_arr) - 1:
                    frac = (0.99 * U_inf - U_xz_avg[i, j]) / (U_xz_avg[i, j+1] - U_xz_avg[i, j] + 1e-30)
                    delta99[i] = yg_arr[j] + frac * (yg_arr[j+1] - yg_arr[j])
                else:
                    delta99[i] = yg_arr[j]
                break
    return delta99

delta99 = compute_delta99(U, x, yg)
print(f'delta99 range: [{delta99.min():.6f}, {delta99.max():.6f}]')

# Reference BL thicknesses
Rex = U_inf * x / nu
delta99_blasius = 4.91 * x / np.sqrt(Rex + 1e-30)
delta99_turb = 0.37 * x * (Rex + 1e-30)**(-0.2)

# Inlet delta99 and scaled x-coordinate (used in all plots)
d99_inlet = delta99_blasius[0]
x_scaled = x / d99_inlet

# ============================================================
# 2. x-z planes of U at y/delta99 heights
# ============================================================
print('Plotting x-z planes of U ...')

# Use heights scaled by inlet delta99
y_targets_d99 = [0.01, 0.1, 0.5, 1.5]  # multiples of delta99_inlet
y_targets = [yt * d99_inlet for yt in y_targets_d99]
# First target: near-wall (first cell)
y_targets[0] = yg[2]
y_labels = [f'near wall (1st cell, y/δ99={yg[2]/d99_inlet:.2f})',
            f'y/δ99={y_targets_d99[1]:.1f}',
            f'y/δ99={y_targets_d99[2]:.1f}',
            f'y/δ99={y_targets_d99[3]:.1f}']

fig, axes = plt.subplots(len(y_targets), 1, figsize=(14, 2.8 * len(y_targets)), constrained_layout=True)
for idx, (yt, yl) in enumerate(zip(y_targets, y_labels)):
    jy = np.argmin(np.abs(yg - yt))
    ax = axes[idx]
    Uslice = U[:, jy, 1:-1].T
    im = ax.pcolormesh(x_scaled, zg[1:-1], Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
    ax.set_title(f'U at y = {yg[jy]:.5f} ({yl})')
    ax.set_xlabel(r'$x / \delta_{99,inlet}$')
    ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, label='U')
fig.suptitle(f'x-z planes of U  (t={t:.4f}, step={nstep})', fontsize=14)
fig.savefig(os.path.join(fig_dir, 'U_xz_planes.png'), dpi=200)
plt.close(fig)

# ============================================================
# 3. x-y plane of U at mid-z (y scaled by delta99)
# ============================================================
print('Plotting x-y plane of U ...')

kz_mid = nzg // 2
fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
Uslice = U[:, :, kz_mid].T
im = ax.pcolormesh(x_scaled, yg / d99_inlet, Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
ax.plot(x_scaled, delta99 / d99_inlet, 'k--', lw=1, label='δ99 (DNS)')
ax.plot(x_scaled, delta99_blasius / d99_inlet, 'k:', lw=0.8, alpha=0.5, label='δ99 (Blasius)')
ax.plot(x_scaled, delta99_turb / d99_inlet, 'r:', lw=0.8, alpha=0.5, label='δ99 (turbulent)')
ax.set_xlabel(r'$x / \delta_{99,inlet}$')
ax.set_ylabel(r'$y / \delta_{99,inlet}$')
ax.set_title(f'U in x-y plane at z = {zg[kz_mid]:.5f} (mid-span)')
ax.set_ylim([0, min(yg[-1] / d99_inlet, max(5, 3 * delta99.max() / d99_inlet))])
ax.legend(fontsize=9)
plt.colorbar(im, ax=ax, label='U')
fig.savefig(os.path.join(fig_dir, 'U_xy_plane.png'), dpi=200)
plt.close(fig)

# ============================================================
# 4. y-z planes of U at several x-locations (y scaled by local delta99)
# ============================================================
print('Plotting y-z planes of U ...')

x_fracs = [0.1, 0.3, 0.5, 0.7, 0.9]
x_targets = [x[0] + f * (x[-1] - x[0]) for f in x_fracs]

fig, axes = plt.subplots(1, len(x_targets), figsize=(4 * len(x_targets), 5), constrained_layout=True)
for idx, xt in enumerate(x_targets):
    ix = np.argmin(np.abs(x - xt))
    d99_local = max(delta99[ix], d99_inlet)
    ax = axes[idx]
    Uslice = U[ix, :, 1:-1]
    im = ax.pcolormesh(zg[1:-1], yg / d99_local, Uslice, shading='auto', cmap='RdBu_r', vmin=-0.1, vmax=1.1)
    ax.set_title(f'x/δ99,inlet = {x[ix] / d99_inlet:.1f}')
    ax.set_xlabel('z')
    ax.set_ylabel(r'$y / \delta_{99}$')
    ax.set_ylim([0, min(yg[-1] / d99_local, 5)])
    plt.colorbar(im, ax=ax, label='U', shrink=0.8)
fig.suptitle(f'y-z planes of U  (t={t:.4f}, step={nstep})', fontsize=14)
fig.savefig(os.path.join(fig_dir, 'U_yz_planes.png'), dpi=200)
plt.close(fig)

# ============================================================
# 5. Cf vs x
# ============================================================
print('Computing Cf ...')

U_wall = U[:, 1, 1:-1].mean(axis=1)
tau_w = nu * U_wall / yg[1]
Cf = 2.0 * tau_w / U_inf**2

Cf_blasius = 0.664 / np.sqrt(Rex + 1e-30)
Cf_turb = 0.0592 / (Rex + 1e-30)**0.2

i_start = 5
fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
ax.plot(x_scaled[i_start:], Cf[i_start:], 'b-', lw=1.5, label='DNS')
ax.plot(x_scaled[i_start:], Cf_blasius[i_start:], 'k--', lw=1, label='Blasius (laminar)')
ax.plot(x_scaled[i_start:], Cf_turb[i_start:], 'r--', lw=1, label=r'$0.0592\,Re_x^{-1/5}$ (turbulent)')
ax.set_xlabel(r'$x / \delta_{99,inlet}$')
ax.set_ylabel(r'$C_f$')
ax.set_title(f'Skin-friction coefficient (t={t:.4f}, step={nstep})')
ax.set_ylim([0, max(Cf[i_start:].max() * 1.3, 0.005)])
ax.legend()
ax.grid(True, alpha=0.3)

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim()[0] * d99_inlet * U_inf / nu, ax.get_xlim()[1] * d99_inlet * U_inf / nu)
ax2.set_xlabel(r'$Re_x$')

fig.savefig(os.path.join(fig_dir, 'Cf_vs_x.png'), dpi=200)
plt.close(fig)

# ============================================================
# 6. delta99 vs x
# ============================================================
print('Plotting delta99 ...')

fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
ax.plot(x_scaled, delta99 / d99_inlet, 'b-', lw=1.5, label='DNS δ99')
ax.plot(x_scaled, delta99_blasius / d99_inlet, 'k--', lw=1, label='Blasius δ99')
ax.plot(x_scaled, delta99_turb / d99_inlet, 'r--', lw=1, label='Turbulent δ99')
ax.axhline(Ly / d99_inlet, color='gray', ls=':', lw=1, label=f'Ly = {Ly:.3f}')
ax.set_xlabel(r'$x / \delta_{99,inlet}$')
ax.set_ylabel(r'$\delta_{99} / \delta_{99,inlet}$')
ax.set_title('Boundary-layer thickness')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(fig_dir, 'delta99_vs_x.png'), dpi=200)
plt.close(fig)

# ============================================================
# Summary
# ============================================================
print(f'\n--- Summary ---')
print(f'Re_x range: {Rex[0]:.0f} — {Rex[-1]:.0f}')
print(f'δ99 inlet (Blasius): {d99_inlet:.5f}')
print(f'δ99 outlet (DNS):    {delta99[-1]:.5f}')
print(f'δ99 outlet (turb):   {delta99_turb[-1]:.5f}')
print(f'Ly / δ99_turb(out):  {Ly / delta99_turb[-1]:.2f}')
print(f'\nAll figures saved to {fig_dir}')
