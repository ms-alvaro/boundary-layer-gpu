"""
Plot initial condition, inflow BC, and domain design for a BL transition simulation.
Reads all parameters from the .turbb config file.

Usage: python3 plot_ic_inflow.py <config_file>
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
import sys
import re

# ============================================================
# Parse config file
# ============================================================
if len(sys.argv) < 2:
    print("Usage: python3 plot_ic_inflow.py <config_file>")
    sys.exit(1)

config_file = sys.argv[1]
params = {}
with open(config_file) as f:
    for line in f:
        line = line.split('!')[0].strip()
        if '=' in line:
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()
            params[key] = val

# Extract parameters
nxyz = params['nxyz'].split()
nx, ny, nz = int(nxyz[0]), int(nxyz[1]), int(nxyz[2])

box = params['boxsize'].split()
Lx, Ly, Lz = float(box[0]), float(box[1]), float(box[2])
alpha_stretch = float(box[3]) if len(box) > 3 else 3.0

nu = float(params['nu'])
cfl_val = float(params['CFL'])
dt = abs(cfl_val)
Amplitude = float(params.get('Amplitude', '0.0'))
inflow_flag = int(params.get('inflow_flag', '1'))
nsteps = int(params['nsteps'])
nsave = int(params['nsave'])

U_inf = 1.0

print(f"Config: {config_file}")
print(f"Grid: {nx} x {ny} x {nz}")
print(f"Domain: Lx={Lx}, Ly={Ly}, Lz={Lz}")
print(f"nu={nu}, dt={dt}")
print(f"inflow_flag={inflow_flag}, Amplitude={Amplitude}")

# ============================================================
# Build grids
# ============================================================
x = np.linspace(0, 1, nx)
x = x[0] * 0 + 1.0 + x * Lx  # x starts at x_inlet

# y: sinh stretching
y_norm = np.linspace(0, 1, ny)
y = np.sinh(alpha_stretch * y_norm) / np.sinh(alpha_stretch) * Ly

z_norm = np.linspace(0, 1, nz)
z = z_norm / z_norm[-2] * Lz if nz > 2 else z_norm * Lz

# cell centers
nyg = ny + 1
ym = 0.5 * (y[:-1] + y[1:])
yg = np.zeros(nyg)
yg[1:-1] = ym
yg[0] = ym[0] - 2 * (ym[0] - y[0])
yg[-1] = ym[-1] + 2 * (y[-1] - ym[-1])

# x_inlet
x_inlet = x[0]
Re_x_inlet = U_inf * x_inlet / nu

# ============================================================
# Blasius solution
# ============================================================
sol = solve_ivp(lambda e, yy: [yy[1], yy[2], -yy[0]*yy[2]], (0, 20), [0, 0, 0.4696],
                t_eval=np.linspace(0, 20, 1000), rtol=1e-12, atol=1e-12)
eta_bl, f_bl, df_bl = sol.t, sol.y[0], sol.y[1]

# BL thicknesses
def delta99_lam(xx):
    return 4.91 * xx / np.sqrt(U_inf * xx / nu + 1e-30)

def delta99_turb(xx):
    Rex = U_inf * xx / nu
    return 0.37 * xx * (Rex + 1e-30)**(-0.2)

def delta_star(xx):
    return 1.7208 * xx / np.sqrt(U_inf * xx / nu + 1e-30)

d99_inlet = delta99_lam(x_inlet)
d99_outlet_lam = delta99_lam(x[-1])
d99_outlet_turb = delta99_turb(x[-1])
ds_inlet = delta_star(x_inlet)
Re_ds_inlet = U_inf * ds_inlet / nu

print(f"\n--- Domain Design ---")
print(f"x_inlet = {x_inlet:.2f}, Re_x_inlet = {Re_x_inlet:.0f}")
print(f"x_outlet = {x[-1]:.2f}, Re_x_outlet = {U_inf * x[-1] / nu:.0f}")
print(f"δ99_inlet (laminar) = {d99_inlet:.5f}")
print(f"δ99_outlet (laminar) = {d99_outlet_lam:.5f}")
print(f"δ99_outlet (turbulent) = {d99_outlet_turb:.5f}")
print(f"δ* inlet = {ds_inlet:.5f}, Re_δ* = {Re_ds_inlet:.0f}")
print(f"Ly / δ99_turb(outlet) = {Ly / d99_outlet_turb:.2f}")
print(f"dx = {Lx/nx:.5f}")

# TS wavelength at most unstable mode (alpha_nd ~ 0.20)
alpha_nd = 0.20
lambda_TS = 2 * np.pi * ds_inlet / alpha_nd
pts_per_TS = lambda_TS / (Lx / nx)
print(f"λ_TS (most unstable) = {lambda_TS:.5f}")
print(f"pts / λ_TS = {pts_per_TS:.1f}")
print(f"Δy_min = {yg[1]:.6f}")
print(f"Δz = {Lz/(nz-2):.6f}")

os.makedirs('figures', exist_ok=True)

# ============================================================
# IC: Blasius profiles
# ============================================================
U_ic = np.ones((nx, nyg))
for ii in range(nx):
    for jj in range(nyg):
        eta = yg[jj] * np.sqrt(U_inf / (2 * nu * x[ii]))
        U_ic[ii, jj] = U_inf * np.interp(eta, eta_bl, df_bl, right=1.0)

V_ic = np.zeros((nx, ny))
for ii in range(nx):
    Rex_loc = U_inf * x[ii] / nu
    for jj in range(ny):
        eta = y[jj] * np.sqrt(U_inf / (2 * nu * x[ii]))
        fp = np.interp(eta, eta_bl, df_bl, right=1.0)
        fv = np.interp(eta, eta_bl, f_bl, right=f_bl[-1])
        V_ic[ii, jj] = U_inf / np.sqrt(2 * Rex_loc) * (eta * fp - fv)

# ---- Plot IC x-y (scaled by delta99) ----
fig, axes = plt.subplots(2, 1, figsize=(14, 7))

im = axes[0].pcolormesh(x, yg / d99_inlet, U_ic.T, cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
axes[0].plot(x, [delta99_lam(xi) / d99_inlet for xi in x], 'k--', lw=1, label='δ99 (Blasius)')
axes[0].set_ylabel(r'$y / \delta_{99,inlet}$')
axes[0].set_title('U — Blasius IC')
axes[0].set_ylim([0, 5])
axes[0].legend()
plt.colorbar(im, ax=axes[0])

im = axes[1].pcolormesh(x, y / d99_inlet, V_ic.T, cmap='RdBu_r', shading='auto')
axes[1].set_xlabel('x')
axes[1].set_ylabel(r'$y / \delta_{99,inlet}$')
axes[1].set_title('V — Blasius IC')
axes[1].set_ylim([0, 5])
plt.colorbar(im, ax=axes[1])

fig.suptitle(f'Initial condition: Blasius (δ99_inlet = {d99_inlet:.5f})', fontsize=14)
fig.tight_layout()
fig.savefig('figures/IC_xy.png', dpi=150)
plt.close()

# ---- Plot IC x-z planes ----
fig, axes = plt.subplots(2, 1, figsize=(14, 5))
for idx, (jj, lab) in enumerate([(2, 'near wall'), (np.argmin(np.abs(yg - 0.5*d99_inlet)), 'y/δ99=0.5')]):
    ax = axes[idx]
    U_plane = np.tile(U_ic[:, jj], (nz, 1))
    im = ax.pcolormesh(x, z, U_plane, cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
    ax.set_title(f'U at y={yg[jj]:.5f} ({lab})')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax)
fig.suptitle('IC: x-z planes (Blasius, uniform in z)', fontsize=14)
fig.tight_layout()
fig.savefig('figures/IC_xz_U.png', dpi=150)
plt.close()

# ============================================================
# Inflow BC
# ============================================================
U_inlet = np.array([U_inf * np.interp(yg[j] * np.sqrt(U_inf / (2 * nu * x[0])),
                    eta_bl, df_bl, right=1.0) for j in range(nyg)])
V_inlet_y = np.array([U_inf / np.sqrt(2 * Re_x_inlet) *
                       (y[j] * np.sqrt(U_inf / (2 * nu * x[0])) *
                        np.interp(y[j] * np.sqrt(U_inf / (2 * nu * x[0])), eta_bl, df_bl, right=1.0) -
                        np.interp(y[j] * np.sqrt(U_inf / (2 * nu * x[0])), eta_bl, f_bl, right=f_bl[-1]))
                       for j in range(ny)])

if inflow_flag == 4:
    # White noise inflow
    np.random.seed(42)
    n_real = 5
    nzg = nz - 1 + 2

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    ax = axes[0]
    ax.plot(U_inlet, yg / d99_inlet, 'k-', lw=2, label='Blasius')
    for r in range(n_real):
        U_noisy = U_inlet + Amplitude * (np.random.rand(nyg) - 0.5)
        ax.plot(U_noisy, yg / d99_inlet, '-', alpha=0.3, label=f'sample {r+1}')
    ax.set_xlabel('U'); ax.set_ylabel(r'$y / \delta_{99,inlet}$')
    ax.set_title(f'Inflow U — {Amplitude*100:.0f}% white noise')
    ax.set_xlim(-0.3, 1.3); ax.set_ylim([0, 5])
    ax.legend(fontsize=7)

    ax = axes[1]
    ax.plot(V_inlet_y, y / d99_inlet, 'k-', lw=2, label='Blasius')
    for r in range(n_real):
        V_noisy = V_inlet_y + Amplitude * (np.random.rand(ny) - 0.5)
        ax.plot(V_noisy, y / d99_inlet, '-', alpha=0.3)
    ax.set_xlabel('V'); ax.set_ylabel(r'$y / \delta_{99,inlet}$')
    ax.set_title('Inflow V'); ax.set_ylim([0, 5])

    ax = axes[2]
    ax.axvline(0, color='k', lw=1, ls='--')
    for r in range(n_real):
        W_noisy = Amplitude * (np.random.rand(nyg) - 0.5)
        ax.plot(W_noisy, yg / d99_inlet, '-', alpha=0.3)
    ax.set_xlabel('W'); ax.set_ylabel(r'$y / \delta_{99,inlet}$')
    ax.set_title('Inflow W'); ax.set_ylim([0, 5])

    fig.suptitle(f'Inflow BC: Blasius + {Amplitude*100:.0f}% white noise (flag={inflow_flag})', fontsize=14)
    fig.tight_layout()
    fig.savefig('figures/inflow_BC_profiles.png', dpi=150)
    plt.close()

    # y-z plane
    fig, ax = plt.subplots(figsize=(8, 5))
    U_yz = np.zeros((nyg, nzg))
    for j in range(nyg):
        for k in range(nzg):
            U_yz[j, k] = U_inlet[j] + Amplitude * (np.random.rand() - 0.5)
    zg_plot = np.linspace(0, Lz, nzg)
    im = ax.pcolormesh(zg_plot, yg / d99_inlet, U_yz, cmap='RdBu_r', shading='auto', vmin=0, vmax=1.1)
    ax.set_xlabel('z'); ax.set_ylabel(r'$y / \delta_{99,inlet}$')
    ax.set_title(f'Inflow U — y-z plane (one realization)')
    ax.set_ylim([0, 5])
    plt.colorbar(im, ax=ax, label='U')
    fig.tight_layout()
    fig.savefig('figures/inflow_BC_yz.png', dpi=150)
    plt.close()

elif inflow_flag == 1:
    # TS modes inflow — just plot the Blasius profiles
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].plot(U_inlet, yg / d99_inlet, 'k-', lw=2)
    axes[0].set_xlabel('U'); axes[0].set_ylabel(r'$y / \delta_{99,inlet}$')
    axes[0].set_title('Inflow U — Blasius base flow')
    axes[0].set_ylim([0, 5])

    axes[1].plot(V_inlet_y, y / d99_inlet, 'k-', lw=2)
    axes[1].set_xlabel('V'); axes[1].set_ylabel(r'$y / \delta_{99,inlet}$')
    axes[1].set_title('Inflow V — Blasius base flow')
    axes[1].set_ylim([0, 5])

    fig.suptitle(f'Inflow BC: Blasius + temporal OS modes (flag={inflow_flag})', fontsize=14)
    fig.tight_layout()
    fig.savefig('figures/inflow_BC_profiles.png', dpi=150)
    plt.close()

# ============================================================
# Domain design plot
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# Top: BL thickness vs x
x_plot = np.linspace(x[0], x[-1], 200)
ax = axes[0]
ax.fill_between(x_plot, 0, [delta99_lam(xi) / d99_inlet for xi in x_plot],
                alpha=0.2, color='blue', label='δ99 laminar')
ax.fill_between(x_plot, 0, [delta99_turb(xi) / d99_inlet for xi in x_plot],
                alpha=0.15, color='red', label='δ99 turbulent')
ax.axhline(Ly / d99_inlet, color='k', ls='-', lw=2, label=f'Ly = {Ly:.3f} ({Ly/d99_inlet:.1f} δ99_inlet)')
ax.set_xlabel('x')
ax.set_ylabel(r'$\delta_{99} / \delta_{99,inlet}$')
ax.set_title('Domain height vs expected BL thickness')
ax.legend()
ax.grid(True, alpha=0.3)

# Bottom: Resolution summary
ax = axes[1]
ax.axis('off')
info = (
    f"Grid: {nx} × {ny} × {nz} = {nx*ny*nz:,} cells\n"
    f"Domain: Lx={Lx}, Ly={Ly}, Lz={Lz}\n"
    f"ν = {nu:.1e}, dt = {dt:.1e}\n"
    f"\n"
    f"Re_x: {Re_x_inlet:.0f} (inlet) → {U_inf * x[-1] / nu:.0f} (outlet)\n"
    f"Re_δ* (inlet) = {Re_ds_inlet:.0f}  (critical ≈ 520)\n"
    f"\n"
    f"δ99_inlet (laminar) = {d99_inlet:.5f}\n"
    f"δ99_outlet (turbulent) = {d99_outlet_turb:.5f}\n"
    f"Ly / δ99_turb(outlet) = {Ly / d99_outlet_turb:.2f}  (need ≥ 2.5)\n"
    f"\n"
    f"dx = {Lx/nx:.5f},  λ_TS = {lambda_TS:.5f},  pts/λ_TS = {pts_per_TS:.1f}  (need ≥ 10)\n"
    f"Δy_min = {yg[1]:.6f}\n"
    f"Δz = {Lz/(nz-2):.6f}\n"
    f"\n"
    f"Inflow: flag={inflow_flag}, Amplitude={Amplitude}\n"
    f"Steps: {nsteps}, save every {nsave} (= {nsave*dt:.3f}s)"
)
ax.text(0.05, 0.95, info, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.suptitle('Domain Design Summary', fontsize=14)
fig.tight_layout()
fig.savefig('figures/domain_design.png', dpi=150)
plt.close()

print(f"\nAll IC/inflow/domain plots saved to figures/")
