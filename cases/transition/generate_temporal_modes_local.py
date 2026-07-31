"""
Generate temporal inflow modes for boundary layer transition.
Solves the Orr-Sommerfeld equation for the Blasius profile to find
unstable TS wave eigenfunctions, then writes the binary file
expected by the BL code (inflow_flag=1).

Binary format (Fortran stream, unformatted):
  ny_inlet        (int32)
  n_modes_inlet   (int32)  -- number of spanwise modes
  m_modes_inlet   (int32)  -- number of temporal modes
  beta_inlet      (float64) -- fundamental spanwise wavenumber
  omega_inlet     (float64) -- fundamental temporal frequency
  ymesh_inlet     (float64 array, ny_inlet)
  zmode_inlet     (float64 array, n_modes_inlet) -- integer indices, scaled by beta later
  tmode_inlet     (float64 array, m_modes_inlet) -- integer indices, scaled by omega later
  qu_inlet_r      (float64 array, ny_inlet x n_modes x m_modes) -- Fortran order
  qu_inlet_i      (float64 array, ...)
  qv_inlet_r      (float64 array, ...)
  qv_inlet_i      (float64 array, ...)
  qw_inlet_r      (float64 array, ...)
  qw_inlet_i      (float64 array, ...)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig
import struct
import os

# ============================================================
# Physical parameters (must match transition.turbb)
# ============================================================
nu = 1e-5
U_inf = 1.0
x_inlet = 1.0       # leading edge distance at inlet
Lz = 0.1083         # spanwise domain size
Ly = 0.41           # wall-normal domain size

Re_x = U_inf * x_inlet / nu   # 100,000
delta_star = 1.7208 * x_inlet / np.sqrt(Re_x)  # displacement thickness
Re_ds = U_inf * delta_star / nu  # Re based on delta*
print(f"Re_x = {Re_x:.0f}")
print(f"delta* = {delta_star:.6f}")
print(f"Re_delta* = {Re_ds:.1f}")

# ============================================================
# 1. Blasius base flow
# ============================================================
def blasius():
    """Solve Blasius: f''' + f*f'' = 0, f(0)=0, f'(0)=0, f''(0)=0.4696"""
    sol = solve_ivp(lambda e, yy: [yy[1], yy[2], -yy[0]*yy[2]],
                    (0, 30), [0, 0, 0.4696],
                    t_eval=np.linspace(0, 30, 2000), rtol=1e-12, atol=1e-12)
    return sol.t, sol.y[0], sol.y[1], sol.y[2]  # eta, f, f', f''

eta_bl, f_bl, fp_bl, fpp_bl = blasius()

# ============================================================
# 2. Chebyshev differentiation matrix
# ============================================================
def cheb_diff(N):
    """Chebyshev differentiation matrix on [-1, 1] with N+1 points."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N+1) / N)
    c = np.ones(N+1)
    c[0] = 2.0; c[N] = 2.0
    c *= (-1.0)**np.arange(N+1)
    X = np.tile(x, (N+1, 1))
    dX = X - X.T
    D = np.outer(c, 1.0/c) / (dX + np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

# ============================================================
# 3. Orr-Sommerfeld eigenvalue solver (temporal)
# ============================================================
def solve_OS(alpha_nd, beta_nd, Re, N=180, y_max_nd=40.0):
    """
    Solve the temporal Orr-Sommerfeld equation for the Blasius profile.

    Non-dimensionalized by delta* and U_inf:
      alpha_nd = alpha * delta*
      beta_nd  = beta * delta*
      Re       = U_inf * delta* / nu
      y_max_nd = y_max / delta*

    Returns: eigenvalues omega (non-dim), eigenvectors v_hat(y)
    """
    k2 = alpha_nd**2 + beta_nd**2

    # Chebyshev grid on [-1, 1]
    D_xi, xi = cheb_diff(N)

    # Map to y in [0, y_max]: y = y_max*(1+xi)/2
    # dy/dxi = y_max/2
    scale = 2.0 / y_max_nd
    D1 = D_xi * scale          # d/dy
    D2 = D1 @ D1               # d²/dy²
    D4 = D2 @ D2               # d⁴/dy⁴

    y_nd = y_max_nd * (1 + xi) / 2  # y/delta*

    # Blasius base flow at these y-points
    eta_scale = 1.7208 / np.sqrt(2)
    eta_pts = y_nd * eta_scale

    U_base = np.interp(eta_pts, eta_bl, fp_bl, right=1.0)
    f_vals = np.interp(eta_pts, eta_bl, f_bl, right=f_bl[-1])
    fpp_vals = np.interp(eta_pts, eta_bl, fpp_bl, right=0.0)
    Upp_base = (-f_vals * fpp_vals) * eta_scale**2  # d²U/dy² (non-dim)

    I = np.eye(N+1)

    # Operators
    L2 = D2 - k2 * I      # D² - k²
    L2sq = L2 @ L2         # (D² - k²)²

    # OS equation: A v = omega B v
    U_diag = np.diag(U_base)
    Upp_diag = np.diag(Upp_base)

    A = 1j * alpha_nd * (U_diag @ L2 - Upp_diag) - (1.0/Re) * L2sq
    B = 1j * L2

    # Boundary conditions: v(0) = v'(0) = 0, v(y_max) = v'(y_max) = 0
    for bc_row in [0, 1, N-1, N]:
        A[bc_row, :] = 0.0
        B[bc_row, :] = 0.0

    A[0, :] = I[0, :]
    A[1, :] = D1[0, :]
    A[N-1, :] = D1[N, :]
    A[N, :] = I[N, :]

    eigvals, eigvecs = eig(A, B)

    return eigvals, eigvecs, y_nd, U_base, Upp_base, D1


def find_unstable_modes(alpha_nd, beta_nd, Re, N=180, y_max_nd=40.0):
    """Find physically meaningful unstable modes."""
    eigvals, eigvecs, y_nd, U_base, Upp_base, D1 = solve_OS(alpha_nd, beta_nd, Re, N, y_max_nd)

    c = eigvals / alpha_nd
    mask = (np.isfinite(eigvals) &
            (np.abs(eigvals) < 100) &
            (c.real > 0) & (c.real < 1.2))

    idx_physical = np.where(mask)[0]
    if len(idx_physical) == 0:
        return [], [], y_nd, U_base, D1

    omega_phys = eigvals[idx_physical]
    sort_idx = np.argsort(-omega_phys.imag)

    return omega_phys[sort_idx], eigvecs[:, idx_physical[sort_idx]], y_nd, U_base, D1


# ============================================================
# 4. Scan for unstable modes
# ============================================================
print("\n--- Scanning for unstable modes ---")
print(f"{'alpha*delta*':>14s} {'omega_r':>10s} {'omega_i':>10s} {'c_r':>8s} {'c_i':>10s}")

best_growth = -999
best_alpha = None
best_omega = None
best_eigvec = None
best_y = None
best_U = None
best_D1 = None

alpha_scan = np.linspace(0.10, 0.40, 31)
for a_nd in alpha_scan:
    omegas, evecs, y_nd, U_base, D1 = find_unstable_modes(a_nd, 0.0, Re_ds)
    if len(omegas) > 0 and omegas[0].imag > 0:
        c = omegas[0] / a_nd
        print(f"  {a_nd:12.4f}   {omegas[0].real:10.6f}   {omegas[0].imag:10.6f}   {c.real:8.4f}   {c.imag:10.6f}  *UNSTABLE*")
        if omegas[0].imag > best_growth:
            best_growth = omegas[0].imag
            best_alpha = a_nd
            best_omega = omegas[0]
            best_eigvec = evecs[:, 0]
            best_y = y_nd
            best_U = U_base
            best_D1 = D1
    else:
        if len(omegas) > 0:
            c = omegas[0] / a_nd
            if a_nd in alpha_scan[::5]:
                print(f"  {a_nd:12.4f}   {omegas[0].real:10.6f}   {omegas[0].imag:10.6f}   {c.real:8.4f}   {c.imag:10.6f}")

if best_alpha is None:
    print("No unstable modes found! Check Re_delta* and alpha range.")
    exit(1)

print(f"\nMost unstable mode:")
print(f"  alpha*delta* = {best_alpha:.4f}")
print(f"  omega = {best_omega.real:.6f} + {best_omega.imag:.6f}i (non-dim)")
print(f"  c = {best_omega.real/best_alpha:.4f} + {best_omega.imag/best_alpha:.6f}i")

# Convert to dimensional
alpha_dim = best_alpha / delta_star
omega_dim = best_omega.real * U_inf / delta_star
F_reduced = omega_dim * nu / U_inf**2
print(f"\nDimensional:")
print(f"  alpha = {alpha_dim:.2f} (1/m)")
print(f"  omega = {omega_dim:.2f} (rad/s)")
print(f"  F = omega*nu/U^2 = {F_reduced:.2e}")
print(f"  wavelength_x = {2*np.pi/alpha_dim:.6f}")


# ============================================================
# 5. Compute mode shapes (u, v, w) from v_hat eigenfunction
# ============================================================
def compute_uvw_from_vhat(v_hat, alpha_nd, beta_nd, omega_nd, Re, y_nd, U_base, D1):
    k2 = alpha_nd**2 + beta_nd**2
    Dv = D1 @ v_hat

    if abs(beta_nd) < 1e-10:
        u_hat = 1j * Dv / alpha_nd
        w_hat = np.zeros_like(v_hat)
    else:
        N = len(v_hat) - 1
        I = np.eye(N+1)
        D2 = D1 @ D1
        L2 = D2 - k2 * I

        eta_scale = 1.7208 / np.sqrt(2)
        eta_pts = y_nd * eta_scale
        Up_base = np.interp(eta_pts, eta_bl, fpp_bl, right=0.0) * eta_scale

        LHS = (-1j*omega_nd + 1j*alpha_nd*np.diag(U_base)) - (1.0/Re)*L2
        RHS = -1j * beta_nd * Up_base * v_hat

        LHS[0, :] = I[0, :]
        RHS[0] = 0.0
        LHS[N, :] = I[N, :]
        RHS[N] = 0.0

        eta_hat = np.linalg.solve(LHS, RHS)

        u_hat = 1j * (alpha_nd * Dv + beta_nd * eta_hat) / k2
        w_hat = 1j * (beta_nd * Dv - alpha_nd * eta_hat) / k2

    return u_hat, v_hat, w_hat


# ============================================================
# 6. Set up modes for the file
# ============================================================
beta_0_dim = 2 * np.pi / Lz
beta_0_nd = beta_0_dim * delta_star
print(f"\nFundamental beta = 2*pi/Lz = {beta_0_dim:.2f} (dim), {beta_0_nd:.4f} (non-dim)")

omega_0_dim = best_omega.real * U_inf / delta_star
print(f"Fundamental omega = {omega_0_dim:.2f} (dim)")

# Spanwise modes: 0 (2D TS), ±1, ±2 (oblique)
z_mode_indices = np.array([0.0, 1.0, -1.0, 2.0, -2.0])
n_modes = len(z_mode_indices)

# Temporal modes: 1, 2, 3 (harmonics of fundamental)
t_mode_indices = np.array([1.0, 2.0, 3.0])
m_modes = len(t_mode_indices)

print(f"\nMode grid: {n_modes} spanwise x {m_modes} temporal = {n_modes*m_modes} modes")

# y-grid for mode shapes
ny_modes = 200
y_modes_nd = np.linspace(0, 40.0, ny_modes)  # non-dim by delta*
y_modes_dim = y_modes_nd * delta_star

# Compute mode shapes
qu_r = np.zeros((ny_modes, n_modes, m_modes), order='F')
qu_i = np.zeros((ny_modes, n_modes, m_modes), order='F')
qv_r = np.zeros((ny_modes, n_modes, m_modes), order='F')
qv_i = np.zeros((ny_modes, n_modes, m_modes), order='F')
qw_r = np.zeros((ny_modes, n_modes, m_modes), order='F')
qw_i = np.zeros((ny_modes, n_modes, m_modes), order='F')

mode_amplitude = 0.002 * U_inf  # 0.2% perturbation per mode

print("\nComputing mode shapes...")
for n_idx in range(n_modes):
    beta_n_nd = z_mode_indices[n_idx] * beta_0_nd
    for m_idx in range(m_modes):
        omega_m_nd = t_mode_indices[m_idx] * best_omega.real

        k2_target = best_alpha**2 + beta_n_nd**2
        alpha_eff = np.sqrt(max(k2_target - beta_n_nd**2, 0.01))

        try:
            omegas_nm, evecs_nm, y_nm, U_nm, D1_nm = find_unstable_modes(
                alpha_eff, abs(beta_n_nd), Re_ds, N=150, y_max_nd=40.0)

            if len(omegas_nm) > 0:
                v_hat = evecs_nm[:, 0]
                omega_os = omegas_nm[0]

                u_hat, v_hat_full, w_hat = compute_uvw_from_vhat(
                    v_hat, alpha_eff, abs(beta_n_nd), omega_os, Re_ds, y_nm, U_nm, D1_nm)

                max_u = np.max(np.abs(u_hat))
                if max_u > 1e-15:
                    u_hat /= max_u
                    v_hat_full /= max_u
                    w_hat /= max_u

                u_hat *= mode_amplitude
                v_hat_full *= mode_amplitude
                w_hat *= mode_amplitude

                if z_mode_indices[n_idx] < 0:
                    u_hat = np.conj(u_hat)
                    v_hat_full = np.conj(v_hat_full)
                    w_hat = -np.conj(w_hat)

                harmonic_decay = 1.0 / t_mode_indices[m_idx]
                u_hat *= harmonic_decay
                v_hat_full *= harmonic_decay
                w_hat *= harmonic_decay

                qu_r[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], u_hat.real[::-1])
                qu_i[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], u_hat.imag[::-1])
                qv_r[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], v_hat_full.real[::-1])
                qv_i[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], v_hat_full.imag[::-1])
                qw_r[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], w_hat.real[::-1])
                qw_i[:, n_idx, m_idx] = np.interp(y_modes_nd, y_nm[::-1], w_hat.imag[::-1])

                status = "unstable" if omega_os.imag > 0 else "stable"
                print(f"  beta={z_mode_indices[n_idx]:+.0f}*beta0, omega={t_mode_indices[m_idx]:.0f}*omega0: "
                      f"alpha_eff={alpha_eff:.3f}, omega_OS=({omega_os.real:.4f},{omega_os.imag:.4f}) [{status}]")
            else:
                print(f"  beta={z_mode_indices[n_idx]:+.0f}*beta0, omega={t_mode_indices[m_idx]:.0f}*omega0: no physical mode found")
        except Exception as e:
            print(f"  beta={z_mode_indices[n_idx]:+.0f}*beta0, omega={t_mode_indices[m_idx]:.0f}*omega0: FAILED ({e})")

# ============================================================
# 7. Write binary file
# ============================================================
outfile = "temporal_modes_inlet.dat"
print(f"\nWriting {outfile}...")

with open(outfile, 'wb') as f:
    f.write(np.int32(ny_modes).tobytes())
    f.write(np.int32(n_modes).tobytes())
    f.write(np.int32(m_modes).tobytes())
    f.write(np.float64(beta_0_dim).tobytes())
    f.write(np.float64(omega_0_dim).tobytes())

    f.write(y_modes_dim.astype(np.float64).tobytes())
    f.write(z_mode_indices.astype(np.float64).tobytes())
    f.write(t_mode_indices.astype(np.float64).tobytes())

    f.write(np.asfortranarray(qu_r).tobytes())
    f.write(np.asfortranarray(qu_i).tobytes())
    f.write(np.asfortranarray(qv_r).tobytes())
    f.write(np.asfortranarray(qv_i).tobytes())
    f.write(np.asfortranarray(qw_r).tobytes())
    f.write(np.asfortranarray(qw_i).tobytes())

print(f"Written {os.path.getsize(outfile)} bytes")

# ============================================================
# 8. Diagnostic plots
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('figures', exist_ok=True)

# Plot the most unstable eigenfunction
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

v_hat = best_eigvec
u_hat, v_hat_full, w_hat = compute_uvw_from_vhat(
    v_hat, best_alpha, 0.0, best_omega, Re_ds, best_y, best_U, best_D1)
max_u = np.max(np.abs(u_hat))
u_hat /= max_u; v_hat_full /= max_u; w_hat /= max_u

ax = axes[0]
ax.plot(u_hat.real, best_y, 'b-', label='Re(û)')
ax.plot(u_hat.imag, best_y, 'r--', label='Im(û)')
ax.set_ylabel('y/δ*'); ax.set_title('û (streamwise)')
ax.legend(); ax.set_ylim(0, 15)

ax = axes[1]
ax.plot(v_hat_full.real, best_y, 'b-', label='Re(v̂)')
ax.plot(v_hat_full.imag, best_y, 'r--', label='Im(v̂)')
ax.set_ylabel('y/δ*'); ax.set_title('v̂ (wall-normal)')
ax.legend(); ax.set_ylim(0, 15)

ax = axes[2]
ax.plot(best_U, best_y, 'k-', lw=2, label='U(y)')
ax.set_ylabel('y/δ*'); ax.set_title('Base flow')
ax.legend(); ax.set_ylim(0, 15); ax.set_xlim(0, 1.1)

fig.suptitle(f'Most unstable TS mode: αδ*={best_alpha:.3f}, ω=({best_omega.real:.4f},{best_omega.imag:.4f}), '
             f'Re_δ*={Re_ds:.0f}', fontsize=12)
fig.tight_layout()
fig.savefig('figures/OS_eigenfunction.png', dpi=150)
plt.close()

# Plot all mode shapes
fig, axes = plt.subplots(n_modes, m_modes, figsize=(4*m_modes, 3*n_modes))
if n_modes == 1: axes = axes[None, :]
if m_modes == 1: axes = axes[:, None]

for n_idx in range(n_modes):
    for m_idx in range(m_modes):
        ax = axes[n_idx, m_idx]
        ax.plot(qu_r[:, n_idx, m_idx], y_modes_nd, 'b-', label='Re(û)')
        ax.plot(qu_i[:, n_idx, m_idx], y_modes_nd, 'r--', label='Im(û)')
        ax.set_ylim(0, 15)
        ax.set_title(f'β={z_mode_indices[n_idx]:+.0f}β₀, ω={t_mode_indices[m_idx]:.0f}ω₀', fontsize=9)
        if n_idx == 0 and m_idx == 0: ax.legend(fontsize=7)
        if m_idx == 0: ax.set_ylabel('y/δ*')

fig.suptitle('Inlet mode shapes û(y) for all (β,ω) combinations', fontsize=12)
fig.tight_layout()
fig.savefig('figures/inlet_mode_shapes.png', dpi=150)
plt.close()

# Neutral curve scan
print("\n--- Neutral curve scan ---")
alpha_range = np.linspace(0.05, 0.50, 46)
growth_rates = np.zeros_like(alpha_range)
frequencies = np.zeros_like(alpha_range)

for i, a_nd in enumerate(alpha_range):
    omegas, _, _, _, _ = find_unstable_modes(a_nd, 0.0, Re_ds, N=150)
    if len(omegas) > 0:
        growth_rates[i] = omegas[0].imag
        frequencies[i] = omegas[0].real

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.plot(alpha_range, growth_rates, 'b-o', ms=3)
ax.axhline(0, color='k', ls='--', lw=0.5)
ax.set_xlabel('αδ*'); ax.set_ylabel('ω_i (growth rate)')
ax.set_title(f'Temporal growth rate (Re_δ*={Re_ds:.0f})')
ax.axvline(best_alpha, color='r', ls='--', alpha=0.5, label=f'most unstable: αδ*={best_alpha:.3f}')
ax.legend()

ax = axes[1]
ax.plot(alpha_range, frequencies, 'b-o', ms=3)
ax.set_xlabel('αδ*'); ax.set_ylabel('ω_r (frequency)')
ax.set_title('TS wave frequency')
ax.axvline(best_alpha, color='r', ls='--', alpha=0.5)

fig.tight_layout()
fig.savefig('figures/OS_neutral_curve.png', dpi=150)
plt.close()

print("\nDone! Generated:")
print(f"  {outfile} — temporal modes binary file")
print(f"  figures/OS_eigenfunction.png — most unstable eigenfunction")
print(f"  figures/inlet_mode_shapes.png — all mode shapes")
print(f"  figures/OS_neutral_curve.png — growth rate vs wavenumber")
