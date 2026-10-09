#!/usr/bin/env python3
"""Mean excess transpiration row for top_potential_vmean (potential-flow matching lid, src/topbc_mod.cuf).

usage: make_top_vmean.py STATS.h5 RESTART OUT.bin [--sigma 0.25]

STATS.h5: statistics of a run of the same case with the old lid (design/compute_stats.py: profiles_x/U_e,
delta_star on grid/x, attrs x_le, nu). RESTART: any restart of the target grid (for x and xm). The excess
displacement D(x) = U_e delta* - delta*_lam(x) (flux-deficit form; delta*_lam = 1.7208 sqrt(nu (x - x_le)))
is smoothed with a Gaussian of width sigma in x (scales below the lid height do not matter at the lid), set to
zero upstream of the plate, interpolated to the V centres and differentiated: dV(xg) = dD/dx, the steady mean
transpiration that the fixed-flux lid must carry in addition to the laminar V_top for the mean free stream to
stay at zero pressure gradient. Stream binary: int32 20261010, int32 nxg, float64 dV(nxg) (ghosts = 0).
"""
import argparse, struct, sys
import numpy as np, h5py
ap = argparse.ArgumentParser(); ap.add_argument("stats"); ap.add_argument("restart"); ap.add_argument("out")
ap.add_argument("--sigma", type=float, default=0.25); a = ap.parse_args()
f = h5py.File(a.stats, "r"); x_le = float(f.attrs["x_le"]); nu = float(f.attrs["nu"])
xs = f["grid/x"][:]; Ue = f["profiles_x/U_e"][:]; ds = f["profiles_x/delta_star"][:]
with open(a.restart, "rb") as g:                      # header: t, nu, magic, step; then x, y, z, xm, ... (int32 n + float64 n)
    g.seek(24); (nx,) = struct.unpack("i", g.read(4)); x = np.frombuffer(g.read(8 * nx))
    for _ in range(2): (n,) = struct.unpack("i", g.read(4)); g.seek(8 * n, 1)
    (nxm,) = struct.unpack("i", g.read(4)); xm = np.frombuffer(g.read(8 * nxm))
xg = np.concatenate([[xm[0] - 2 * (xm[0] - x[0])], xm, [xm[-1] + 2 * (x[-1] - xm[-1])]]); nxg = xg.size
D = Ue * ds - 1.7208 * np.sqrt(np.maximum(nu * (xs - x_le), 0.0)); D[~np.isfinite(D)] = 0.0; D[xs < x_le + 1.0] = 0.0
# Gaussian smoothing on the (uniform) stats grid
dx = xs[1] - xs[0]; n = int(4 * a.sigma / dx); k = np.exp(-0.5 * (np.arange(-n, n + 1) * dx / a.sigma) ** 2); k /= k.sum()
Dp = np.pad(D, n, mode="edge"); Ds = np.convolve(Dp, k, mode="valid")
Dg = np.interp(xg, xs, Ds); dv = np.gradient(Dg, xg); dv[0] = 0.0; dv[-1] = 0.0
with open(a.out, "wb") as g:
    g.write(struct.pack("2i", 20261010, nxg)); g.write(dv.astype("<f8").tobytes())
print(f"{a.out}: nxg {nxg}; D(outlet) {Ds[-1]:.5f}, net dV flux / L {dv[1:-1].mean():.3e}, max |dV| {np.abs(dv).max():.3e}, "
      f"min {dv.min():.3e}; stations Re_x 2e5/3e5/4e5/6e5/8e5: " + " ".join(f"{np.interp(R*nu + x_le, xg, dv):+.2e}" for R in (2e5, 3e5, 4e5, 6e5, 8e5)))
