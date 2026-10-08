#!/usr/bin/env python3
"""Reference row for top_potential_ref (potential-flow matching lid, src/topbc_mod.cuf).

usage: make_top_reference.py RESTART OUT.bin

RESTART must be a steady laminar solution of the SAME grid computed with the old lid (top_flag = 4 without
top_potential), e.g. the converged state of a zero-perturbation run. The file holds the span mean of U on the
top-cell centre row (index nyg-1, all x faces): stream binary int32 20261008, int32 nx, float64 Uref(nx).
With it the matching acts on the deviation from that discrete base flow, so the steady state is preserved exactly.
"""
import struct, sys
import numpy as np

def main(src, out):
    with open(src, "rb") as f:
        t, nu = struct.unpack("2d", f.read(16))
        magic, step = struct.unpack("2i", f.read(8))
        assert magic == -73, "not a boundary-layer-gpu snapshot"
        n = {}
        for name in ("x", "y", "z", "xm", "ym", "zm"):
            (m,) = struct.unpack("i", f.read(4)); n[name] = m; f.seek(8 * m, 1)
        nx, nyg, nzg = struct.unpack("3i", f.read(12))
        off = f.tell()
    U = np.memmap(src, np.float64, "r", offset=off, shape=(nx, nyg, nzg - 1), order="F")
    row = np.asarray(U[:, nyg - 2, 1:nzg - 2])          # top-cell row, unique z planes 2..nzg-2
    spread = np.abs(row - row.mean(axis=1, keepdims=True)).max()
    uref = row.mean(axis=1)
    with open(out, "wb") as g:
        g.write(struct.pack("2i", 20261008, nx)); g.write(uref.astype("<f8").tobytes())
    print(f"{out}: nx {nx}, t {t:.3f}, step {step}; Uref(1) {uref[0]:.15f}, Uref(nx) {uref[-1]:.15f}; "
          f"max spanwise spread of the source row {spread:.1e} (should be 0 for a 2-D laminar state)")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
