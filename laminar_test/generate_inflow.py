import numpy as np
import struct
from scipy.interpolate import splev, splrep
import sys, os

# Add parent pyfiles to path for Re_theta profiles
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pyfiles'))

xi, Ui = 1., 1.

n, kt, kd = 1./7., .016, .16
Ret2Rex = lambda R_: (R_ / kt)**(1 / (1-n))
Rex2Ret = lambda R_: kt * R_**(1-n)

def generate_profile(y_, Reti_in, profdir):
    Rexi = Ret2Rex(Reti_in)
    di = xi * kd / Rexi**n

    ym = .5 * (y_[:-1] + y_[1:])
    yg = np.zeros(ym.size + 2)
    yg[0], yg[1:-1], yg[-1] = -ym[0], ym, 2*y_[-1] - ym[-1]

    datprof = np.loadtxt(os.path.join(profdir, f'Re_theta.{Reti_in}.prof'), skiprows=32)
    yp, Up, Vp = datprof[:, 0], datprof[:, 6], datprof[:, 7]
    Up /= Up[-1]
    yp *= di

    iyg = np.where(yg > yp[-1])[0][0]
    iyp = np.where(y_ > yp[-1])[0][0]

    slu, slv = splrep(yp, Up), splrep(yp, Vp)
    Uint = splev(yg, slu)
    Vint = splev(y_, slv)

    Uint[iyg:] = Uint[iyg - 1]
    Vint[iyp:] = Vint[iyp - 1]

    # Output filename (native endian for gfortran)
    ny = y_.size
    fout = f'Mean_profile_Retheta{Reti_in}_ny{ny}_Ly{y_[-1]:.2f}.dat'
    with open(fout, 'wb') as f:
        f.write(struct.pack('i', ny))  # native endian
        for ui in Uint:
            f.write(struct.pack('d', ui))
        for vi in Vint:
            f.write(struct.pack('d', vi))

    print(f'Wrote {fout} with ny={ny}')
    return fout, Uint, Vint, yg


if __name__ == "__main__":
    Reti = 1100
    Ly = 0.40
    ny = 32
    alpy = 2.5

    y = np.linspace(0., 1., ny)
    y = np.sinh(alpy * y) / np.sinh(alpy) * Ly

    profdir = os.path.join(os.path.dirname(__file__), '..', 'pyfiles')
    fout, Uint, Vint, yg = generate_profile(y, Reti, profdir)
