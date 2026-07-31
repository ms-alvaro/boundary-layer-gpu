import numpy as np
import struct
from scipy.interpolate import splev, splrep

xi, Ui = 1., 1. # normalizing  variables

n, kt, kd  = 1./7., .016, .16
Ret2Rex = lambda R_: ( R_ / kt )**( 1 / (1-n) )
Rex2Ret = lambda R_: kt * R_**(1-n)

def generate_profile( y_, Reti_in ):
    
    Rexi = Ret2Rex( Reti_in )
    di   = xi * kd / Rexi**n  

    ym    = .5 * ( y_[:-1] + y_[1:] )
    yg    = np.zeros( ym.size+2 )
    yg[0], yg[1:-1], yg[-1] = -ym[0], ym, 2*y[-1]-ym[-1]
    
    datprof = np.loadtxt( f'Re_theta.{Reti_in}.prof', skiprows=32)
    yp, Up, Vp = datprof[:,0], datprof[:,6], datprof[:,7]
    Up /= Up[-1]    # velocity at infinity
    yp *= di        # data from profile is normalized with deltai
    
    iyg = np.where(yg>yp[-1])[0][0]
    iyp = np.where(y >yp[-1])[0][0]
    
    slu, slv = splrep( yp, Up ), splrep( yp, Vp )
    Uint = splev( yg, slu ) 
    Vint = splev( y_ , slv ) 
    
    Uint[iyg:] = Uint[iyg-1]
    Vint[iyp:] = Vint[iyp-1]
    
    # output filename
    fout = f'Mean_profile_Retheta{Reti}.dat'
    
    ny = y.size
    with open(fout,'wb') as f:
        f.write( struct.pack( '>i', ny ) )
        for ui in Uint: f.write( struct.pack( '>d', ui ) )
        for vi in Vint: f.write( struct.pack( '>d', vi ) )
    


if __name__ == "__main__":
    ################################################################################
    ## Define input parameters
    Reti    = 1100      # Retau at inlet (necessary for data reading)
    Ly      = 0.5       # Height of domain 

    ny   = 192          # Number of grid points in y
    alpy = 2.5          # Stretching
    
    ################################################################################
    
    y     = np.linspace( 0., 1., ny )
    y     = np.sinh(alpy*y)/np.sinh(alpy) * Ly
    
    generate_profile( y, Reti )

