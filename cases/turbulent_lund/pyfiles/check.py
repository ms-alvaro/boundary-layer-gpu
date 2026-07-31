import BL
import numpy as np
from pathlib import Path
from NGlib import post as ngpost
import matplotlib.pyplot as plt

datapath = Path( '../data' )
nf_flow  = ngpost.NumberedFile( datapath, 'test_swap', ext = None )
nf_stats = ngpost.NumberedFile( datapath, 'test_swap', ext = 'stats.txt' )

ifrl_f = nf_flow.getifrl()
ifrl_s = nf_stats.getifrl()

Xp, Xm, (U, V, W) = BL.readflow( nf_flow.getfullfilename(ifrl_f[-1]), readnup=False )

print( W.mean() )

case_stats = BL.stats( nf_stats )

cols = plt.get_cmap('Reds')(np.linspace(.2,1.,len(ifrl_s[::1])))
for k_, ifr in enumerate(ifrl_s[::1]):

    case_stats.updatedata( ifr )
    case_stats.getCf()
    plt.plot( case_stats.x, case_stats.Cf, c = cols[k_] )
