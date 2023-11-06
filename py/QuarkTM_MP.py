import numpy as np
from itertools import product
from multiprocessing.pool import ThreadPool, Pool
import QuarkTM
import TMQGP as tm

cpu = 4



def get_S_test(omrange, T, iTM, iR, eps1, eps2):
    # map_S = lambda z: tm.sigma_ff_onshell(z, 1e-3, T, iTM, iR, eps1, eps2)

    with Pool(cpu) as p:
        out = p.starmap(tm.sigma_ff_onshell, [[o, 1e-3, T, iTM, iR, eps1, eps2] for o in omrange])

    return np.array(out)