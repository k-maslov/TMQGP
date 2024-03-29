import argparse
import pickle
import os
import h5py
import matplotlib
from timeit import default_timer as timer

import TMQGP as tm
import QuarkTM
from QuarkTM import Channel, Particle

import tqdm
import numpy as np
from scipy import signal
from scipy.integrate import quad
from syntax_sugar import END, pipe
from syntax_sugar import process_syntax as p
from syntax_sugar import thread_syntax as t
NTHR = 18


out_folder = '.'

Trange = [0.16, 0.2, 0.3, 0.4]

mGs = [1.4, 1.25, 1., 0.9]
mQs = [0.54, 0.52, 0.48, 0.46]

screen = 0.06
suppress = 0.8

G = 6.0
G1 = 8.
L = 0.5

##### Iteration logic ######
init_arg = ''

force_iterate = 0
force_thermo = 1

for T, mQ, mG in zip(Trange, mQs, mGs):
    fname = 'data_single_%i.hdf5'%(int(1e3*T))
    exists = 0
    ### continue broken execution
    if os.path.exists(fname):
        f = h5py.File(fname, 'r')
        if f.attrs['status'] == 'DONE':
            f.close()
            exists = 1
            init_arg = '--init data_single_%i.hdf5'%(int(1e3*T))
        else:
            f.close()
    
    if not exists or force_iterate:
        cmd = f'python3 -m tmqgp_iterate_single {T} {mQ} {mG} {G} {G1} {L} {screen} {suppress} --showtime '
        cmd += init_arg
        cmd += ' --save_iter'
        print('Running ' + cmd)
        ret_code = os.system(cmd)
        if ret_code != 0:
            exit()

    init_arg = '--init data_single_%i.hdf5'%(int(1e3*T))

    ## calculate pressure
    if not exists or force_thermo:
        cmd_th = f'python3 -m tmqgp_thermo_single {fname}'
        print('Running ' + cmd_th)
        os.system(cmd_th)