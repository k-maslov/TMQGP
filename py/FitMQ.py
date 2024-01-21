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
import pandas as pd
from scipy import signal
from scipy.integrate import quad
from syntax_sugar import END, pipe
from syntax_sugar import process_syntax as p
from syntax_sugar import thread_syntax as t
from scipy.interpolate import interp1d
NTHR = 18
import RunHelpers

from scipy.optimize import leastsq


out_folder = '.'

Trange = [0.16, 0.2, 0.3, 0.4]

mGs = [1.4, 1.25, 1., 0.9]

mQs = [0.61, 0.6, 0.56, 0.49]

iMQ = interp1d(Trange, mQs, kind='cubic')
iMG = interp1d(Trange, mGs, kind='cubic')

lat = pd.read_csv(os.path.join(os.path.dirname(__file__), "PT.csv"))

iLat = interp1d(lat.x, lat.PT_lat, kind='cubic')

Trange_fit = lat.x[(lat.x > 0.16) & (lat.x < 0.4)].values[::2]
# print(Trange_fit)
# exit()

screen = 0.0065
suppress = 0.8

G = 14.0
G1 = 14.5
L = 0.2

# exit()

class TargetMemory:
    def __init__(self):
        self.init = ''

    def reset(self):
        self.init = ''

    def target(self, T, mQ, mG, G, G1, L, screen, suppress):
        fname = 'fit_%.3f_%.16e.h5py'%(T, mQ)
        ret_code = RunHelpers.iterate(fname, T, float(mQ), mG, G, G1, L, screen, suppress, 
                    init=self.init, mode='HI')
        ret_th = RunHelpers.thermo(fname)

        df = h5py.File(fname, 'r')
        pressure = df.attrs['Ptot'] / df.attrs['T']**4

        self.init = fname
        return pressure - iLat(df.attrs['T'])

i = 1

# exit()

t = TargetMemory()

# t.init = './data_single_200.hdf5'



# print(t.target(Trange[i], mQs[i], mGs[i], G, G1, L, screen, suppress))
init = 0.61
for T in Trange_fit:
    print(iLat(T))

    mG = iMG(T)
    sol = leastsq(lambda z: t.target(T, z, mG, G, G1, L, screen, suppress), init,
            epsfcn=1e-6)

    init = sol[0]

    print(sol)
    np.savetxt('sol_%.3f.dat'%T, sol[0])