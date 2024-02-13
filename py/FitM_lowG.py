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


parser = argparse.ArgumentParser()
parser.add_argument('T', type=float)

args = parser.parse_args()

out_folder = '.'

T = args.T

lat = pd.read_csv(os.path.join(os.path.dirname(__file__), "PT.csv_keep"))

iLat = interp1d(lat.x, lat.PT_lat, kind='cubic')
P_target = iLat(T)
# print(Trange_fit)
# exit()

screen = 0.01
suppress = 0.8

G = 6.3
G1 = 8
L = 0.5

# exit()

class TargetMemory:
    def __init__(self):
        self.init = ''

    def reset(self):
        self.init = ''

    def target(self, T, mQ, mG, G, G1, L, screen, suppress):
        fname = 'fit_%.3f_%.16e.h5py'%(T, mQ)
        ret_code = RunHelpers.iterate(fname, T, float(mQ), mG, G, G1, L, screen, suppress, 
                    init=self.init, mode='LO')
        ret_th = RunHelpers.thermo(fname)

        df = h5py.File(fname, 'r')
        pressure = df.attrs['Ptot'] / df.attrs['T']**4

        phi_frac = df.attrs['P_Phi'] / df.attrs['Ptot']

        print('P = ', pressure, 'Plat = ', iLat(df.attrs['T']), 'Phi_frac = ', phi_frac)
        self.init = fname

        return pressure - iLat(df.attrs['T'])


# exit()

t = TargetMemory()

# t.init = './data_single_200.hdf5'



# print(t.target(Trange[i], mQs[i], mGs[i], G, G1, L, screen, suppress))
init = 0.61

mG = 1.4

sol = leastsq(lambda z: t.target(T, z, mG, G, G1, L, screen, suppress), init,
        epsfcn=1e-3)

init = sol[0]

print(sol)
np.savetxt('sol_%.3f.dat'%T, sol[0])