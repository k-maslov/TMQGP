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

fname_fit = '/home/rfrgroup/k.maslov/Numerics/TMQGP/ipy/TMQGP/run/run_delta/fit_fix_G_HD/result.hdf5'
print('Loading fit from ', fname_fit)
df_fit = h5py.File(fname_fit, 'r')

Trange = df_fit.attrs['Trange']
mGs = df_fit.attrs['mGs']
mQs = df_fit.attrs['mQs']

out_folder = '.'

