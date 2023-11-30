import argparse
import pickle
import os
import h5py
import matplotlib
from timeit import default_timer as timer
import pandas as pd
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
from numpy import log, sqrt, exp, pi, real, imag, array
from multiprocess import pool


def get_S(qrange, erange, T, func, iImLT, p2, iEps1, iEps2):
    ST = array([
            pipe(erange) | p[lambda z: func(z, q, T, iImLT, p2.R, 
                                        iEps1, iEps2)] * NTHR | END
                    for q in tqdm.tqdm(qrange)])

    return ST