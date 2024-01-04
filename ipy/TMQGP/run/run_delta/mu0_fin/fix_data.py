import h5py
from numpy import *


Tkeys = [160, 200, 300, 400]

for Tkey in Tkeys:
    fname = 'data_single_%i.hdf5'%Tkey
    df = h5py.File(fname, 'r+')
    df.attrs.update({
        'mu' : 1.*Tkey * 1e-3
    })
    df.close()