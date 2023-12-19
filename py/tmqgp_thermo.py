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
import TM_thermo
# parse the arguments

parser = argparse.ArgumentParser()
parser.add_argument("folder", type=str, help='Folder to process')

args = parser.parse_args()

folder = args.folder
if not os.path.exists(folder):
    raise ValueError('Folder not found')

df = h5py.File(os.path.join(folder, 'data.hdf5'), 'r')
Tkeys = list(df.keys())

Trange = array([1e-3*int(T) for T in Tkeys])
qrange = df.attrs['qrange']
erange = df.attrs['erange']

mQ = df.attrs['mQ']
mG = df.attrs['mG']

pQs = []
pGs = []

for T, Tkey in zip(Trange, Tkeys):
    pQ = Particle(df.attrs['mQ'], qrange, erange, Gtab=df[Tkey]['Q']['G'])
    pG = Particle(df.attrs['mG'], qrange, erange, Gtab=df[Tkey]['G']['G'], stat='b', d=16)

    pQs += [pQ]
    pGs += [pG]

ps_Q = np.array([tm.OmQ_F(T, pt.iImG, pt.iReG) for T, pt in zip(Trange, pQs)])
ps_G = np.array([tm.OmQ_B(T, pt.iImG, pt.iReG) for T, pt in zip(Trange, pGs)])

ps_free_Q = np.array([quad(lambda z: z*z*T*log(1 + exp(-sqrt(mQ**2 + z**2)/T)) / 2/pi**2, 0, 5)[0] for T in Trange])
ps_free_G = np.array([quad(lambda z: -z*z*T*log(1 - exp(-sqrt(mG**2 + z**2)/T)) / 2/pi**2, 0, 5)[0] for T in Trange])

ps_S_Q = []
ps_S_G = []

for T, Tkey, pt in zip(Trange, Tkeys, pQs):
    sigma = df[Tkey]['Q']['S']
    
    iReS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(sigma)))
    iImS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(sigma)))
    
    iImG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(pt.Gtab)))
    iReG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(pt.Gtab)))
    
    ps_S_Q += [tm.OmS_F(T, iImG, iReG, iImS, iReS)]

for T, Tkey, pt in zip(Trange, Tkeys, pGs):
    sigma = df[Tkey]['G']['S']

    iReS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(sigma)))
    iImS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(sigma)))
    
    iImG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(pt.Gtab)))
    iReG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(pt.Gtab)))
    
    ps_S_G += [tm.OmS_B(T, iImG, iReG, iImS, iReS)]

ps_S_Q = np.array(ps_S_Q)
ps_S_G = np.array(ps_S_G)

################################## Calculating Phi #################################
NFS = {
    'qa1' : 3, 'qa8' : 3, 'qq3' : 3, 'qq6' : 3, 'qg3' : 1, 'qg6' : 1, 'qg15' : 1,
    'gq3' : 3, 'gq6' : 3, 'gq15' : 3, 'gg1' : 1, 'gg16' : 1, 'gg27' : 1
}

keys_QQ = ['qa1', 'qa8', 'qq3', 'qq6']
keys_QG = ['qg3', 'qg6', 'qg15']
keys_GG = ['gg1', 'gg16', 'gg27']
keys_GQ = ['gq3', 'gq6', 'gq15']

f_phi = h5py.File(os.path.join(folder, 'Phi.h5py'), 'w')
lmax = 1

LTs = [dict() for T in Trange]
LTs_Q = []
for T, Tkey, LT in zip(Trange, Tkeys, LTs):
    LT_Q = {'QQ':0, 'QG':0}
    
    for keys, _ in zip([keys_QQ, keys_QG], ['QQ', 'QG']):
        for k in keys:
            ds = df[Tkey]['V'][k]['0'].attrs['ds']
            da = df[Tkey]['V'][k]['0'].attrs['da']
            Fa = df[Tkey]['V'][k]['0'].attrs['Fa']
            # print(k, ds, da, Fa, NFS[k])
    
            _lt = 0
            lts = []
            LT[k] = 0
            for l in range(lmax + 1):
                x = array(df[Tkey]['X'][k]['%i'%l])
                v = array(df[Tkey]['V'][k]['%i'%l])
                
                _lt = (2*l + 1) * np.sign(v[1])*4*pi*NFS[k]*ds * da / 6 * v**2 / x * log(1 - x)
                LT_Q[_] += _lt
                LT[k] += _lt
                f_phi.create_dataset(f'{Tkey}/LogT/{k}/{l}', data=_lt)
    LTs_Q += [LT_Q]

LTs_G = []
for T, Tkey, LT in zip(Trange, Tkeys, LTs):
    LT_G = {'GG':0, 'GQ':0}
    
    for keys, _ in zip([keys_GG, keys_GQ], ['GG', 'GQ']):
        for k in keys:
            ds = df[Tkey]['V'][k]['0'].attrs['ds']
            da = df[Tkey]['V'][k]['0'].attrs['da']
            Fa = df[Tkey]['V'][k]['0'].attrs['Fa']
    
            _lt = 0
            lts = []
            LT[k] = 0
            for l in range(lmax + 1):
                x = array(df[Tkey]['X'][k]['%i'%l])
                v = array(df[Tkey]['V'][k]['%i'%l])
                
                _lt = (2*l + 1) * np.sign(v[1])*4*pi*NFS[k]*ds * da / 16 * v**2 / x * log(1 - x)
                LT_G[_] += _lt
                LT[k] += _lt
                f_phi.create_dataset(f'{Tkey}/LogT/{k}/{l}', data=_lt)
    LTs_G += [LT_G]

################################ Calculating Phi ###############################

################################ Log S for quarks ##############################
LogSs_Q = []
Phis_Q = []

for i, T, Tkey, LTs in zip(range(len(Trange)), Trange, Tkeys, LTs_Q):
    ImLogSs_Q = []
    ReLogSs_Q = []
    # iterating over types of diagrams -- QQ and QG
    for key, func, p2 in zip(['QQ', 'QG'], [tm.sigma_ff_onshell, tm.sigma_fb_onshell], [pQs[i], pGs[i]]):
        LT = LTs[key]

        iImLT = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(LT)))

        iEps1 = tm.Interpolator(qrange, pQs[i].om0(qrange), 'cubic')
        iEps2 = tm.Interpolator(qrange, p2.om0(qrange), 'cubic')
        
        ST = array([
            pipe(erange) | p[lambda z: func(z, q, T, iImLT, p2.R, 
                                           iEps1, iEps2)] * NTHR | END
                    for q in tqdm.tqdm(qrange)])
    
        ST = ST.transpose()
    
        ReST = []
    
        for res in (ST.transpose()):
            iImSigma = tm.Interpolator(erange, np.ascontiguousarray(res), 'cubic')
            ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in erange]
            ReST += [ReSigma]

        ReST = np.array(ReST).transpose()

        ImLogSs_Q += [ST]
        ReLogSs_Q += [ReST]

        f_phi.create_dataset(f'{Tkey}/LogS/{key}', data=1j*ST + ReST)

    LogSs_Q += [1j*np.array(ImLogSs_Q) + np.array(ReLogSs_Q)]

################################# Calculating Phi for gluons ###########################

Phis_Q = []

for pt, LogSs, T in zip(pQs, LogSs_Q, Trange):
    iImST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(imag(LogSs), axis=0)))
    iReST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(real(LogSs), axis=0)))
    
    Phi = 0.5 * tm.OmS_F(T, pt.iImG, pt.iReG, iImST, iReST)
    Phis_Q += [Phi]

Phis_Q = array(Phis_Q)
################################ Log S for gluons ##############################

LogSs_G = []
Phis_G = []

for i, T, Tkey, LTs in zip(range(len(Trange)), Trange, Tkeys, LTs_G):
    ImLogSs_G = []
    ReLogSs_G = []
    # iterating over types of diagrams -- GG and GQ
    for key, func, p2 in zip(['GG', 'GQ'], [tm.sigma_bb_onshell, tm.sigma_bf_onshell], [pGs[i], pQs[i]]):
        LT = LTs[key]

        iImLT = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(LT)))

        iEps1 = tm.Interpolator(qrange, pGs[i].om0(qrange), 'cubic')
        iEps2 = tm.Interpolator(qrange, p2.om0(qrange), 'cubic')
        
        ST = array([
            pipe(erange) | p[lambda z: func(z, q, T, iImLT, p2.R, 
                                           iEps1, iEps2)] * NTHR | END
                    for q in tqdm.tqdm(qrange)])
    
        ST = ST.transpose()
    
        ReST = []
    
        for res in (ST.transpose()):
            iImSigma = tm.Interpolator(erange, np.ascontiguousarray(res), 'cubic')
            ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in erange]
            ReST += [ReSigma]

        ReST = np.array(ReST).transpose()

        ImLogSs_G += [ST]
        ReLogSs_G += [ReST]

        f_phi.create_dataset(f'{Tkey}/LogS/{key}', data=1j*ST + ReST)

    LogSs_G += [1j*np.array(ImLogSs_G) + np.array(ReLogSs_G)]


################################# Calculating Phi for gluons ###########################

Phis_G = []

for pt, LogSs, T in zip(pGs, LogSs_G, Trange):
    iImST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(imag(LogSs), axis=0)))
    iReST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(real(LogSs), axis=0)))
    
    Phi = 0.5 * tm.OmS_B(T, pt.iImG, pt.iReG, iImST, iReST)
    Phis_G += [Phi]

Phis_G = array(Phis_G)
################################## Write out the total pressure #########################

Nc = 3
Nf = 3
Ns = 2
Na = 2

P_QP_G = (Nc**2 - 1) * Ns * (ps_G - ps_S_G)
P_QP_Q = Nc * Nf * Ns * Na * (ps_Q - ps_S_Q)

P_Q_G = (Nc**2 - 1) * Ns * (ps_G)
P_Q_Q = Nc * Nf * Ns * Na * (ps_Q)

P_S_G = (Nc**2 - 1) * Ns * (- ps_S_G)
P_S_Q = Nc * Nf * Ns * Na * (- ps_S_Q)

P_Phi_Q = Nc * Nf * Ns * Na * Phis_Q 
P_Phi_G = (Nc**2 - 1) * Ns * Phis_G 

P_Phi = P_Phi_G + P_Phi_Q
P_tot = P_QP_G + P_QP_Q + P_Phi


df_P = pd.DataFrame(array([Trange, P_tot, P_Q_Q, P_Q_G, P_S_Q, P_S_G, P_Phi, P_Phi_Q, P_Phi_G]).transpose(), 
               columns=['T', 'Ptot', 'P_Q_Q', 'P_Q_G', 'P_S_Q', 'P_S_G', 'P_Phi', 'P_Phi_Q', 'P_Phi_G'])

df_P.to_csv(os.path.join(folder, 'pressure.csv'))

f_phi.attrs.update(
    dict(zip(['T', 'Ptot', 'P_Q_Q', 'P_Q_G', 'P_S_Q', 'P_S_G', 'P_Phi', 'P_Phi_Q', 'P_Phi_G'], 
            array([Trange, P_tot, P_Q_Q, P_Q_G, P_S_Q, P_S_G, P_Phi, P_Phi_Q, P_Phi_G])))
)

f_phi.close()

from matplotlib import pyplot as plt
import matplotlib
matplotlib.style.use('publication')

lp, = plt.plot(Trange, P_Phi/Trange**4, label=r'$\frac{1}{ 2 }\Phi$')
lQ, = plt.plot(Trange, P_QP_Q/Trange**4, label='quarks')
plt.plot(Trange, 3*3*2*2*ps_free_Q/Trange**4, ls=':', c=lQ.get_c())
lG, = plt.plot(Trange, P_QP_G/Trange**4, label='gluons')
plt.plot(Trange, 8*2*ps_free_G/Trange**4, ls=':', c=lG.get_c(), label='free')
plt.plot(Trange, (P_QP_Q + P_QP_G)/Trange**4, ls='-.', c='black')
plt.plot(Trange, P_tot/Trange**4, c='black', label='total')
plt.ylim(-0.5, 5)
plt.legend(ncol=2, fontsize=14)

lat = pd.read_csv(os.path.join(os.path.dirname(__file__), "PT.csv"))
plt.plot(lat.x, lat.PT_lat, ls='none', marker='o')
plt.axhline(0, lw=1, ls=':', c='black')
plt.savefig(os.path.join(folder, 'PT.pdf'), bbox_inches='tight')
plt.xlabel('T [GeV]')
plt.ylabel(r'P/T$^4$')