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
from matplotlib import pyplot as plt
import matplotlib
# matplotlib.style.use('publication')
# parse the arguments

# print(os.path.dirname(__file__))
# exit()

parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help='Data file name')

args = parser.parse_args()

file = args.file
if not os.path.exists(file):
    raise ValueError('Data file not found')

df = h5py.File(file, 'r+')

df_phi = h5py.File('th_'+file, 'w')

Trange = [df.attrs['T']]
mu = df.attrs['mu']

qrange = df.attrs['qrange']
erange = df.attrs['erange']

qrange2b = df.attrs['qrange2b']
erange2b = df.attrs['erange2b']

# mQ = df.attrs['mQ']
mQs = [df.attrs['mQ']]
# mG = df.attrs['mG']
mGs = [df.attrs['mG']]

pQs = []
# pGs = []
pAs = []

for T, mQ, mG in zip(Trange, mQs, mGs):
    pQ = Particle(mQ, qrange, erange, Gtab=np.array(df['Q']['G']), S=np.array(df['Q']['S']), mu=mu)
    pA = Particle(mQ, qrange, erange, Gtab=np.array(df['A']['G']), S=np.array(df['Q']['S']), mu=-mu)
    # pG = Particle(mG, qrange, erange, Gtab=df['G']['G'], stat='b', d=16)

    pQs += [pQ]
    # pGs += [pG]
    pAs += [pA]

ps_Q = np.array([tm.OmQ_F(T, pt.iImG, pt.iReG) for T, pt in zip(Trange, pQs)])
# print(pQs[0].Gtab[:, 0])
# print(ps_Q)
ps_A = np.array([tm.OmQ_F(T, pt.iImG, pt.iReG) for T, pt in zip(Trange, pAs)])
# ps_G = np.array([tm.OmQ_B(T, pt.iImG, pt.iReG) for T, pt in zip(Trange, pGs)])

ps_free_Q = np.array([quad(lambda z: z*z*T*log(1 + exp(-(sqrt(mQ**2 + z**2) - mu)/T)) / 2/pi**2, 0, 5)[0] for T, mQ in zip(Trange, mQs)])
ps_free_A = np.array([quad(lambda z: z*z*T*log(1 + exp(-(sqrt(mQ**2 + z**2) + mu)/T)) / 2/pi**2, 0, 5)[0] for T, mQ in zip(Trange, mQs)])
# ps_free_G = np.array([quad(lambda z: -z*z*T*log(1 - exp(-sqrt(mG**2 + z**2)/T)) / 2/pi**2, 0, 5)[0] for T, mG in zip(Trange, mGs)])

ps_S_Q = []
ps_S_A = []
# ps_S_G = []

for T, pt in zip(Trange, pQs):
    sigma = df['Q']['S']
    
    iReS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(sigma)))
    iImS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(sigma)))
    
    # iImG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(pt.Gtab)))
    # iReG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(pt.Gtab)))
    
    ps_S_Q += [tm.OmS_F(T, pt.iImG, pt.iReG, iImS, iReS)]

for T, pt in zip(Trange, pAs):
    sigma = df['A']['S']
    
    iReS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(sigma)))
    iImS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(sigma)))
    
    # iImG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(pt.Gtab)))
    # iReG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(pt.Gtab)))
    
    ps_S_A += [tm.OmS_F(T, pt.iImG, pt.iReG, iImS, iReS)]


# for T, pt in zip(Trange, pGs):
#     sigma = df['G']['S']

#     iReS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(sigma)))
#     iImS = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(sigma)))
    
#     iImG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(pt.Gtab)))
#     iReG = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(real(pt.Gtab)))
    
#     ps_S_G += [tm.OmS_B(T, iImG, iReG, iImS, iReS)]

ps_S_Q = np.array(ps_S_Q)
ps_S_A = np.array(ps_S_A)
# ps_S_G = np.array(ps_S_G)

################################## Calculating Phi #################################
NFS = {
    'qa1' : 3, 'qa8' : 3, 'qq3' : 3, 'qq6' : 3, 'qg3' : 1, 'qg6' : 1, 'qg15' : 1,
    'gq3' : 3, 'gq6' : 3, 'gq15' : 3, 'gg1' : 1, 'gg16' : 1, 'gg27' : 1,
    'aq1' : 3, 'aq8' : 3, 'aa3' : 3, 'aa6' : 3, 'ag3' : 1, 'ag6' : 1, 'ag15' : 1,
    'ga3' : 3, 'ga6' : 3, 'ga15' : 3,
}

keys_QQ = ['qq3']
keys_QA = ['qa1']
# keys_QG = ['qg3', 'qg6', 'qg15']

keys_AA = ['aa3']
keys_AQ = ['aq1']
# keys_AG = ['ag3', 'ag6', 'ag15']

# keys_GG = ['gg1', 'gg16', 'gg27']
# keys_GQ = ['gq3', 'gq6', 'gq15']
# keys_GA = ['ga3', 'ga6', 'ga15']

# f_phi = h5py.File(os.path.join(folder, 'Phi.h5py'), 'w')
lmax = 1

###################### Populating the LogT array ##########################

LTs = [dict() for T in Trange]
LTs_Q = []
for T, LT in zip(Trange, LTs):
    LT_Q = {'QQ':0, 'QA': 0}
    
    for keys, _ in zip([keys_QQ, keys_QA], ['QQ', 'QA']):
        for k in keys:
            ds = df['V'][k]['0'].attrs['ds']
            da = df['V'][k]['0'].attrs['da']
            Fa = df['V'][k]['0'].attrs['Fa']
            # print(k, ds, da, Fa, NFS[k])
    
            _lt = 0
            lts = []
            LT[k] = 0
            for l in range(lmax + 1):
                x = array(df['X'][k]['%i'%l])
                v = array(df['V'][k]['%i'%l])
                
                _lt = (2*l + 1) * np.sign(v[1])*4*pi*NFS[k]*ds * da / 6 * v**2 / x * log(1 - x)
                LT_Q[_] += _lt
                LT[k] += _lt

                if k == 'qa1' and l == 0:
                    print(np.max(real(x)))
                    print(np.max(real(_lt)))
                    print(k, ds, da, Fa, NFS[k])

                df_phi.create_dataset(f'LogT/{k}/{l}', data=_lt)
    LTs_Q += [LT_Q]

LTs_A = []
for T, LT in zip(Trange, LTs):
    LT_A = {'AQ':0, 'AA' : 0, 'AG':0}
    
    for keys, _ in zip([keys_AQ, keys_AA], ['AQ', 'AA']):
        for k in keys:
            ds = df['V'][k]['0'].attrs['ds']
            da = df['V'][k]['0'].attrs['da']
            Fa = df['V'][k]['0'].attrs['Fa']
            # print(k, ds, da, Fa, NFS[k])
    
            _lt = 0
            lts = []
            LT[k] = 0
            for l in range(lmax + 1):
                x = array(df['X'][k]['%i'%l])
                v = array(df['V'][k]['%i'%l])
                
                _lt = (2*l + 1) * np.sign(v[1])*4*pi*NFS[k]*ds * da / 6 * v**2 / x * log(1 - x)
                LT_A[_] += _lt
                LT[k] += _lt

                df_phi.create_dataset(f'LogT/{k}/{l}', data=_lt)
    LTs_A += [LT_A]

LTs_G = []
for T, LT in zip(Trange, LTs):
    LT_G = {'GG':0, 'GA' : 0, 'GQ':0}
    
    # for keys, _ in zip([keys_GG, keys_GQ, keys_GA], ['GG', 'GQ', 'GA']):
    #     for k in keys:
    #         ds = df['V'][k]['0'].attrs['ds']
    #         da = df['V'][k]['0'].attrs['da']
    #         Fa = df['V'][k]['0'].attrs['Fa']
    
    #         _lt = 0
    #         lts = []
    #         LT[k] = 0
    #         for l in range(lmax + 1):
    #             x = array(df['X'][k]['%i'%l])
    #             v = array(df['V'][k]['%i'%l])
                
    #             _lt = (2*l + 1) * np.sign(v[1])*4*pi*NFS[k]*ds * da / 16 * v**2 / x * log(1 - x)
    #             LT_G[_] += _lt
    #             LT[k] += _lt
    #             df_phi.create_dataset(f'LogT/{k}/{l}', data=_lt)
    LTs_G += [LT_G]

################################ Calculating Phi ###############################

################################ Log S for quarks ##############################
LogSs_Q = []
Phis_Q = []

for i, T, LTs in zip(range(len(Trange)), Trange, LTs_Q):
    ImLogSs_Q = []
    ReLogSs_Q = []
    # iterating over types of diagrams -- QQ and QG
    for key, func, p2 in zip(['QQ', 'QA'], 
                [tm.sigma_ff_onshell, tm.sigma_ff_onshell], 
                [pQs[i], pAs[i]]):
        LT = LTs[key]

        iImLT = tm.Interpolator2D(qrange2b, erange2b, np.ascontiguousarray(imag(LT)))

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

        df_phi.create_dataset(f'LogS/{key}', data=1j*ST + ReST)

    LogSs_Q += [1j*np.array(ImLogSs_Q) + np.array(ReLogSs_Q)]

################################ Log S for antiquarks ##############################
LogSs_A = []
Phis_A = []

for i, T, LTs in zip(range(len(Trange)), Trange, LTs_A):
    ImLogSs_A = []
    ReLogSs_A = []
    # iterating over types of diagrams -- QQ and QG
    for key, func, p2 in zip(['AA', 'AQ'], 
                [tm.sigma_ff_onshell, tm.sigma_ff_onshell], 
                [pAs[i], pQs[i]]):
        LT = LTs[key]

        iImLT = tm.Interpolator2D(qrange2b, erange2b, np.ascontiguousarray(imag(LT)))

        iEps1 = tm.Interpolator(qrange, pAs[i].om0(qrange), 'cubic')
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

        ImLogSs_A += [ST]
        ReLogSs_A += [ReST]

        df_phi.create_dataset(f'LogS/{key}', data=1j*ST + ReST)

    LogSs_A += [1j*np.array(ImLogSs_A) + np.array(ReLogSs_A)]



################################# Calculating Phi for quarks and antiquarks ###########################

Phis_Q = []

for pt, LogSs, T in zip(pQs, LogSs_Q, Trange):
    iImST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(imag(LogSs), axis=0)))
    iReST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(real(LogSs), axis=0)))
    
    Phi = 0.5 * tm.OmS_F(T, pt.iImG, pt.iReG, iImST, iReST)
    Phis_Q += [Phi]

Phis_Q = array(Phis_Q)

Phis_A = []

for pt, LogSs, T in zip(pAs, LogSs_A, Trange):
    iImST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(imag(LogSs), axis=0)))
    iReST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(real(LogSs), axis=0)))
    
    Phi = 0.5 * tm.OmS_F(T, pt.iImG, pt.iReG, iImST, iReST)
    Phis_A += [Phi]

Phis_A = array(Phis_A)
################################ Log S for gluons ##############################

# LogSs_G = []
# Phis_G = []

# for i, T, LTs in zip(range(len(Trange)), Trange, LTs_G):
#     ImLogSs_G = []
#     ReLogSs_G = []
#     # iterating over types of diagrams -- GG and GQ
#     for key, func, p2 in zip(['GG', 'GQ', 'GA'], 
#     [tm.sigma_bb_onshell, tm.sigma_bf_onshell, tm.sigma_bf_onshell], 
#     [pGs[i], pQs[i], pAs[i]]):
#         LT = LTs[key]

#         iImLT = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(imag(LT)))

#         iEps1 = tm.Interpolator(qrange, pGs[i].om0(qrange), 'cubic')
#         iEps2 = tm.Interpolator(qrange, p2.om0(qrange), 'cubic')
        
#         ST = array([
#             pipe(erange) | p[lambda z: func(z, q, T, iImLT, p2.R, 
#                                            iEps1, iEps2)] * NTHR | END
#                     for q in tqdm.tqdm(qrange)])
    
#         ST = ST.transpose()
    
#         ReST = []
    
#         for res in (ST.transpose()):
#             iImSigma = tm.Interpolator(erange, np.ascontiguousarray(res), 'cubic')
#             ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in erange]
#             ReST += [ReSigma]

#         ReST = np.array(ReST).transpose()

#         ImLogSs_G += [ST]
#         ReLogSs_G += [ReST]

#         df_phi.create_dataset(f'LogS/{key}', data=1j*ST + ReST)

#     LogSs_G += [1j*np.array(ImLogSs_G) + np.array(ReLogSs_G)]


# ################################# Calculating Phi for gluons ###########################

# Phis_G = []

# for pt, LogSs, T in zip(pGs, LogSs_G, Trange):
#     iImST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(imag(LogSs), axis=0)))
#     iReST = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.sum(real(LogSs), axis=0)))
    
#     Phi = 0.5 * tm.OmS_B(T, pt.iImG, pt.iReG, iImST, iReST)
#     Phis_G += [Phi]

Phis_G = array([0.])
################################## Write out the total pressure #########################

Nc = 3
Nf = 3
Ns = 2
Na = 1 ### antiquarks taken explicitly

# P_QP_G = (Nc**2 - 1) * Ns * (ps_G - ps_S_G)
P_QP_Q = Nc * Nf * Ns * Na * (ps_Q - ps_S_Q)
P_QP_A = Nc * Nf * Ns * Na * (ps_A - ps_S_A)


# P_Q_G = (Nc**2 - 1) * Ns * (ps_G)
P_Q_Q = Nc * Nf * Ns * Na * (ps_Q)
P_Q_A = Nc * Nf * Ns * Na * (ps_A)

# P_S_G = (Nc**2 - 1) * Ns * (- ps_S_G)
P_S_Q = Nc * Nf * Ns * Na * (- ps_S_Q)
P_S_A = Nc * Nf * Ns * Na * (- ps_S_A)

P_Phi_Q = Nc * Nf * Ns * Na * Phis_Q 
P_Phi_A = Nc * Nf * Ns * Na * Phis_A 
P_Phi_G = (Nc**2 - 1) * Ns * Phis_G 

P_Phi = P_Phi_G + P_Phi_Q + P_Phi_A
P_tot = P_QP_Q + P_QP_A + P_Phi


df_P = pd.DataFrame(array([Trange, P_tot, P_Q_Q, P_Q_A, 
                           P_S_Q, P_S_A, 
                           P_Phi, P_Phi_Q, P_Phi_A, P_Phi_G]).transpose(), 
               columns=['T', 'Ptot', 'P_Q_Q', 'P_Q_A', 
                        'P_S_Q', 'P_S_A', 
                        'P_Phi', 'P_Phi_Q', 'P_Phi_A', 'P_Phi_G'])



# df_P.to_csv('pressure.csv')


df.attrs.update(
     dict(zip(['Ptot', 'P_Q_Q', 'P_Q_A',  'P_S_Q', 'P_S_A',  'P_Phi', 'P_Phi_Q', 'P_Phi_A', 'P_Phi_G'], 
            array([P_tot, P_Q_Q, P_Q_A,  P_S_Q, P_S_A,  P_Phi, P_Phi_Q, P_Phi_A, P_Phi_G])))
)

df_phi.attrs.update(
dict(zip(['Ptot', 'P_Q_Q', 'P_Q_A',  'P_S_Q', 'P_S_A',  'P_Phi', 'P_Phi_Q', 'P_Phi_A', 'P_Phi_G'], 
            array([P_tot, P_Q_Q, P_Q_A,  P_S_Q, P_S_A,  P_Phi, P_Phi_Q, P_Phi_A, P_Phi_G])))
            )

df.close()
df_phi.close()
# from matplotlib import pyplot as plt
# import matplotlib
# matplotlib.style.use('publication')

# lp, = plt.plot(Trange, P_Phi/Trange**4, label=r'$\frac{1}{ 2 }\Phi$')
# lQ, = plt.plot(Trange, P_QP_Q/Trange**4, label='quarks')
# plt.plot(Trange, 3*3*2*2*ps_free_Q/Trange**4, ls=':', c=lQ.get_c())
# lG, = plt.plot(Trange, P_QP_G/Trange**4, label='gluons')
# plt.plot(Trange, 8*2*ps_free_G/Trange**4, ls=':', c=lG.get_c(), label='free')
# plt.plot(Trange, (P_QP_Q + P_QP_G)/Trange**4, ls='-.', c='black')
# plt.plot(Trange, P_tot/Trange**4, c='black', label='total')
# plt.ylim(-0.5, 5)
# plt.legend(ncol=2, fontsize=14)

# lat = pd.read_csv(os.path.join(os.path.dirname(__file__), "PT.csv"))
# plt.plot(lat.x, lat.PT_lat, ls='none', marker='o')
# plt.axhline(0, lw=1, ls=':', c='black')
# plt.xlabel('T [GeV]')

# # try:
# #     plt.title(df.attrs['comment'])
# # except:
# #     pass
# plt.ylabel(r'P/T$^4$')


# plt.savefig(os.path.join(folder, 'PT.pdf'), bbox_inches='tight')
