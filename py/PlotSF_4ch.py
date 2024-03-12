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
matplotlib.style.use('publication')
import PyPDF2

parser = argparse.ArgumentParser()

parser.add_argument('folder')
args = parser.parse_args()

folder = args.folder
df = h5py.File(os.path.join(args.folder, 'data.hdf5'), 'r')

nplots = len(list(df.keys()))

keys = list(df.keys())

erange = df[keys[0]].attrs['erange']
trange = df.attrs['Trange']

fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

for i, key in enumerate(sorted(list(df.keys()))):
    ax[0][i].plot(erange, imag(df[key]['TM']['qa1']['0'][:, 0]))
    ax[1][i].plot(erange, imag(df[key]['TM']['qq3']['0'][:, 0]))
    ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

for _ in ax.flatten():
    _.set_xlim(0, 2.99)

for _ in ax[1]:
    _.set_xlabel(r'E [GeV]')

ax[0][0].set_ylabel(r'$T_{\bar q q_1}$ [GeV$^{-2}$]')
ax[1][0].set_ylabel(r'$T_{q q_3}$ [GeV$^{-2}$]')

plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig(os.path.join(folder, 'T_Q_l0.pdf'), bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

for i, key in enumerate(sorted(list(df.keys()))):
    ax[0][i].plot(erange, imag(df[key]['TM']['aq1']['0'][:, 0]))
    ax[1][i].plot(erange, imag(df[key]['TM']['aa3']['0'][:, 0]))
    ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

for _ in ax.flatten():
    _.set_xlim(0, 2.99)

for _ in ax[1]:
    _.set_xlabel(r'E [GeV]')

ax[0][0].set_ylabel(r'$T_{\bar q q_1}$ [GeV$^{-2}$]')
ax[1][0].set_ylabel(r'$T_{\bar q \bar q_3}$ [GeV$^{-2}$]')

plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig(os.path.join(folder, 'T_A_l0.pdf'), bbox_inches='tight')
plt.close()


fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

for i, key in enumerate(sorted(list(df.keys()))):
    ax[0][i].plot(erange, imag(df[key]['TM']['qa1']['1'][:, 10]))
    ax[1][i].plot(erange, imag(df[key]['TM']['qq3']['1'][:, 10]))
    ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

for _ in ax.flatten():
    _.set_xlim(0, 2.99)

for _ in ax[1]:
    _.set_xlabel(r'E [GeV]')

ax[0][0].set_ylabel(r'$T_{\bar q q_1}$ [GeV$^{-2}$]')
ax[1][0].set_ylabel(r'$T_{q q_3}$ [GeV$^{-2}$]')

plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig(os.path.join(folder, 'T_Q_l1.pdf'), bbox_inches='tight')
plt.close()

# fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

# for i, key in enumerate(sorted(list(df.keys()))):
#     # ax[0][i].plot(erange, imag(df[key]['TM']['gg1']['0'][:, 0]))
#     # ax[1][i].plot(erange, imag(df[key]['TM']['gg27']['0'][:, 0]))
#     ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

# for _ in ax.flatten():
#     _.set_xlim(0, 4.99)

# for _ in ax[1]:
#     _.set_xlabel(r'E [GeV]')

# ax[0][0].set_ylabel(r'$T_{gg_1}$ [GeV$^{-2}$]')
# ax[1][0].set_ylabel(r'$T_{gg_{27}}$ [GeV$^{-2}$]')

# plt.subplots_adjust(hspace=0, wspace=0)

# plt.savefig(os.path.join(folder, 'T_G_l0.pdf'), bbox_inches='tight')

# plt.close()
fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

for i, key in enumerate(sorted(list(df.keys()))):
    ax[0][i].plot(erange, imag(df[key]['Q']['S'][:, 0]))
    ax[0][i].plot(erange, real(df[key]['Q']['S'][:, 0]), ls='--')

    ax[1][i].plot(erange, (df[key]['Q']['R'][:, 0]))
    ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

for _ in ax.flatten():
    _.set_xlim(-2, 3.99)

for _ in ax[1]:
    _.set_xlabel(r'E [GeV]')

ax[0][0].set_ylabel(r'Im$\Sigma_q$ [GeV]')
ax[1][0].set_ylabel(r'$\rho_q$ [GeV$^{-1}$]')

plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig(os.path.join(folder, 'Rho_Q.pdf'), bbox_inches='tight')

plt.close()

fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

for i, key in enumerate(sorted(list(df.keys()))):
    ax[0][i].plot(erange, imag(df[key]['A']['S'][:, 0]))
    ax[1][i].plot(erange, (df[key]['A']['R'][:, 0]))
    ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

for _ in ax.flatten():
    _.set_xlim(-2, 3.99)

for _ in ax[1]:
    _.set_xlabel(r'E [GeV]')

ax[0][0].set_ylabel(r'Im$\Sigma_q$ [GeV]')
ax[1][0].set_ylabel(r'$\rho_q$ [GeV$^{-1}$]')

plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig(os.path.join(folder, 'Rho_A.pdf'), bbox_inches='tight')

plt.close()

# fig, ax = plt.subplots(2, nplots, figsize=(4*nplots, 8), sharey='row', sharex='all')

# for i, key in enumerate(sorted(list(df.keys()))):
#     ax[0][i].plot(erange, imag(df[key]['G']['S'][:, 0]))
#     ax[1][i].plot(erange, (df[key]['G']['R'][:, 0]))
#     ax[0][i].set_title('T = %.2f'%trange[i])
# for i in irange:

# for _ in ax.flatten():
#     _.set_xlim(-2, 3.99)

# for _ in ax[1]:
#     _.set_xlabel(r'E [GeV]')

# ax[0][0].set_ylabel(r'Im$\Sigma_g$ [GeV]')
# ax[1][0].set_ylabel(r'$\rho_g$ [GeV$^{-1}$]')

# plt.subplots_adjust(hspace=0, wspace=0)

# plt.savefig(os.path.join(folder, 'Rho_G.pdf'), bbox_inches='tight')

# plt.close()


resonances = ['qa1', 'qq3']

mress = []

for r in resonances:
    mres = []
    for key in df.keys():
        tm = imag(df[key]['TM'][r]['0'][:, 0])
        # plt.plot(erange, tm)
        # print(erange[np.argmin(tm)])
        mres += [erange[np.argmin(tm)]]
    mress += [mres]

fig, ax = plt.subplots(1, 2, figsize=(12, 6))


ax[0].plot(trange, df.attrs['mQs'], label='quark')
ax[0].plot(trange, df.attrs['mGs'], label='gluon')
ax[0].legend()
for mres, r in zip(mress, resonances):
    ax[1].plot(trange, mres, label=r)

plt.legend(fontsize=14)
plt.ylim(0, 3)
ax[0].set_ylim(0, 2)
ax[0].set_ylabel('m [GeV]')
plt.xlabel('T [GeV]')
plt.ylabel(r'$m_R$ [GeV]')

plt.savefig(os.path.join(folder, 'mres.pdf'), bbox_inches='tight')

d = dict(df.attrs)

# del d['erange']
# del d['qrange']

fig = plt.figure()
fig.text(0, 0, str(d))
plt.savefig(os.path.join(folder, 'caption.pdf'), bbox_inches='tight')

merger = PyPDF2.PdfWriter()

for f in ['caption.pdf', 'mres.pdf', 'Rho_Q.pdf', 'Rho_A.pdf', 'T_Q_l0.pdf', 'T_A_l0.pdf']:
    pdf = os.path.join(folder, f)
    merger.append(pdf)
merger.write(os.path.join(folder, 'report.pdf'))
merger.close()