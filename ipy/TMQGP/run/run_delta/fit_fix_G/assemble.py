import h5py
import os
import pandas as pd
import numpy as np
import QuarkTM

# iLat = interp1d(lat.x, lat.PT_lat, kind='cubic')
lat = pd.read_csv(os.path.join(os.path.dirname(QuarkTM.__file__), "PT.csv"))
Trange_fit = lat.x[(lat.x > 0.16) & (lat.x < 0.4)].values[::2]

df_out = h5py.File('result.hdf5', 'w')
th_out = h5py.File('thermo.hdf5', 'w')

mQs = []
mGs = []
for T in Trange_fit:
    mQ = float(np.loadtxt('sol_%.3f.dat'%T))
    fname = 'fit_%.3f_%.16e.h5py'%(T, mQ)

    df = h5py.File(fname, 'r')
    df_th = h5py.File('th_' + fname, 'r')

    Tkey = str(int(1e3*T))
    # df_out.create_group(f'{Tkey}/')
    df.copy(df, df_out, Tkey)
    df_th.copy(df_th, th_out, Tkey)
    # df_out.create_dataset(f'{Tkey}/', data=df)
    mQs += [df.attrs['mQ']]
    mGs += [df.attrs['mG']]

    df.close()
    df_th.close()

print(mQs)
print(mGs)
print(Trange_fit)

df_out.attrs.update({'mQs' : mQs, 'mGs' : mGs, 'Trange' : Trange_fit})

keys = ['P_Phi', 'P_Phi_G', 'P_Phi_Q', 'P_Q_G', 'P_Q_Q', 'P_S_G', 'P_S_Q', 'Ptot']

d_p = {k : [float(d.attrs[k]) for _, d in df_out.items()] for k in keys}

df_p = pd.DataFrame(d_p)

print(df_p)



df_out.close()

from matplotlib import pyplot as plt
import matplotlib
matplotlib.style.use('publication')
Trange = Trange_fit

P_QP_Q = df_p.P_Q_Q + df_p.P_S_Q
P_QP_G = df_p.P_Q_G + df_p.P_S_G

lp, = plt.plot(Trange, df_p.P_Phi/Trange**4, label=r'$\frac{1}{ 2 }\Phi$')
lQ, = plt.plot(Trange, P_QP_Q/Trange**4, label='quarks')
# plt.plot(Trange, 3*3*2*2*ps_free_Q/Trange**4, ls=':', c=lQ.get_c())
lG, = plt.plot(Trange, P_QP_G/Trange**4, label='gluons')
# plt.plot(Trange, 8*2*ps_free_G/Trange**4, ls=':', c=lG.get_c(), label='free')
plt.plot(Trange, (P_QP_Q + P_QP_G)/Trange**4, ls='-.', c='black')
plt.plot(Trange, df_p.Ptot/Trange**4, c='black', label='total')
plt.ylim(-0.5, 5)
plt.legend(ncol=2, fontsize=14)

# lat = pd.read_csv(os.path.join(os.path.dirname(Quark), "PT.csv"))
plt.plot(lat.x, lat.PT_lat, ls='none', marker='o')
plt.axhline(0, lw=1, ls=':', c='black')
plt.xlabel('T [GeV]')

# try:
#     plt.title(df.attrs['comment'])
# except:
#     pass
plt.ylabel(r'P/T$^4$')


plt.savefig('PT.pdf', bbox_inches='tight')



