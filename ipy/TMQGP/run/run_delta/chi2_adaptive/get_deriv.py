import h5py

from numpy import *
from matplotlib import pyplot as plt

Ps = []
dmu = 0.2

mus_range = array([0, 0.2, 0.4, 0.6])

for muscale in [0, 0.2, 0.4, 0.6]:
    print (muscale)
    df = h5py.File('%.2f/data.hdf5'%muscale, 'r')
    Ps += [array(df.attrs['Ptot']) / array(df.attrs['Trange'])**4]

Trange = df.attrs["Trange"]

# f = (mus_range - 1)**3
#
# d2f = (2*f[0] - 5*f[1] + 4*f[2] - f[3]) / dmu**2
# print(d2f)

chi2 = (2*Ps[0] - 5*Ps[1] + 4*Ps[2] - Ps[3]) / dmu**2 / 9

savetxt(
    'chi2.txt',
    array([df.attrs['Trange'], chi2]).transpose()
)

plt.plot(df.attrs['Trange'], chi2)

plt.savefig('chi2.pdf', bbox_inches='tight')

plt.close()

for P in Ps:
    plt.plot(Trange, P, marker='.')

plt.savefig('Ps.pdf', bbox_inches='tight')

plt.close()

for P in Ps[1:]:
    plt.plot(Trange, P/Ps[0])

plt.savefig('Ps_ratio.pdf', bbox_inches='tight')
