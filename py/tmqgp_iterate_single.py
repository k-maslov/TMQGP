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
from scipy import signal
from scipy.integrate import quad
from syntax_sugar import END, pipe
from syntax_sugar import process_syntax as p
from syntax_sugar import thread_syntax as t
NTHR = 18

################## Reading the parameter file ######################

parser = argparse.ArgumentParser()
parser.add_argument('--mode', type=str, default='LO', choices=['LO', 'HI', 'XHI'])
parser.add_argument('T', type=float)
parser.add_argument('mQ', type=float)
parser.add_argument('mG', type=float)
parser.add_argument('G', type=float, default=10)
parser.add_argument('G1', type=float, default=10)
parser.add_argument('L',  type=float, default=0.2)
parser.add_argument('screen',  type=float, default=0.01)
parser.add_argument('suppress',  type=float, default=1)
parser.add_argument('--init', type=str, default='')
parser.add_argument('--folder', type=str, default='')
parser.add_argument('--save_iter', action='store_true')
parser.add_argument('-f', '--force', action='store_true')
parser.add_argument('--showtime', action='store_true')
parser.add_argument('--file', type=str, default='')
parser.add_argument('--max_iter', type=int, default=25)


args = parser.parse_args()
print(args)
mode = args.mode

T = args.T
mQ = args.mQ
mG = args.mG
G = args.G
G1 = args.G1
L = args.L
screen = args.screen

lmax = 1

suppress = args.suppress

if args.folder == '':
    out_folder = '.'
else:
    out_folder = args.folder

# comment = 'Both masses decrease, m_g 1.4 to 1.0\nscreen  0.9'

save_iterations = args.save_iter

if out_folder != '.':
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    else:
        if not args.force:
            raise ValueError(f"Not gonna overwrite folder {out_folder}")
        else:
            pass

####################################################################
######################## Setting up the parameters #################
if mode == 'XHI':
    erange = np.linspace(-5, 5, 401)
    qrange = np.linspace(0, 5, 101)   
elif mode == 'HI':
    erange = np.linspace(-5, 5, 201)
    qrange = np.linspace(0, 5, 51)    

elif mode == 'HI_Q':
    erange = np.linspace(-5, 5, 201)
    qrange = np.linspace(0, 5, 201)    
    
elif mode == 'LO':
    erange = np.linspace(-5, 5, 101)
    qrange = np.linspace(0, 5, 21)
else:
    raise ValueError

eps = 0.05

params = {'G' : G, 'L' : L, 'screen' : screen}
params1 = {'G' : G1, 'L' : L, 'screen' : screen}

params_QG = {'G' : G, 'L': suppress*L, 'screen': screen}
params_QG1 = {'G' : G1, 'L' : suppress*L, 'screen' : screen}

params_GG = {'G' : G, 'L': suppress**2 * L, 'screen': screen}
params_GG1 = {'G' : G1, 'L' : suppress**2 * L, 'screen' : screen}
   
params_rep = params.copy()
params_rep['G'] = -params_rep['G']
params_rep1 = params1.copy()
params_rep1['G'] = -params_rep1['G']
params_rep_QG = params_QG.copy()
params_rep_QG['G'] = -params_QG['G']
params_rep_QG1 = params_QG1.copy()
params_rep_QG1['G'] = -params_QG1['G']
params_rep_GG = params_GG.copy()
params_rep_GG['G'] = -params_GG['G']
params_rep_GG1 = params_GG1.copy()
params_rep_GG1['G'] = -params_GG1['G']

############################# Setting up the IO ##############

Tlabel = '%i'%(T*1e3)

fname = f'data_single_{Tlabel}.hdf5'
if args.file != '':
    fname = args.file

f = h5py.File(os.path.join(out_folder, fname), 'w')
f.attrs.update({
    'G' : G,
    'G1' : G1,
    "screen" : screen,
    'L' : L,
    'T' : T,
    'mQ' : mQ,
    'mG' : mG,
    'qrange' : qrange,
    'erange' : erange,
    'suppress' : suppress,
    'status' : 'RUN'
    # 'comment' : comment
})

########################### Iteration logic ###################

if args.init == '':
    quark_run = Particle(mQ, qrange, erange, eps=eps)
    gluon_run = Particle(mG, qrange, erange, eps=eps, stat='b', d=16, propagator=1)
else:
    print('Loading init from ' + args.init)
    f_init = h5py.File(args.init, 'r')
    quark_run = Particle(mQ, qrange, erange, eps=eps, Gtab=np.array(f_init['Q']['G']))
    gluon_run = Particle(mG, qrange, erange, eps=eps, stat='b',
                         d=16, propagator=1, Gtab=np.array(f_init['G']['G']))
    f_init.close()

print('Running G = ', G, 'T = ', T)

chss = []
pts = []

IMAGss = []
REALss = []

delta = 1

Nf = 3

n_iter = 0
current_iter = 0
thr = 1e-2

while abs(delta) > thr:
    channels_QQ = QuarkTM.ChannelGroup()
    channels_QG = QuarkTM.ChannelGroup()

    pss = [params, params1]
    pss_rep = [params_rep, params_rep1]

    pss_QG = [params_QG, params_QG1]
    pss_rep_QG = [params_rep_QG, params_rep_QG1]

    pss_GG = [params_GG, params_GG1]
    pss_rep_GG = [params_rep_GG, params_rep_GG1]
    
    
    channels_QQ.addChannel(
        QuarkTM.ChannelL('qq3', lmax, quark_run, quark_run, T, pss, ds=4, da=3, Fa=1/2, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qa1', lmax, quark_run, quark_run, T, pss, ds=4, da=1, Fa=1, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qq6', lmax, quark_run, quark_run, T, pss_rep, ds=4, da=6, Fa=1/4, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qa8', lmax, quark_run, quark_run, T, pss_rep, ds=4, da=8, Fa=1/8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg3', lmax, quark_run, gluon_run, T, pss_QG, ds=4, da=3, Fa=9./8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg6', lmax, quark_run, gluon_run, T, pss_QG, ds=4, da=6, Fa=3./8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg15', lmax, quark_run, gluon_run, T, pss_rep_QG, ds=4, da=15, Fa=3./8, )
    )

    channels_GQ = QuarkTM.ChannelGroup()
    channels_GG = QuarkTM.ChannelGroup()
# 
    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq3', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=3, Fa=9/8)
    )

    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq6', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=6, Fa=3/8)
    )

    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq15', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=15, Fa=3/8)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg1', lmax, gluon_run, gluon_run, T, pss_GG, ds=4, da=1, Fa=9/4)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg16', lmax, gluon_run, gluon_run, T, pss_GG, ds=4, da=16, Fa=9/8)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg27', lmax, gluon_run, gluon_run, T, pss_rep_GG, ds=4, da=27, Fa=3/4)
    )

    if n_iter == 0:
        for k, ch_l in (list(channels_QQ.channels.items()) + list(channels_QG.channels.items()) 
                    + list(channels_GQ.channels.items()) + list(channels_GG.channels.items())):
            for l in range(lmax + 1):
                ch = ch_l.chs[l]
                # print(k, ch.Nf)
                lbl = f'V/{k}/{l}'
                # print(lbl)
                f.create_dataset(lbl, data=ch.v(ch.qrange))
                f[f'V/{k}/{l}'].attrs.update({'ds': ch.ds, 'da' : ch.da, 'Fa' : ch.Fa})

    # exit()

    for chg in [channels_QQ, channels_QG, channels_GQ, channels_GG]:
        TM = chg.get_T()

    IMAGs = dict()
    REALs = dict()
    keys = ['QQ', 'QG', 'GQ', 'GG']
    funcs = [tm.sigma_ff_onshell, tm.sigma_fb_onshell, tm.sigma_bf_onshell, tm.sigma_bb_onshell]

    for key, func, channels in zip(keys, funcs, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
        TM_tot = channels.get_T()

        # lbl = f'TMtot/{key}'
        # print(lbl)
        # f.create_dataset(lbl, data=TM_tot)

        # break
        
        # plt.plot(erange, imag(TM_tot[:, 0]))
        ch = list(channels.channels.items())[0][1] #ch = list(channels.items())[0][1] # take any of the channels since SFs are the same
        iImTM_tot = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.imag(TM_tot)))
        eps1 = tm.Interpolator(qrange, ch.p_i.om0(qrange), 'cubic')
        eps2 = tm.Interpolator(qrange, ch.p_j.om0(qrange), 'cubic')
        Ntot = len(erange)*len(qrange)
        pairs = np.array([[[q, e] for e in erange] for q in qrange]).reshape(1, Ntot, 2)[0]
        
        # for q in tqdm.tqdm(qrange):
        
            # res = pipe(erange) | p[lambda z: func(z, q, ch.T, iImTM_tot, 
            #                                                 ch.p_j.R, eps1, eps2)]*(NTHR//1) | END
        start = timer()
        ress = np.array(pipe(pairs) | p[lambda z: func(z[1], z[0], ch.T, iImTM_tot, ch.p_j.R, eps1, eps2)]*(NTHR//1) | END)
        end = timer()
        if args.showtime:
            print('Channel ', key, 'time = ', end-start, ' s')
        ress = ress.reshape(len(qrange), len(erange)).transpose()

        # ress = np.array(ress).transpose()
        IMAGs[key] = ress

        ReSigmas = []

        for res in ress.transpose():
            iImSigma = tm.Interpolator(ch.erange, np.ascontiguousarray(res), 'cubic')
            ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in ch.erange]
            ReSigmas += [ReSigma]

        REALs[key] = np.array(ReSigmas).transpose()            

    IMAGss += [IMAGs]
    REALss += [REALs]
    
    ImS_Q = IMAGs['QQ'] + IMAGs['QG']
    ReS_Q = REALs['QQ'] + REALs['QG']

    om0_k = np.array([np.sqrt(mQ**2 + qrange**2) for e in quark_run.erange])
    arrE = np.array([quark_run.erange for q in quark_run.qrange]).transpose()

    G_Q_new = 1/(arrE - om0_k + 0*1j*quark_run.eps - (ReS_Q + 1j*ImS_Q))

    quark_new = Particle(mQ, qrange, erange, eps=quark_run.eps, Gtab=G_Q_new)
    
    quark_new.S = ReS_Q + 1j*ImS_Q
    
    ImS_G = IMAGs['GG'] + IMAGs['GQ']
    ReS_G = REALs['GG'] + REALs['GQ']
    
    # om0_k = np.array([gluon_run.om0(gluon_run.qrange) for e in gluon_run.erange])
    om0_k = np.array([np.sqrt(mG**2 + qrange**2) for e in gluon_run.erange])
    arrE = np.array([gluon_run.erange for q in gluon_run.qrange]).transpose()

#     G_G_new = 1/(arrE**2 - om0_k**2 + 2*1j*gluon_run.eps*arrE - (ReS_G + 1j*ImS_G))
    G_G_new = 1/(arrE - om0_k + 0*1j*gluon_run.eps - (ReS_G + 1j*ImS_G))
#     gluon_new = Particle(gluon_run.m, qrange, erange, eps=2e-2, stat='b', d=16, R=-2*imag(G_G_new))
    gluon_new = Particle(mG, qrange, erange, eps=gluon_run.eps, stat='b', d=16, Gtab=G_G_new)
    
    gluon_new.S = ReS_G + 1j*ImS_G


#         delta = 0
#         delta += np.sqrt(np.sum((quark_new.Rtab - quark_run.Rtab)**2)) / len(erange) / len(qrange)
# #         delta += sum(quark_new.Rtab / quark_run.Rtab) / len(erange) / len(qrange)
#         delta += np.sqrt(np.sum((gluon_new.Rtab - gluon_run.Rtab)**2)) / len(erange) / len(qrange)
#         delta /= 2
#         print(delta)  

    delta = 0
    wh = quark_run.Rtab != 0
    delta += np.max(np.abs(quark_new.Rtab[wh] - quark_run.Rtab[wh]))
    wh = gluon_run.Rtab != 0
    delta += np.max(np.abs(gluon_new.Rtab[wh] - gluon_run.Rtab[wh]))
    delta /= 2


        
    if abs(delta) < thr and n_iter == 0:
        quark_run.S = quark_new.S
        gluon_run.S = gluon_new.S
    
    chss += [[channels_QQ, channels_QG, channels_GQ, channels_GG]]
    pts += [[quark_run, gluon_run]]

    quark_bup = quark_run

    quark_run = quark_new
    gluon_run = gluon_new


    sumQ = np.trapz(quark_new.Rtab[:, 0], x=erange)
    sumG = np.trapz(gluon_new.Rtab[:, 0], x=erange)
    # print('delta = ', delta)
    # print('sum Q = ', sumQ)
    # print('sum G = ', sumG)

    print('T = %.3f iteration %i: delta = %f, sum Q = %f, sum G = %f'%(T, n_iter,
    delta, sumQ, sumG))

    if save_iterations:
        folder_iter = os.path.join(out_folder, Tlabel)
        if not os.path.exists(folder_iter):
            os.mkdir(folder_iter)
        fname_iter = folder_iter + '/iter_%i.hdf5'%n_iter
        print(fname_iter)
        f_iter = h5py.File(fname_iter, 'w')
        # f.create_dataset('%i/TM'%n_iter)
        kk = ['Q', 'Q', 'G', 'G']
        for _k, gr in zip(kk, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
            for k, chl in gr.channels.items():
                for l in range(chl.lmax + 1):
                    c = chl.chs[l]
                    f_iter.create_dataset(f'TM/{k}/{l}', data=c.TM)
                    f_iter.create_dataset(f'X/{k}/{l}', data=c.XS[0])

        for k in IMAGs.keys():
            f_iter.create_dataset(f'S/{k}', data=REALs[k] + 1j*IMAGs[k])
        f_iter.create_dataset(f'Q/G', data=quark_new.Gtab)
        f_iter.create_dataset(f'G/G', data=gluon_new.Gtab)
        f_iter.create_dataset(f'Q/S', data=quark_new.S)
        f_iter.create_dataset(f'G/S', data=gluon_new.S)
        for key, func, channels in zip(keys, funcs, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
            TM_tot = channels.get_T()

            lbl = f'TMtot/{key}'
            print(lbl)
            f_iter.create_dataset(lbl, data=TM_tot)
        f_iter.close()


    if n_iter > 1:
        if any([abs(sumQ - 1) > 1e-1, abs(sumG - 1) > 1e-1]):
            # break
            raise ValueError('Sum rule not fulfilled, STOP')
        

    n_iter += 1

    if n_iter > args.max_iter:
        print('\n########## max iter reached ############\n')
        break

IMAGs = IMAGss[-1]
REALs = REALss[-1]

kk = ['Q', 'Q', 'G', 'G']
for _k, gr in zip(kk, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
    for k, chl in gr.channels.items():
        for l in range(chl.lmax + 1):
            c = chl.chs[l]
            f.create_dataset(f'TM/{k}/{l}', data=c.TM)
            f.create_dataset(f'X/{k}/{l}', data=c.XS[0])


for k in IMAGs.keys():
    f.create_dataset(f'S/{k}', data=REALs[k] + 1j*IMAGs[k])

f.create_dataset(f'Q/R', data=pts[-1][0].Rtab)
f.create_dataset(f'G/R', data=pts[-1][1].Rtab)
# print(pts[-1][0].Gtab)

f.create_dataset(f'Q/G', data=pts[-1][0].Gtab)
f.create_dataset(f'G/G', data=pts[-1][1].Gtab)

f.create_dataset(f'Q/S', data=pts[-1][0].S)
f.create_dataset(f'G/S', data=pts[-1][1].S)

f.attrs.update({'status' : 'DONE'})

f.close()



########################## Thermodynamics ##########################