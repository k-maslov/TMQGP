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
mode = 'LO'
Trange = [0.16, 0.2, 0.3, 0.4, 0.5]
mQ = 0.6
mG = 1.8
G = 14
G1 = 14.5
L = .2
screen = .01

lmax = 1

suppress = 1

out_folder = './output/TestSuppressIter_' +mode+'_G=(%.2f,%.2f)L=%.3fMQ=%.2fMG=%.2fscreen=%.3fsuppress=%.2f/'%(G, G1, L, mQ, mG, screen, suppress)

save_iterations = 0

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

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
quark_run = Particle(mQ, qrange, erange, eps=eps)
gluon_run = Particle(mG, qrange, erange, eps=eps, stat='b', d=16, propagator=1)

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

f = h5py.File(out_folder + 'data.hdf5', 'w')
f.attrs.update({
    'G' : G,
    'G1' : G1,
    "screen" : screen,
    'L' : L,
    'Trange' : Trange,
    'mQ' : mQ,
    'mG' : mG,
    'qrange' : qrange,
    'erange' : erange,
    'suppress' : suppress
})

########################### Iteration logic ###################

for T in Trange[:]:    
    print('Running G = ', G, 'T = ', T)
    Tlabel = '%i'%(T*1e3)
    chss = []
    pts = []

    IMAGss = []
    REALss = []

    delta = 1
    
    Nf = 3

    n_iter = 0
    current_iter = 0
    
    while abs(delta) > 1e-4:


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
    
                    lbl = f'{Tlabel}/V/{k}/{l}'
                    
                    f.create_dataset(lbl, data=ch.v(ch.qrange))
                    f[f'{Tlabel}/V/{k}/{l}'].attrs.update({'ds': ch.ds, 'da' : ch.da, 'Fa' : ch.Fa})
                    
                    np.savetxt(out_folder + f'{k}_vq_l=%i_T=%.3f.dat'%(l, ch.T), ch.v(ch.qrange))
                np.savetxt(out_folder + f'{k}_weights', np.array([ch.ds, ch.da, ch.Fa]))

        for chg in [channels_QQ, channels_QG, channels_GQ, channels_GG]:
            TM = chg.get_T()

        IMAGs = dict()
        REALs = dict()
        keys = ['QQ', 'QG', 'GQ', 'GG']
        funcs = [tm.sigma_ff_onshell, tm.sigma_fb_onshell, tm.sigma_bf_onshell, tm.sigma_bb_onshell]

        for key, func, channels in zip(keys, funcs, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
            TM_tot = channels.get_T()

            ch = list(channels.channels.items())[0][1] #ch = list(channels.items())[0][1] # take any of the channels since SFs are the same
            iImTM_tot = tm.Interpolator2D(qrange, erange, np.ascontiguousarray(np.imag(TM_tot)))
            eps1 = tm.Interpolator(qrange, ch.p_i.om0(qrange), 'cubic')
            eps2 = tm.Interpolator(qrange, ch.p_j.om0(qrange), 'cubic')
            Ntot = len(erange)*len(qrange)
            pairs = np.array([[[q, e] for e in erange] for q in qrange]).reshape(1, Ntot, 2)[0]

            start = timer()
            ress = np.array(pipe(pairs) | p[lambda z: func(z[1], z[0], ch.T, iImTM_tot, ch.p_j.R, eps1, eps2)]*(NTHR//1) | END)
            end = timer()
            print('Channel ', key, 'time = ', end-start, ' s')
            ress = ress.reshape(len(qrange), len(erange)).transpose()

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

        om0_k = np.array([quark_run.om0(quark_run.qrange) for e in quark_run.erange])
        arrE = np.array([quark_run.erange for q in quark_run.qrange]).transpose()

        G_Q_new = 1/(arrE - om0_k + 0*1j*quark_run.eps - (ReS_Q + 1j*ImS_Q))

        quark_new = Particle(quark_run.m, qrange, erange, eps=quark_run.eps, Gtab=G_Q_new)
        
        quark_new.S = ReS_Q + 1j*ImS_Q
        
        ImS_G = IMAGs['GG'] + IMAGs['GQ']
        ReS_G = REALs['GG'] + REALs['GQ']
        
        om0_k = np.array([gluon_run.om0(gluon_run.qrange) for e in gluon_run.erange])
        arrE = np.array([gluon_run.erange for q in gluon_run.qrange]).transpose()

        G_G_new = 1/(arrE - om0_k + 0*1j*gluon_run.eps - (ReS_G + 1j*ImS_G))

        gluon_new = Particle(gluon_run.m, qrange, erange, eps=gluon_run.eps, stat='b', d=16, Gtab=G_G_new)
        
        gluon_new.S = ReS_G + 1j*ImS_G
        
        chss += [[channels_QQ, channels_QG, channels_GQ, channels_GG]]
        pts += [[quark_run, gluon_run]]

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
        print('delta = ', delta)

        quark_bup = quark_run

        quark_run = quark_new
        gluon_run = gluon_new


        sumQ = np.trapz(quark_new.Rtab[:, 0], x=erange)
        sumG = np.trapz(gluon_new.Rtab[:, 0], x=erange)
        print('sum Q = ', sumQ)
        print('sum G = ', sumG)

        if any([abs(sumQ - 1) > 1e-1, abs(sumG - 1) > 1e-1]):
            raise ValueError('Sum rule not fulfilled, STOP')

        if save_iterations:
            folder_iter = out_folder + Tlabel
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
            f_iter.close()

        n_iter += 1

    IMAGs = IMAGss[-1]
    REALs = REALss[-1]

    kk = ['Q', 'Q', 'G', 'G']
    for _k, gr in zip(kk, [channels_QQ, channels_QG, channels_GQ, channels_GG]):
        for k, chl in gr.channels.items():
            for l in range(chl.lmax + 1):
                c = chl.chs[l]
                f.create_dataset(f'{Tlabel}/TM/{k}/{l}', data=c.TM)
                f.create_dataset(f'{Tlabel}/X/{k}/{l}', data=c.XS[0])
                np.savetxt(out_folder + _k + '_TM_'+k+'_l=%i'%l + '_T=%.3f.dat'%T, c.TM)
                np.savetxt(out_folder + _k + '_X_'+k+'_l=%i'%l + '_T=%.3f.dat'%T, c.XS[0])


    for k in IMAGs.keys():
        f.create_dataset(f'{Tlabel}/S/{k}', data=REALs[k] + 1j*IMAGs[k])
        np.savetxt(out_folder + k+'_S_T=%.3f.dat'%T, 
                  REALs[k] + 1j*IMAGs[k])

    f.create_dataset(f'{Tlabel}/Q/R', data=pts[-1][0].Rtab)
    np.savetxt(out_folder + 'Q_Rho_T=%.3f.dat'%T, pts[-1][0].Rtab)
    f.create_dataset(f'{Tlabel}/G/R', data=pts[-1][1].Rtab)
    np.savetxt(out_folder + 'G_Rho_T=%.3f.dat'%T, pts[-1][1].Rtab)

    f.create_dataset(f'{Tlabel}/Q/G', data=pts[-1][0].Gtab)
    np.savetxt(out_folder + 'Q_G_T=%.3f.dat'%T, pts[-1][0].Gtab)
    f.create_dataset(f'{Tlabel}/G/G', data=pts[-1][1].Gtab)
    np.savetxt(out_folder + 'G_G_T=%.3f.dat'%T, pts[-1][1].Gtab)

    f.create_dataset(f'{Tlabel}/Q/S', data=pts[-1][0].S)
    np.savetxt(out_folder + 'Q_S_T=%.3f.dat'%T, pts[-1][0].S)
    f.create_dataset(f'{Tlabel}/G/S', data=pts[-1][1].S)
    np.savetxt(out_folder + 'G_S_T=%.3f.dat'%T, pts[-1][1].S)

    np.savetxt(out_folder + 'erange_T=%.3f.dat'%T, erange)
    np.savetxt(out_folder + 'qrange_T=%.3f.dat'%T, qrange)
    
    # pickle.dump(chss, open(out_folder + "chss_T=%.3f.p"%T, "wb" ))

    # pickle.dump(pts, open(out_folder + "pts_T=%.3f.p"%T, "wb" ))


f.close()