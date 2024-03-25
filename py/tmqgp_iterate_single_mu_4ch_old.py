import argparse
import pickle
import os
import h5py
import matplotlib
from timeit import default_timer as timer

import TMQGP as tm
import QuarkTM
from QuarkTM import Channel, Particle
from matplotlib import pyplot as plt
import tqdm
import numpy as np
from scipy import signal
from scipy.integrate import quad
from syntax_sugar import END, pipe
from syntax_sugar import process_syntax as p
from syntax_sugar import thread_syntax as t

import warnings
warnings.filterwarnings('ignore')
NTHR = 18

################## Reading the parameter file ######################

parser = argparse.ArgumentParser()
parser.add_argument('--mode', type=str, default='LO', choices=['LO', 'HI', 'XHI'])
parser.add_argument('mu', type=float)
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

mu = args.mu
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

save_iterations = 1 #args.save_iter

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
    qrange = np.linspace(0, 5, 201)   
elif mode == 'HI':
    erange = np.linspace(-5, 5, 201)
    qrange = np.linspace(0, 5, 51)

elif mode == 'HI_Q':
    erange = np.linspace(-5, 5, 201)
    qrange = np.linspace(0, 5, 201)    
    
elif mode == 'LO':
    erange = np.linspace(-5, 5, 101)
    qrange = np.linspace(0, 5, 51)
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
    'status' : 'RUN',
    'mu' : mu
    # 'comment' : comment
})

########################### Iteration logic ###################

if args.init == '':
    quark_run = Particle(mQ, qrange, erange, eps=eps, mu=mu)
    aquark_run = Particle(mQ, qrange, erange, eps=eps, mu=-mu)
    # gluon_run = Particle(mG, qrange, erange, eps=eps, stat='b', d=16, propagator=1)
else:
    print('Loading init from ' + args.init)
    f_init = h5py.File(args.init, 'r')
    quark_run = Particle(mQ, qrange, erange, eps=eps, Gtab=np.array(f_init['Q']['G']), mu=mu)
    aquark_run = Particle(mQ, qrange, erange, eps=eps, Gtab=np.array(f_init['A']['G']), mu=-mu)
    # gluon_run = Particle(mG, qrange, erange, eps=eps, stat='b',
    #                      d=16, propagator=1, Gtab=np.array(f_init['G']['G']))
    f_init.close()


# plt.plot(qrange, quark_run.peaks)
# plt.show()

# plt.plot(erange, [quark_run.R(0, e) for e in erange])
# plt.show()


print('Running G = ', G, 'T = ', T)

chss = []
pts = []

IMAGss = []
REALss = []

delta = 1

n_iter = 0
current_iter = 0
thr = 1e-4

erange2b = np.linspace(0, 10, 1501)
qrange2b = np.linspace(0, 10, 501)

def getReKK(erange, ims):
    iIm = tm.Interpolator(erange, ims, 'cubic')
    
    res = np.array([tm.ReSigmaKK(e, iIm, Lambda=10) for e in erange])
    return res


chps = {'G2_mode' : 0, 'expand' : 0}

while abs(delta) > thr:
    channels_QQ = QuarkTM.ChannelGroup(mu0=False)
    # channels_QG = QuarkTM.ChannelGroup(mu0=False)
    channels_QA = QuarkTM.ChannelGroup(mu0=False)

    pss = [params, params1]
    pss_rep = [params_rep, params_rep1]

    pss_QG = [params_QG, params_QG1]
    pss_rep_QG = [params_rep_QG, params_rep_QG1]

    pss_GG = [params_GG, params_GG1]
    pss_rep_GG = [params_rep_GG, params_rep_GG1]



    ##### Parallel calculation of G_QQ
    Ntot = len(qrange2b) * len(erange2b)
    pairs = np.array([[[q, e] for e in erange2b] for q in qrange2b]).reshape(1, Ntot, 2)[0]

    

    def getIm(q, pt1, pt2):
        return np.array([-np.pi*tm.G2_conv_ff(e, q, T, pt1.R, pt2.R) for e in erange2b])


    # print('Calculating Im')
    # ImG2_QQ = np.array([
    #     [-np.pi*tm.G2_conv_ff(e, q, T, quark_run.R, quark_run.R, Lambda=10) for e in erange2b]
    # for q in (qrange2b)]).transpose()

    
    # ImG2_QQ = np.array(pipe(pairs) | p[lambda z: -np.pi*tm.G2_conv_ff(z[1], z[0], T, quark_run.R, quark_run.R, Lambda=10)]*(NTHR//1) | END)
    # start = timer()
    # ImG2_QQ = np.array(pipe(qrange2b) | p[lambda z: getIm(z, quark_run, quark_run)]*(NTHR//1) | END)
    # ImG2_QQ = ImG2_QQ.transpose()
    # # ImG2_QQ = ImG2_QQ.reshape(len(qrange2b), len(erange2b)).transpose()
    # # print('Done')

    # # print('Calculating Re')

    # ReG2_QQ = np.array(pipe(ImG2_QQ.transpose()) | p[lambda z: getReKK(erange2b, z)]*(NTHR//1) | END)
    # ReG2_QQ = ReG2_QQ.transpose()
    # # print('done')

    # end = timer()

    # if args.showtime:
    #     print('Channel QQ time = ', end-start, ' s')


    # # ch_test = QuarkTM.ChannelL('qq3', lmax, quark_run, quark_run, T, pss, ds=4, da=3, Fa=1/2, mu=mu)
    # plt.plot(erange2b, ImG2_QQ[:, 0])
    # plt.plot(erange2b, np.imag(ch_test.G2[:, 0]))

    # plt.show()

    # # exit()


    channels_QQ.addChannel(
        QuarkTM.ChannelL('qq3', lmax, quark_run, quark_run, T, pss, ds=4, da=3, Fa=1/2, mu=mu, **chps)
    )

    f.attrs.update({'erange2b': channels_QQ['qq3'].chs[0].erange, 
                    'qrange2b': channels_QQ['qq3'].chs[0].qrange})

    ch = channels_QQ['qq3'].chs[0]

    # plt.cla()
    # # plt.plot(ch.erange2b, [quark_run.iImG(0, e) for e in ch.erange2b])
    # # plt.plot(ch.erange2b, np.imag(channels_QQ['qq3'].chs[0].TM[:, 0]))
    # plt.plot(ch.erange2b, np.imag(channels_QQ['qq3'].chs[0].G2[:, 0]))
    # plt.plot(ch.erange2b, np.real(channels_QQ['qq3'].chs[0].G2[:, 0]))
    # plt.show()

    # channels_QQ.addChannel(
    #     QuarkTM.ChannelL('qq6', lmax, quark_run, quark_run, T, pss_rep, ds=4, da=6, Fa=1/4, mu=mu)
    # )

    # channels_QA.addChannel(
    #     QuarkTM.ChannelL('qa8', lmax, quark_run, aquark_run, T, pss_rep, ds=4, da=8, Fa=1/8, mu=mu)
    # )

    # start = timer()
    # ImG2_QA = np.array(pipe(qrange2b) | p[lambda z: getIm(z, quark_run, aquark_run)]*(NTHR//1) | END)
    # ImG2_QA = ImG2_QA.transpose()
    # # ImG2_QQ = ImG2_QQ.reshape(len(qrange2b), len(erange2b)).transpose()
    # # print('Done')

    # # print('Calculating Re')

    # ReG2_QA = np.array(pipe(ImG2_QQ.transpose()) | p[lambda z: getReKK(erange2b, z)]*(NTHR//1) | END)
    # ReG2_QA = ReG2_QA.transpose()
    # print('done')

    end = timer()

    if args.showtime:
        print('Channel QA time = ', end-start, ' s')

    channels_QA.addChannel(
        QuarkTM.ChannelL('qa1', lmax, quark_run, aquark_run, T, pss, ds=4, da=1, Fa=1, mu=mu, G2=ReG2_QA + 1j*ImG2_QA)
    )

    # channels_QG.addChannel(
    #     QuarkTM.ChannelL('qg3', lmax, quark_run, gluon_run, T, pss_QG, ds=4, da=3, Fa=9./8, mu=mu)
    # )

    # channels_QG.addChannel(
    #     QuarkTM.ChannelL('qg6', lmax, quark_run, gluon_run, T, pss_QG, ds=4, da=6, Fa=3./8, mu=mu)
    # )

    # channels_QG.addChannel(
    #     QuarkTM.ChannelL('qg15', lmax, quark_run, gluon_run, T, pss_rep_QG, ds=4, da=15, Fa=3./8, mu=mu)
    # )


    channels_AA = QuarkTM.ChannelGroup(mu0=False)
    # channels_AG = QuarkTM.ChannelGroup(mu0=False)
    channels_AQ = QuarkTM.ChannelGroup(mu0=False)


    start = timer()
    ImG2_AA = np.array(pipe(qrange2b) | p[lambda z: getIm(z, aquark_run, aquark_run)]*(NTHR//1) | END)
    ImG2_AA = ImG2_AA.transpose()
    # ImG2_QQ = ImG2_QQ.reshape(len(qrange2b), len(erange2b)).transpose()
    # print('Done')

    # print('Calculating Re')

    ReG2_AA = np.array(pipe(ImG2_QQ.transpose()) | p[lambda z: getReKK(erange2b, z)]*(NTHR//1) | END)
    ReG2_AA = ReG2_AA.transpose()
    # print('done')

    end = timer()

    if args.showtime:
        print('Channel AA time = ', end-start, ' s')


    channels_AA.addChannel(
        QuarkTM.ChannelL('aa3', lmax, aquark_run, aquark_run, T, pss, ds=4, da=3, Fa=1/2, mu=mu, G2=ReG2_AA + 1j*ImG2_AA)
    )

    # channels_AA.addChannel(
    #     QuarkTM.ChannelL('aa6', lmax, aquark_run, aquark_run, T, pss_rep, ds=4, da=6, Fa=1/4, mu=mu)
    # )

    # channels_AQ.addChannel(
    #     QuarkTM.ChannelL('aq8', lmax, aquark_run, quark_run, T, pss_rep, ds=4, da=8, Fa=1/8, mu=mu)
    # )

    channels_AQ.addChannel(
        QuarkTM.ChannelL('aq1', lmax, aquark_run, quark_run, T, pss, ds=4, da=1, Fa=1, mu=mu, G2=ReG2_QA + 1j*ImG2_QA)
    )

    # channels_AG.addChannel(
    #     QuarkTM.ChannelL('ag3', lmax, aquark_run, gluon_run, T, pss_QG, ds=4, da=3, Fa=9./8, mu=mu)
    # )

    # channels_AG.addChannel(
    #     QuarkTM.ChannelL('ag6', lmax, aquark_run, gluon_run, T, pss_QG, ds=4, da=6, Fa=3./8, mu=mu)
    # )

    # channels_AG.addChannel(
    #     QuarkTM.ChannelL('ag15', lmax, aquark_run, gluon_run, T, pss_rep_QG, ds=4, da=15, Fa=3./8, mu=mu)
    # )


#     channels_GQ = QuarkTM.ChannelGroup(mu0=False)
#     channels_GA = QuarkTM.ChannelGroup(mu0=False)
#     channels_GG = QuarkTM.ChannelGroup(mu0=False)
# # 
#     channels_GQ.addChannel(
#         QuarkTM.ChannelL('gq3', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=3, Fa=9/8, mu=mu)
#     )

#     channels_GQ.addChannel(
#         QuarkTM.ChannelL('gq6', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=6, Fa=3/8, mu=mu)
#     )

#     channels_GQ.addChannel(
#         QuarkTM.ChannelL('gq15', lmax, gluon_run, quark_run, T, pss_QG, ds=4, da=15, Fa=3/8, mu=mu)
#     )

#     channels_GA.addChannel(
#         QuarkTM.ChannelL('ga3', lmax, gluon_run, aquark_run, T, pss_QG, ds=4, da=3, Fa=9/8, mu=mu)
#     )

#     channels_GA.addChannel(
#         QuarkTM.ChannelL('ga6', lmax, gluon_run, aquark_run, T, pss_QG, ds=4, da=6, Fa=3/8, mu=mu)
#     )

#     channels_GA.addChannel(
#         QuarkTM.ChannelL('ga15', lmax, gluon_run, aquark_run, T, pss_QG, ds=4, da=15, Fa=3/8, mu=mu)
#     )

#     channels_GG.addChannel(
#         QuarkTM.ChannelL('gg1', lmax, gluon_run, gluon_run, T, pss_GG, ds=4, da=1, Fa=9/4, mu=mu)
#     )

#     channels_GG.addChannel(
#         QuarkTM.ChannelL('gg16', lmax, gluon_run, gluon_run, T, pss_GG, ds=4, da=16, Fa=9/8, mu=mu)
#     )

#     channels_GG.addChannel(
#         QuarkTM.ChannelL('gg27', lmax, gluon_run, gluon_run, T, pss_rep_GG, ds=4, da=27, Fa=3/4, mu=mu)
#     )
    

    erange_dense = np.linspace(-5, 5, 1001)
    # plt.plot(erange_dense, [quark_run.R(0, e) for e in erange_dense])
    # plt.show()

    if n_iter == 0:
        for k, ch_l in (list(channels_QQ.channels.items()) + 
                    list(channels_QA.channels.items()) + list(channels_AA.channels.items()) + 
                    list(channels_AQ.channels.items())):
            for l in range(lmax + 1):
                ch = ch_l.chs[l]
                # print(k, l, ch.Nf)
                lbl = f'V/{k}/{l}'
                # print(lbl)
                f.create_dataset(lbl, data=ch.v(ch.qrange))
                f[f'V/{k}/{l}'].attrs.update({'ds': ch.ds, 'da' : ch.da, 'Fa' : ch.Fa})
    

    # exit()
    for chg in [channels_QQ, channels_QA, channels_AQ, channels_AA]:
        TM = chg.get_T()

    ch = channels_QA['qa1'].chs[0]

    # plt.plot(ch.erange, np.imag(ch.TM[:, 0]))
    # plt.show()

    IMAGs = dict()
    REALs = dict()
    keys = ['QQ', 'QA', 
            'AQ', 'AA']

    funcs = [tm.sigma_ff_onshell, tm.sigma_ff_onshell, 
             tm.sigma_ff_onshell, tm.sigma_ff_onshell]

    chg_list = [channels_QQ, channels_QA,
                channels_AQ, channels_AA]

    for key, func, channels in zip(keys, funcs, chg_list):
        TM_tot = channels.get_T()


        # print(key)

        # # break
        # lbl = f'TMtot/{key}'
        # print(lbl)
        # f.create_dataset(lbl, data=TM_tot)


        # plt.plot(erange, imag(TM_tot[:, 0]))
        ch = list(channels.channels.items())[0][1] #ch = list(channels.items())[0][1] # take any of the channels since SFs are the same
        iImTM_tot = tm.Interpolator2D(ch.chs[0].qrange, ch.chs[0].erange, np.ascontiguousarray(np.imag(TM_tot)))


        # plt.plot(erange_dense, [iImTM_tot(0, e) for e in erange_dense])
        # plt.show()
        


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
            iImSigma = tm.Interpolator(erange, np.ascontiguousarray(res), 'cubic')
            ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in erange]
            ReSigmas += [ReSigma]

        REALs[key] = np.array(ReSigmas).transpose()            

    IMAGss += [IMAGs]
    REALss += [REALs]
    
    ImS_Q = IMAGs['QQ']  + IMAGs['QA']
    ReS_Q = REALs['QQ']  + REALs['QA']

    om0_k = np.array([np.sqrt(mQ**2 + qrange**2) for e in quark_run.erange])
    arrE = np.array([quark_run.erange for q in quark_run.qrange]).transpose()

    G_Q_new = 1/(arrE - om0_k + 0*1j*quark_run.eps - (ReS_Q + 1j*ImS_Q) + mu)


    quark_new = Particle(mQ, qrange, erange, eps=quark_run.eps, Gtab=G_Q_new, mu=mu, S=ReS_Q + 1j*ImS_Q)
    
    # quark_new.S = ReS_Q + 1j*ImS_Q
##########
    ImS_A = IMAGs['AA'] + IMAGs['AQ']
    ReS_A = REALs['AA'] + REALs['AQ']

    om0_k = np.array([np.sqrt(mQ**2 + qrange**2) for e in aquark_run.erange])
    arrE = np.array([aquark_run.erange for q in aquark_run.qrange]).transpose()

    G_A_new = 1/(arrE - om0_k + 0*1j*quark_run.eps - (ReS_A + 1j*ImS_A) - mu)
    

    aquark_new = Particle(mQ, qrange, erange, eps=aquark_run.eps, Gtab=G_A_new, mu=-mu, S=ReS_A + 1j*ImS_A)
    
    # aquark_new.S = ReS_A + 1j*ImS_A
####    
#     ImS_G = 0.05#IMAGs['GG'] + IMAGs['GQ'] + IMAGs['GA']
#     ReS_G = 0# REALs['GG'] + REALs['GQ'] + REALs['GA']
    
#     # om0_k = np.array([gluon_run.om0(gluon_run.qrange) for e in gluon_run.erange])
#     om0_k = np.array([np.sqrt(mG**2 + qrange**2) for e in gluon_run.erange])
#     arrE = np.array([gluon_run.erange for q in gluon_run.qrange]).transpose()

# #     G_G_new = 1/(arrE**2 - om0_k**2 + 2*1j*gluon_run.eps*arrE - (ReS_G + 1j*ImS_G))
#     G_G_new = 1/(arrE - om0_k + 1j*gluon_run.eps - (ReS_G + 1j*ImS_G))
# #     gluon_new = Particle(gluon_run.m, qrange, erange, eps=2e-2, stat='b', d=16, R=-2*imag(G_G_new))
#     gluon_new = Particle(mG, qrange, erange, eps=gluon_run.eps, stat='b', d=16, Gtab=G_G_new)
    
#     gluon_new.S = ReS_G + 1j*ImS_G


# #         delta = 0
#         delta += np.sqrt(np.sum((quark_new.Rtab - quark_run.Rtab)**2)) / len(erange) / len(qrange)
# #         delta += sum(quark_new.Rtab / quark_run.Rtab) / len(erange) / len(qrange)
#         delta += np.sqrt(np.sum((gluon_new.Rtab - gluon_run.Rtab)**2)) / len(erange) / len(qrange)
#         delta /= 2
#         print(delta)  

    delta = 0
    wh = quark_run.Rtab != 0
    delta += np.max(np.abs(quark_new.Rtab[wh] - quark_run.Rtab[wh]))
    wh = aquark_run.Rtab != 0
    delta += np.max(np.abs(aquark_new.Rtab[wh] - aquark_run.Rtab[wh]))
    # wh = gluon_run.Rtab != 0
    # delta += np.max(np.abs(gluon_new.Rtab[wh] - gluon_run.Rtab[wh]))
    delta /= 2


        
    if abs(delta) < thr and n_iter == 0:
        quark_run.S = quark_new.S
        aquark_run.S = aquark_new.S
        # gluon_run.S = gluon_new.S
    
    chss += [[channels_QQ]]
    pts += [[quark_run, aquark_run]]

    quark_bup = quark_run

    quark_run = quark_new
    aquark_run = aquark_new
    # gluon_run = gluon_new


    sumQ = np.trapz(quark_new.Rtab[:, 0], x=erange)
    sumA = np.trapz(aquark_new.Rtab[:, 0], x=erange)
    # sumG = np.trapz(gluon_new.Rtab[:, 0], x=erange)
    # print('delta = ', delta)
    # print('sum Q = ', sumQ)
    # print('sum G = ', sumG)

    print('T = %.3f iteration %i: delta = %f, sum Q = %f, sum A = %f'%(T, n_iter,
    delta, sumQ, sumA))

    if save_iterations:
        folder_iter = os.path.join(out_folder, Tlabel)
        if not os.path.exists(folder_iter):
            os.mkdir(folder_iter)
        fname_iter = folder_iter + '/iter_%i.hdf5'%n_iter
        # print(fname_iter)
        f_iter = h5py.File(fname_iter, 'w')
        # f.create_dataset('%i/TM'%n_iter)
        # kk = ['Q', 'Q', 'G', 'G']
        for _k, gr in zip(keys, chg_list):
            for k, chl in gr.channels.items():
                for l in range(chl.lmax + 1):
                    c = chl.chs[l]
                    f_iter.create_dataset(f'TM/{k}/{l}', data=c.TM)
                    f_iter.create_dataset(f'X/{k}/{l}', data=c.XS[0])

        for k in IMAGs.keys():
            f_iter.create_dataset(f'S/{k}', data=REALs[k] + 1j*IMAGs[k])
        f_iter.create_dataset(f'Q/G', data=quark_new.Gtab)
        f_iter.create_dataset(f'A/G', data=aquark_new.Gtab)
        # f_iter.create_dataset(f'G/G', data=gluon_new.Gtab)
        f_iter.create_dataset(f'Q/S', data=quark_new.S)
        f_iter.create_dataset(f'A/S', data=aquark_new.S)
        for _k, gr in zip(keys, chg_list):
            f_iter.create_dataset(f'G2/{_k}', data=list(gr.channels.items())[0][1].chs[0].G2)

        # f_iter.create_dataset(f'G/S', data=gluon_new.S)
        for key, func, channels in zip(keys, funcs, chg_list):
            TM_tot = channels.get_T()

            lbl = f'TMtot/{key}'
            # print(lbl)
            f_iter.create_dataset(lbl, data=TM_tot)
        f_iter.close()


    if n_iter > 1:
        if any([abs(sumQ - 1) > 1e-1]):
            # break
            raise ValueError('Sum rule not fulfilled, STOP')
        

    n_iter += 1

    if n_iter > args.max_iter:
        print('\n########## max iter reached ############\n')
        break

IMAGs = IMAGss[-1]
REALs = REALss[-1]

# kk = ['Q', 'Q', 'G', 'G']
for _k, gr in zip(keys, chg_list):
    for k, chl in gr.channels.items():
        for l in range(chl.lmax + 1):
            c = chl.chs[l]
            f.create_dataset(f'TM/{k}/{l}', data=c.TM)
            f.create_dataset(f'X/{k}/{l}', data=c.XS[0])

    f.create_dataset(f'G2/{_k}', data=list(gr.channels.items())[0][1].chs[0].G2)


for k in IMAGs.keys():
    f.create_dataset(f'S/{k}', data=REALs[k] + 1j*IMAGs[k])

f.create_dataset(f'Q/R', data=pts[-1][0].Rtab)
f.create_dataset(f'A/R', data=pts[-1][1].Rtab)
# f.create_dataset(f'G/R', data=pts[-1][2].Rtab)
# print(pts[-1][0].Gtab)

f.create_dataset(f'Q/G', data=pts[-1][0].Gtab)
f.create_dataset(f'A/G', data=pts[-1][1].Gtab)
# f.create_dataset(f'G/G', data=pts[-1][2].Gtab)

f.create_dataset(f'Q/S', data=pts[-1][0].S)
f.create_dataset(f'A/S', data=pts[-1][1].S)
# f.create_dataset(f'G/S', data=pts[-1][2].S)

f.attrs.update({'status' : 'DONE'})

f.close()

########################## Thermodynamics ##########################