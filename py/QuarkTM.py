import numpy as np
from scipy import signal
import TMQGP as tm
import tqdm.notebook as tqdm
from syntax_sugar import pipe, END

from syntax_sugar import thread_syntax as t, process_syntax as p

class IterController:
    def __init__(self) -> None:
        pass

#TODO class for spectral functions and T-matrix

# class Runner:
#     def __init__(self, T, channel='ff', Fa=1, G=8, L=0.3, m=0.5, eps=5e-2, erange=None, qrange=None, inters=None, screen=0):
#         self.m = m
#         self.T = T
#         self.L = L
#         self.G = G
#         self.eps = eps
#         self.screen = screen
#         self.Fa = Fa
#         self.channel = channel

#         if len(set(channel)) == 2:
#             self.statistics = 'f'
#         else:
#             self.statistics = 'b'

#         if channel == 'bb':
#             self.func = tm.sigma_bb
#         elif channel == 'bf':
#             self.func = tm.sigma_bf
#         elif channel == 'ff':
#             self.func = tm.sigma_ff
#         elif channel == 'fb':
#             self.func = tm.sigma_fb
        
#         if erange is None:
#             self.erange = np.linspace(-5, 5, 100)
#         else: self.erange = erange
#         if qrange is None:
#             self.qrange = np.linspace(0, 5, 50)
#         else: self.qrange = qrange
#         if inters is None:
#             self.iV = tm.Interpolator(self.qrange, self.v(self.qrange), 'linear')
#             self.iOm = tm.Interpolator(self.qrange, self.om0(self.qrange), 'linear')
#             arrGqq = np.array([[self.Gqq0(e, k) for k in self.qrange] for e in self.erange])
#             self.iReGqq = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.real(arrGqq)))
#             self.iImGqq = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag(arrGqq)))
#             self.inters = [self.iV, self.iOm, self.iReGqq, self.iImGqq]
#             arrGq = np.array([[-2*np.imag(self.Gq0(_e, _k)) for _k in self.qrange] for _e in self.erange])
#             self.iGq = tm.Interpolator2D(self.qrange, self.erange, (arrGq))
#         else:
#             self.inters = inters
#             self.iV, self.iOm, self.iReGqq, self.iImGqq, self.iGq = inters

#     def v(self, q):
#         # return self.G * self.L**2 / (q**2 + self.L**2)#exp(-k**2 / L**2)
#         mult = 1 / (1 + self.screen*self.T**2)
#         return self.Fa * self.G * np.exp(-q**2 / (self.L * mult)**2)

#     def V(self, q, q1):
#         return self.v(q) * self.v(q1)

#     def om0(self, k, m=None):
#         if m == None:
#             m = self.m
#         return np.sqrt(m**2 + k**2)

#     def Gqq0(self, E, k):
#         return 1/2/(E/2 - self.om0(k) + 1j*self.eps)

#     def Gq0(self, E, k):    
#         return 1 / (E - self.om0(k) + 1j*self.eps)

#     def populate_T(self):
#         self.TM = np.array([[tm.T_solve(E, k, k, self.T, self.iV, self.iOm, self.iReGqq, self.iImGqq) for k in (self.qrange)]
#                 for E in tqdm.tqdm(self.erange)])

#         self.iImT = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag((self.TM))))

#     def populate_S(self):
#         ress_cm = np.array([[self.func(e, q, self.T, self.iImT, self.iGq) for q in self.qrange]
#         for e in tqdm.tqdm(self.erange)])
#         self.ImS = ress_cm

#         ReSigmas = []

#         for res in tqdm.tqdm(ress_cm.transpose()):
#             iImSigma = tm.Interpolator(self.erange, np.ascontiguousarray(res), 'cubic')
#             ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in self.erange]
#             ReSigmas += [ReSigma]

#         self.ReS = np.array(ReSigmas).transpose()

#     def nf(self, e, T):
#         return 1/(np.exp(e/T) + 1)

#     def nb(self, e, T):
#         # if e == 0:
#         #     return 10
#         res = 1/(np.exp(e/T) - 1)
#         res[e == 0] = 10
#         return res

#     def next_step(self):
#         om0_k = np.array([self.om0(self.qrange) for e in self.erange])
#         arrE = np.array([self.erange for q in self.qrange]).transpose()
#         Gq_res = 1/(arrE - om0_k + 1j*self.eps - (self.ReS + 1j*self.ImS))
#         self.rhoQ_new = -2*np.imag(Gq_res)
#         de = self.erange[1] - self.erange[0]
#         self.ImGqq_new = -np.array([
#             signal.convolve(_rho, _rho, mode='same')* de / 4 / np.pi for _rho in self.rhoQ_new.transpose()
#         ])

#         if self.channel[0] == 'f':
#             cf1 = lambda z, T: self.nf(z, T)
#         else:
#             cf1 = lambda z, T: -self.nb(z, T)

#         if self.channel[1] == 'f':
#             cf2 = lambda z, T: self.nf(z, T)
#         else:
#             cf2 = lambda z, T: -self.nb(z, T)

#         self.ImGqq_new += np.array([
#             signal.convolve(_rho, _rho*cf1(self.erange, self.T), mode='same')* de / 4 / np.pi for _rho in self.rhoQ_new.transpose()
#         ])

#         self.ImGqq_new += np.array([
#             signal.convolve(_rho, _rho*cf2(self.erange, self.T), mode='same')* de / 4 / np.pi for _rho in self.rhoQ_new.transpose()
#         ])

#         ReGqq_new = []

#         for _im in self.ImGqq_new:
#             iImNew = tm.Interpolator(self.erange, _im, 'linear')
#             ReNew = np.array([tm.ReSigmaKK(e, iImNew) for e in self.erange])
#             ReGqq_new += [ReNew]

#         self.ReGqq_new = np.array(ReGqq_new)
        

#         iImGqq_new = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(self.ImGqq_new.transpose()))
#         iReGqq_new = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(self.ReGqq_new.transpose()))

#         iGq_new = tm.Interpolator2D(self.qrange, self.erange, self.rhoQ_new)

#         r_new = Runner(self.T, self.channel, self.Fa, self.G, self.L, self.m, self.eps, self.erange, self.qrange, 
#                     inters=[self.iV, self.iOm, iReGqq_new, iImGqq_new, iGq_new])

#         return r_new


        

class Particle:
    def __init__(self, m, qrange, erange, R=None, stat='f', eps=5e-2, d=6, 
            propagator=3, Gtab=None):
        self.m = m
        self.eps = eps
        self.R = None
        self.Rtab = None
        self.qrange = qrange
        self.erange = erange
        self.stat = stat
        self.propagator = propagator
        self.d = d
        if Gtab is not None:
            self.Gtab = Gtab
            self.set_R(-2*np.imag(Gtab))
            self.iImG = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag(Gtab)))
            self.iReG = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.real(Gtab)))
        else:
            if R is None:
                self.Gtab =np.array([[self.G0(_e, _k)
                                for _k in self.qrange]
                                for _e in self.erange])
                self.iImG = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag(self.Gtab)))
                self.iReG = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.real(self.Gtab)))
                arrR0 = - 2 * np.imag(self.Gtab)
                self.set_R(arrR0)
            else:
                self.set_R(R)
    
    def set_R(self, Rtab):
        if Rtab.shape[0] != len(self.erange) or Rtab.shape[1] != len(self.qrange):
            raise
        self.Rtab = Rtab
        self.R = tm.Interpolator2D(self.qrange, self.erange, 
        np.ascontiguousarray(Rtab))

        
    def om0(self, q):
        return np.sqrt(self.m**2 + q**2)
        
    def G0(self, E, q):
        # if self.propagator == 'Th':
        if self.stat == 'f':
            # return 1 / (E - self.om0(q) + 1j*self.eps*(1 + np.tanh(E/0.001))/2)
            return 1 / (E - self.om0(q) + 1j*self.eps*(E))
        # elif self.propagator == 'BBS':
        elif self.stat == 'b':
            if self.propagator == 1:
                # return 1 / (E - self.om0(q) + 1j*self.eps*E)
                return 1 / (E - self.om0(q) + 1j*self.eps*(np.tanh(E/0.5)**3))
            elif self.propagator == 2:
                return 1/((E)**2 - self.om0(q)**2 + 2j*self.eps*np.sign(E))
            elif self.propagator == 3:
                res = 1 / (E - self.om0(q) + 1j*self.eps*(np.tanh(E/0.5)**3)) * (E>0)
                res += -1 / (-E - self.om0(q) + 1j*self.eps*(np.tanh(-E/0.5)**3)) * (E<0)
                return res  
                    
        else:
            raise

    def __copy__(self):
        return Particle(self.m, self.qrange, self.erange, self.Rtab, 
                       self.stat, self.eps)

class Channel:
    def __init__(self, p_i: Particle, p_j: Particle, T : float,  
                 Fa=1, G=6, L=0.5, screen=4, ds=1, da=1, calc=2, 
                 do_rel=1, parallel=-1,
                 test_potential=0):
        self.p_i = p_i
        self.p_j = p_j
        
        self.test_potential = test_potential

        if p_i.stat == 'b' and p_j.stat == 'b':
            if calc == 1:
                self.func = tm.sigma_bb
            elif calc == 2:
                self.func = tm.sigma_bb2
            elif calc == 3:
                self.func = tm.sigma_bb3
            else:
                raise

            self.Tfunc = tm.T_solve_BB
            
        elif p_i.stat == 'b' and p_j.stat == 'f':
            self.func = tm.sigma_bf
            self.Tfunc = tm.T_solve_BF
            
        elif p_i.stat == 'f' and p_j.stat == 'f':
            self.func = tm.sigma_ff_onshell
            self.Tfunc = tm.T_solve

        elif p_i.stat == 'f' and p_j.stat == 'b':
            self.func = tm.sigma_fb
            self.Tfunc = tm.T_solve_BF
        else:
            raise
            

        if self.test_potential:
            self.Tfunc = tm.T_solve_explicit
        self.ds = ds
        self.da = da
        self.erange = p_i.erange ### TODO: dirty
        self.qrange = p_i.qrange
        self.screen = screen
        self.T = T
        self.Fa = Fa
        self.G = G
        self.L = L
        self.parallel = parallel
        self.iV = tm.Interpolator(self.qrange, self.v(self.qrange, do_rel=do_rel), 'linear')
        self.iOm = tm.Interpolator(self.qrange, self.p_i.om0(self.qrange), 'linear')
        self.set_G2()
        
    def v(self, q, do_rel=1):
        mult = 1 / (1 + self.screen*self.T**2)

        ### Relativistic correction to the vertex: R_S(q, q') from Liu&Rapp
        rel_factor = 1
        if do_rel:
            rel_factor = np.sqrt(self.p_i.m * self.p_j.m / self.p_i.om0(q) / self.p_j.om0(q))

        return rel_factor * np.sqrt(self.Fa) * self.G * np.exp(-q**2 / (self.L * mult)**2)

    
    
    def nf(self, e, T):
        return 1/(np.exp(e/T) + 1)

    def nb(self, e, T):
        res = np.real(1/(np.exp((e + 1e-15j)/T) - 1))
        # res[e == 0] = 10
        return res

    def set_G2(self):
        de = self.p_i.erange[1] - self.p_i.erange[0]



        self.ImG2 = -np.array([
            signal.convolve(r1, r2, mode='same') * de / 4 / np.pi
            for r1, r2 in zip(self.p_i.Rtab.transpose(), self.p_j.Rtab.transpose())
        ]).transpose()
        
        

        if self.p_i.stat == 'f':
            cf1 = lambda z, T: self.nf(z, T)
        else:
            cf1 = lambda z, T: -self.nb(z, T)

        if self.p_j.stat == 'f':
            cf2 = lambda z, T: self.nf(z, T)
        else:
            cf2 = lambda z, T: -self.nb(z, T)
            
            
        
        self.ImG2 += np.array([
            signal.convolve(r1*cf1(self.erange, self.T), r2, 
                            mode='same')* de / 4 / np.pi 
            for r1, r2 in zip(self.p_i.Rtab.transpose(), self.p_j.Rtab.transpose())
        ]).transpose()
        
        self.ImG2 += np.array([
            signal.convolve(r1, r2*cf2(self.erange, self.T), 
                            mode='same')* de / 4 / np.pi
            for r1, r2 in zip(self.p_i.Rtab.transpose(), self.p_j.Rtab.transpose())
        ]).transpose()
        
            
            
        ReG2 = []
            
        for _im in self.ImG2.transpose():
            iImNew = tm.Interpolator(self.erange, _im, 'linear')
            ReNew = np.array([tm.ReSigmaKK(e, iImNew) for e in self.erange])
            ReG2 += [ReNew]

        self.ReG2 = np.array(ReG2).transpose()
        
        self.iReG2 = tm.Interpolator2D(self.qrange, self.erange, 
                            np.ascontiguousarray(self.ReG2))
        self.iImG2 = tm.Interpolator2D(self.qrange, self.erange, 
                            np.ascontiguousarray(self.ImG2))
        
    def G20(self, E, k):
        return 1/2/(E/2 - self.p_i.om0(k) + 1j*self.p_i.eps)
    
    def populate_T_old(self):
        # if self.p_i.stat == 'b' and self.p_j.stat == 'b':
        self.TM = np.array([[self.Tfunc(E, k, k, self.T, self.iV, self.iOm, self.iReG2, self.iImG2, 5,
            int(np.sign(self.G))) for k in self.qrange]
            for E in tqdm.tqdm(self.erange)])
        # else:
        #     self.TM = np.array([[tm.T_solve(E, k, k, self.T, self.iV, self.iOm, self.iReG2, self.iImG2) for k in (self.qrange)]
        #         for E in tqdm.tqdm(self.erange)])

        self.iImT = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag((self.TM))))
        
    def populate_T(self):
        # if self.p_i.stat == 'b' and self.p_j.stat == 'b':

        if not self.test_potential:
            self.X = np.array([[tm.x_solve(E, k, k, self.T, self.iV, self.iOm, self.iReG2, self.iImG2, 5,
                int(np.sign(self.G))) for k in self.qrange]
                for E in tqdm.tqdm(self.erange)])
        # else:
        #     self.TM = np.array([[tm.T_solve(E, k, k, self.T, self.iV, self.iOm, self.iReG2, self.iImG2) for k in (self.qrange)]
        #         for E in tqdm.tqdm(self.erange)])
            v1v2 = np.sign(self.G)*self.v(self.qrange)**2
            self.TM = - 4*np.pi*v1v2 / (1 - self.X)

        else:
            print('XXXX')
            self.X = np.array([[tm.J_solve_explicit(E, k, k, self.T, self.iV, self.iOm, self.iReG2, self.iImG2, 5,
                int(np.sign(self.G))) for k in self.qrange]
                for E in tqdm.tqdm(self.erange)])

            v1v2 = self.v(self.qrange)**2 * np.sign(self.G)
            self.TM = - v1v2 / (np.array([self.erange**2 for q in self.qrange]).transpose() - 0.7**2 - self.X)


        self.iImT = tm.Interpolator2D(self.qrange, self.erange, np.ascontiguousarray(np.imag((self.TM))))
        
    def get_S_q(self, q):
        res = pipe(self.erange) | p[lambda z: self.func(z, q, self.T, self.iImT,
                                             self.p_j.iImG)]*(self.parallel//1) | END
        return np.array(res)
        
    def populate_S(self):
        # print(self.func(0, 0, self.T, self.iImT, self.p_j.R))
        if self.parallel < 2:
            if self.func != tm.sigma_ff_onshell:
                ress_cm = np.array([[self.func(e, q, self.T, self.iImT, self.p_j.R) for q in self.qrange]
                    for e in tqdm.tqdm(self.erange)])
            else:
                eps1 = tm.Interpolator(self.qrange, self.p_i.om0(self.qrange))
                eps2 = tm.Interpolator(self.qrange, self.p_j.om0(self.qrange))
                ress_cm = np.array([[self.func(e, q, self.T, self.iImT, self.p_j.R,
                eps1, eps2) for q in self.qrange]
                    for e in tqdm.tqdm(self.erange)])
            self.ImS = ress_cm
        else:
            ress_cm = np.array([self.get_S_q(q) for q in tqdm.tqdm(self.qrange)])
            
            self.ImS = ress_cm.transpose()

        # Antisymmetrize the BF case 

        if self.p_i.stat == 'b' and self.p_j.stat == 'f':
            ress_new = []
            for res in ress_cm.transpose():
                iImSigma = tm.Interpolator(self.erange, np.ascontiguousarray(res), 'cubic')
                ImS_new = []
                for e in self.erange:
                    if e < 0:
                        ImS_new += [-iImSigma(-e)]
                    else:
                        ImS_new += [iImSigma(e)]
                ress_new += [ImS_new]
            
            ress_cm = np.array(ress_new).transpose()
            self.ImS = ress_cm

        
        # Calculate the real part

        ReSigmas = []

        for res in tqdm.tqdm(ress_cm.transpose()):
            iImSigma = tm.Interpolator(self.erange, np.ascontiguousarray(res), 'cubic')
            ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in self.erange]
            ReSigmas += [ReSigma]

        self.ReS = np.array(ReSigmas).transpose()

def get_S_q(ch, q):
    res = pipe(ch.erange) | p[lambda z: ch.func(z, q, ch.T, ch.iImT, ch.p_j.iImG)]*(ch.parallel//1) | END
    return res


# TODO: Very goddamn dirty
def set_S_q(ch):
    ress_cm = np.array([get_S_q(ch, q) for q in tqdm.tqdm(ch.qrange)])

    ch.ImS = ress_cm

    ReSigmas = []

    for res in tqdm.tqdm(ress_cm.transpose()):
        iImSigma = tm.Interpolator(ch.erange, np.ascontiguousarray(res), 'cubic')
        ReSigma = [tm.ReSigmaKK(e, iImSigma) for e in ch.erange]
        ReSigmas += [ReSigma]

    ch.ReS = np.array(ReSigmas).transpose()
