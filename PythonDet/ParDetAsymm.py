from copy import deepcopy
from IPython.core.debugger import Tracer
from IPython.display import HTML
from itertools import repeat
from numpy import (abs, all, argmax, array_equal, ceil, copy, cos, flipud,
                   float64, linspace, maximum, mean, ones, pi, r_, random,
                   round, sqrt, sum, zeros)
from matplotlib.pyplot import (cla, figure, gca, gcf, Normalize, plot, show,
                               subplots, subplot, xlabel, ylabel)
from matplotlib import animation, rc
from multiprocessing import Pool
from pickle import dump, HIGHEST_PROTOCOL
from scipy import arcsin
from scipy.special import ellipe
import time

class ParSim(object):

    def __init__(self, bc='NEU', dt=0.01, grid_size=100, param_dict=None,
                 ss_prec=1.001, T=1000, store_interval=10):

        if param_dict is None:
            param_dict = {'alpha': 1, 'beta': 2, 'dA': 0.28, 'dP': 0.15,
                          'kAP': 0.19, 'kPA': 2, 'koffA': 0.0054,
                          'koffP': 0.0073, 'konA': 0.00858, 'konP': 0.0474,
                          'Ptot': 1, 'ratio': 1.56, 'sys_size': 134.6/2}
        self.dt = dt  # time step
        self.grid_size = grid_size  # length of grid
        self.ss_prec = ss_prec
        self.store_interval = store_interval
        self.T = T  # wall time
        self.dx = param_dict['sys_size']/self.grid_size  # space step
        self.n = int(self.T/self.dt)
        self.finished_in_time = 0

        # PAR system
        self.alpha = param_dict['alpha']
        self.beta = param_dict['beta']
        self.dA = param_dict['dA']
        self.dP = param_dict['dP']
        self.kAP = param_dict['kAP']
        self.kPA = param_dict['kPA']
        self.konA = param_dict['konA']
        self.konP = param_dict['konP']
        self.koffA = param_dict['koffA']
        self.koffP = param_dict['koffP']
        self.Ptot = param_dict['Ptot']
        self.ratio = param_dict['ratio']
        self.sys_size = param_dict['sys_size']
        if bc is 'NEU':
            # Mutliply by two, because s_to_v takes the entire circumference,
            # not half of it.
            self.StoV = s_to_v('Circumference', [self.sys_size*2, 15/27])
        # self.StoV = 0.174
        self.Atot = self.ratio * self.Ptot

        # Wave pinning
        # self.Da = 0.0
        # self.Db = 10
        # self.delta = 1
        # self.gamma = 1
        # self.k0 = 0.067
        # self.a = 0.2683312
        # self.b = 2.0
        # self.K = 1

    def pickle_data(self, fname):
        '''
        Pickle data to store on disk.

        path     save pickled object elsewhere other than in the directory
                 containing the imaging data.
        postfix  add a postfix to the filename
        '''
        with open(fname + '.pickle', 'wb') as f:
            # Pickle self using the highest protocol available.
            dump(self, f, HIGHEST_PROTOCOL)

    def plot_steady_state(self, c1='salmon', c2='cornflowerblue',
                          lab='_nolabel_', alpha=1.0, lw=1.5, linestyle='-', norm=False,
                          printLast=False):
        # x = linspace(0, self.grid_size*self.dx, self.grid_size)
        x = linspace(-self.grid_size*self.dx/2, self.grid_size*self.dx/2,
                     self.grid_size)
        # for ii, i in enumerate(round(linspace(1, self.A.shape[1], 10))):
        if printLast:
            t_points = [1]
        else:
            t_points = round(linspace(1, self.A.shape[1], 10))

        for ii, i in enumerate(t_points):
            ya = self.A[:, -int(i)]
            yp = self.P[:, -int(i)]
            if norm:
                ya = (ya-min(ya))/(max(ya)-min(ya))
                yp = (yp-min(yp))/(max(yp)-min(yp))
            plot(x, ya, c1, label='_nolabel_', alpha=alpha/(ii+1)/3+0.34, lw=lw, linestyle=linestyle)
            plot(x, yp, c2, label=lab, alpha=alpha/(ii+1)/3+0.34, lw=lw, linestyle=linestyle)
        if self.t[-1] >= self.T:
            print('Steady state not reached, plotting last time point.')
        # show()

    def set_init_profile(self):
        self.A = ones((self.grid_size, int(1.2*self.T/self.store_interval)))
        self.P = ones((self.grid_size, int(1.2*self.T/self.store_interval)))
        self.A[int(round(self.grid_size/2)):, 0] = 0
        self.P[0:int(round(self.grid_size/2)), 0] = 0

        # self.P = 4*ones((self.grid_size, int(1.2*self.T/self.store_interval)))
        # self.P[0:int(round(0.8*self.grid_size)), 0] = 0
        # self.A = 4*ones((self.grid_size, int(1.2*self.T/self.store_interval)))
        # self.A[int(round(0.2*self.grid_size)):, 0] = 0

    def set_init_profile_wave_pin(self):
        o = ones((self.grid_size, int(10*self.T/self.store_interval)))
        self.A = self.a*o
        self.P = self.b*o
        self.stim = zeros(self.grid_size)
        ind = int(self.grid_size/10)
        self.stim[0:ind] = 1 + cos(linspace(0, pi, ind))

    def save_movie(self, fname=None, dpi=200, everynth=10):
        size_factor = self.sys_size/(self.grid_size*self.dx)
        fig = figure(figsize=(5, 4.0), facecolor='black')
        ax = gca()
        ax.set_xlim((0, self.grid_size*self.dx))
        ax.set_ylim((0, 1.1*max(self.A.max(), self.P.max())))
        line, = ax.plot([], [], lw=2)
        lineA, = ax.plot([], [], lw=2, color='red')
        lineP, = ax.plot([], [], lw=2, color='dodgerblue')
        yellow = (242/255, 204/255, 32/255)
        t_text = ax.text(0.1*self.grid_size*self.dx,
                         0.2, r't = ' + str(0) + ' s',
                         fontsize=12, color=yellow)
        ax.tick_params(axis='x', colors=yellow,
                       which='both', labelsize=10)
        ax.tick_params(axis='y', colors=yellow,
                       which='both', labelsize=10)
        xlabel(r'Position $[\mu m]$', color=yellow,
               fontsize=12)
        ylabel(r'Concentration $[a.u.]$', color=yellow,
               fontsize=12)
        ax.spines['left'].set_color(yellow)
        ax.spines['bottom'].set_color(yellow)
        ax.set_facecolor((0, 0, 0))
        gcf().subplots_adjust(bottom=0.19, left=0.16)

        # initialization function: plot the background of each frame
        def init():
            # fig.set_facecolor((0, 0, 0))
            line.set_data([], [])
            return (line,)

        # animation function. This is called sequentially
        def animate(i):
            x = linspace(0, self.grid_size*self.dx, self.A.shape[0])
            yA = self.A[:, i]
            yP = self.P[:, i]
            lineA.set_data(x, yA)
            lineP.set_data(x, yP)
            t_text.set_text(r't = ' + str(round(self.t_stored[i], 0)) + ' s')
            return (line,)

        # call the animator
        i_e = self.A.shape[1]
        a = animation.FuncAnimation(fig, animate, init_func=init,
                                    frames=range(0, i_e, everynth),
                                    interval=100, blit=False)
        if fname is None:
            a.save('lines.mp4', savefig_kwargs={'facecolor': 'black'}, dpi=dpi)
        else:
            a.save(fname + '.mp4', savefig_kwargs={'facecolor': 'black'},
                   dpi=dpi)

    def simulate(self):

        def neu(An, Pn):
            '''
            Von Neumann conditions, setting the derivative at the
            boundaries to zero.
            '''
            delA = laplacianNEU(An, self.dx)
            delP = laplacianNEU(Pn, self.dx)
            Acy = self.Atot - self.StoV*sum(flipud(An), 0)/self.grid_size
            Pcy = self.Ptot - self.StoV*sum(Pn, 0)/self.grid_size
            # Defining Ra and Rp separately is necessary in order to not
            # update An before Pn (which would then have an effect on Pn
            # in the same iteration, making the system asymmetric)
            Ra = (self.dA*delA-self.koffA*An[1:-1]+self.konA*Acy -
                  self.kAP*Pn[1:-1]**self.alpha*An[1:-1])
            Rp = (self.dP*delP-self.koffP*Pn[1:-1]+self.konP*Pcy -
                  self.kPA*An[1:-1]**self.beta*Pn[1:-1])
            # Return arrays with first and last two elements equal
            # respecitvely, to impose zero derivatives
            return r_[Ra[0], Ra, Ra[-1]], r_[Rp[0], Rp, Rp[-1]]

        # Adaptive step size parameters
        atol = 0.000001
        rtol = 0.000001

        # 5TH ORDER RK COEFFICIENTS for Dormand-Prince
        a21, a31, a32, a41, a42, a43 = 1/5, 3/40, 9/40, 44/45, -56/15, 32/9
        a51, a52, a53, a54 = 19372/6561, -25360/2187, 64448/6561, -212/729
        a61, a62, a63 = 9017/3168, -355/33, 46732/5247
        a64, a65 = 49/176, -5103/18656
        a71, a72, a73, a74 = 35/384, 0, 500/1113, 125/192
        a75, a76 = -2187/6784, 11/84

        b1, b2, b3, b4, b5 = 35/384, 0, 500/1113, 125/192, -2187/6784
        b6, b7 = 11/84, 0

        bs1, bs2, bs3, bs4 = 5179/57600, 0, 7571/16695, 393/640
        bs5, bs6, bs7 = -92097/339200, 187/2100, 1/40

        c2, c3, c4, c5, c6, c7 = 1/5, 3/10, 4/5, 8/9, 1, 1

        # Set initial profile
        ################################################################
        self.set_init_profile()
        # self.set_init_profile_wave_pin()
        A0 = self.A[:, 0]
        P0 = self.P[:, 0]
        self.t = [0]
        self.t_stored = [0]
        tnext = 0
        i = 1
        # self.totalerror = [0]
        # self.errRatio = []
        self.reject = 0
        self.no_reject = 0
        while self.t[-1] < self.T:
            # Calculate increments for RK45
            if (self.t[-1] == 0) or not (Pn_new == P0[1]).any():
                A1, P1 = neu(A0, P0)
            else:
                A1, P1 = A7, P7
            A2, P2 = neu(A0+self.dt*(a21*A1), P0+self.dt*(a21*P1))
            A3, P3 = neu(A0+self.dt*(a31*A1+a32*A2),
                         P0+self.dt*(a31*P1+a32*P2))
            A4, P4 = neu(A0+self.dt*(a41*A1+a42*A2+a43*A3),
                         P0+self.dt*(a41*P1+a42*P2+a43*P3))
            A5, P5 = neu(A0+self.dt*(a51*A1+a52*A2+a53*A3+a54*A4),
                         P0+self.dt*(a51*P1+a52*P2+a53*P3+a54*P4))
            A6, P6 = neu(A0+self.dt*(a61*A1+a62*A2+a63*A3+a64*A4+a65*A5),
                         P0+self.dt*(a61*P1+a62*P2+a63*P3+a64*P4+a65*P5))
            A7, P7 = neu(A0+self.dt*(a71*A1+a73*A3+a74*A4+a75*A5+a76*A6),
                         P0+self.dt*(a71*P1+a73*P3+a74*P4+a75*P5+a76*P6))
            # Update concentrations using A1-A6 and P1-P6, coefficient for
            # A7 and P7 is 0.
            An_new = A0+self.dt*(b1*A1+b3*A3+b4*A4+b5*A5+b6*A6)  # b2/7=0
            Pn_new = P0+self.dt*(b1*P1+b3*P3+b4*P4+b5*P5+b6*P6)  # b2/7=0
            # Compute difference between fourth and fifth order
            deltaAnerr = max(abs((b1-bs1)*A1+(b3-bs3)*A3+(b4-bs4)*A4 +
                             (b5-bs5)*A5+(b6-bs6)*A6-bs7*A7))  # b7 is zero
            deltaPnerr = max(abs((b1-bs1)*P1+(b3-bs3)*P3+(b4-bs4)*P4 +
                             (b5-bs5)*P5+(b6-bs6)*P6-bs7*P7))  # b7 is zero
            # Get maximum concentrations for An and Pn
            yAn = maximum(max(abs(An_new)), max(abs(A0)))
            yPn = maximum(max(abs(Pn_new)), max(abs(P0)))
            # Get error scale, combining relative and absolute error
            scaleAn = atol+yAn*rtol
            scalePn = atol+yPn*rtol
            # Compute total error as norm of maximum errors for each
            # species scaled by the error scale
            totalerror = sqrt(1/2*((deltaAnerr/scaleAn)**2 +
                              (deltaPnerr/scalePn)**2))
            # Compute new timestep
            dtnew = 0.8*self.dt*abs(1/totalerror)**(1/5)
            # Upper and lower bound for timestep to avoid changing too fast
            if dtnew > 10*self.dt:
                dtnew = 10*self.dt
            elif dtnew < self.dt/5:
                dtnew = self.dt/5
            # Set timestep for next round
            self.dt = dtnew
            # Accept step if error is on the order of error scale or below
            if totalerror < 1:
                self.t.append(self.t[-1]+self.dt)
                # self.errRatio.append(max(A0/An_new))
                # self.totalerror.append(totalerror)
                if self.t[-1] > tnext:
                    self.A[:, i] = An_new
                    self.P[:, i] = Pn_new
                    self.t_stored.append(self.t[-1])
                    i += 1
                    tnext = tnext + self.store_interval
                # Break if things change by less than 0.1% over
                # the course of 1 min or maximum difference between
                # concentration of a species is less than 5%.
                # if ((max(An_new)/min(An_new) < 1.05) and
                #    (max(Pn_new)/min(Pn_new) < 1.05)):
                #     print('Unpolarized, not necessarily close to SS!')
                #     self.finished_in_time = 1
                #     break
                # elif (max(max(abs(A0/An_new))**(60/dtnew),
                #           max(abs(P0/Pn_new))**(60/dtnew))) < self.ss_prec:
                #     print('Steady state reached.')
                #     self.finished_in_time = 2
                #     break
                A0 = An_new
                P0 = Pn_new
                self.no_reject += 1
            else:
                self.reject += 1

            if (An_new[0] < Pn_new[0]) or (Pn_new[-1] < An_new[-1]):
                break
        self.A[:, i] = An_new
        self.P[:, i] = Pn_new
        self.t_stored.append(self.t[-1])
        # Cut A and P to not include any zeros from preallocation
        self.A = self.A[:, ~all(self.A == 1, axis=0)]
        self.P = self.P[:, ~all(self.P == 1, axis=0)]

        if (self.A[0, -1] < self.P[0, -1]):
            self.finished_in_time = 0
        elif (self.P[-1, -1] < self.A[-1, -1]):
            self.finished_in_time = 1
        else:
            self.finished_in_time = 2


class Sim_Container:

    def __init__(self, dicts='None', no_workers=8):
        self.dicts = dicts
        self.no_workers = no_workers

    def init_simus(self):
        self.simList = []
        for dic in self.dicts:
            self.simList.append(ParSim(param_dict=dic))

    def pickle_data(self, fname):
        '''
        Pickle data to store on disk.

        path     save pickled object elsewhere other than in the directory
                 containing the imaging data.
        postfix  add a postfix to the filename
        '''
        with open(fname + '.pickle', 'wb') as f:
            # Pickle self using the highest protocol available.
            dump(self, f, HIGHEST_PROTOCOL)

    def plot_all_ss(self):
        for i in self.simus:
            i.plot_steady_state()
            show()

    def run_simus(self):

        if self.no_workers == 1:
            self.simus = []
            for i in self.simList:
                self.simus.append(sim_indiv(i))
        else:
            # Create pool, use starmap to pass more than one parameter, do work
            pool = Pool(processes=self.no_workers)
            res = pool.map(sim_indiv, self.simList)
            self.simus = res
            pool.close()
            pool.join()

class Dosage_CPSS:
    '''
    Sweep dosage cpss space for solution where system is in either of the 
    following states:
        # 0 - unpolarized P won
        # 1 - unpolarized A won
        # 2 - polarized

    '''
    def __init__(self, ratio_above, ratio_below, precision, sz_range,
                 base_set=None, max_iter=10):
        self.ratio_above = ratio_above
        self.ratio_below = ratio_below
        self.max_iter = max_iter
        self.precision = precision
        self.sz_range = sz_range

        if base_set is None:
            self.base_set = {'alpha': 2, 'beta': 2, 'dA': 0.1, 'dP': 0.1,
                     'kAP': 1, 'kPA': 1, 'koffA': 0.005, 'koffP': 0.005,
                     'konA': 0.006, 'konP': 0.006, 'Ptot': 1, 'ratio': 1.0,
                     'sys_size': 30}
        else:
            self.base_set = base_set

    def pickle_data(self, fname):
        '''
        Pickle data to store on disk.

        path     save pickled object elsewhere other than in the directory
                 containing the imaging data.
        postfix  add a postfix to the filename
        '''
        with open(fname + '.pickle', 'wb') as f:
            # Pickle self using the highest protocol available.
            dump(self, f, HIGHEST_PROTOCOL)

    def simulate(self, dom='A'):

        if dom == 'A':
            large = 0
            small = 1
        elif dom == 'P':
            large = 1
            small = 0

        sz = linspace(self.sz_range[0], self.sz_range[1], self.precision)

        sym = [deepcopy(self.base_set) for _ in sz]

        for i, dic in enumerate(sym):
            sym[i]['sys_size'] = sz[i]

        ratio_below = [self.ratio_below for _ in sym]
        ratio_above = [self.ratio_above for _ in sym]
        ratio = [self.ratio_above for _ in sym]
        outcome = [3 for _ in sym]

        start_time = time.time()

        for j in range(self.max_iter):
            for i in range(len(sym)):
                sym[i]['ratio'] = ratio[i]
            s = Sim_Container(sym, no_workers=8)
            s.init_simus()
            s.run_simus()
            ratio_temp = copy(ratio)
            for i in range(len(sym)):
                outcome[i] = s.simus[i].finished_in_time
                if (dom=='A') and ((outcome[i] == 2) or (outcome[i] == 0)):
                    ratio[i] = ratio[i] + (ratio_above[i]-ratio[i])/2
                    ratio_below[i] = ratio_temp[i]
                elif (dom == 'A') and (outcome[i] == 1):
                    ratio[i] = ratio[i] - (ratio[i]-ratio_below[i])/2
                    ratio_above[i] = ratio_temp[i]
                elif (dom=='P') and ((outcome[i] == 2) or (outcome[i] == 1)):
                    # Tracer()()
                    ratio[i] = ratio[i] - (ratio[i]-ratio_below[i])/2
                    ratio_above[i] = ratio_temp[i]
                elif (dom == 'P') and (outcome[i] == 0):
                    ratio[i] = ratio[i] + (ratio_above[i]-ratio[i])/2
                    ratio_below[i] = ratio_temp[i]
            # print(ratio)

        self.time = time.time() - start_time

        if dom == 'A':
            self.ratioA = ratio
        elif dom == 'P':
            self.ratioP = ratio

        
            # if abs(ratio_above[i]-ratio_below[i])/ratio_below[i] < 0.002:
            #     break

        self.simus = s

def laplacianNEU(Z, dx):
    Zleft = Z[0:-2]
    Zright = Z[2:]
    Zcenter = Z[1:-1]
    return (Zleft + Zright - 2*Zcenter) / dx**2


def s_to_v(mode, describers):
    ''''Calculate surface area to volume ratio for ellipsoid'''

    if mode is 'Circumference':
        circ = describers[0]
        aspRatio = describers[1]
        e = (1 - aspRatio**2)
        EllInt2nd = ellipe(e)
        lA = circ/(4*EllInt2nd)
        sA = lA * aspRatio
    elif mode is 'Axes':
        sA = describers[0]
        lA = describers[1]
    else:
        print('s_to_v(): Mode ' + str(mode) + 'not known, returning.')
        return

    e = sqrt(1 - sA**2/lA**2)
    S = 2*3.14159265*sA**2*(1 + lA/(sA*e)*arcsin(e))
    V = 4/3*3.14159265*sA**2*lA
    return S/V


def sim_indiv(Simu):
    Simu.simulate()
    return Simu
