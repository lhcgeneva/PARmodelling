from IPython.core.debugger import Tracer
from IPython.display import HTML
from itertools import repeat
from numpy import (argmax, array_equal, ceil, copy, flipud, float64, linspace,
                   mean, ones, r_, random, round, sqrt, sum)
from matplotlib.pyplot import (cla, figure, gca, gcf, Normalize, plot, show,
                               subplots, subplot, xlabel, ylabel)
from matplotlib import animation, rc
from multiprocessing import Pool
from pickle import dump, HIGHEST_PROTOCOL
from scipy import arcsin
from scipy.special import ellipe


class ParSim(object):

    def __init__(self, applyNoise=False, alpha=1, bc='NEU', beta=2,
                 dA=0.28, dP=0.15, dt=0.05, grid_size=100, kAP=0.19, kPA=2,
                 koffA=0.0054, koffP=0.0073, konA=0.00858, konP=0.0474, Ptot=1,
                 ratio=1.56, save_nth=200, ss_prec=0.0001, StoV=0.174,
                 sys_size=134.6/2, T=60000):

        self.applyNoise = applyNoise
        self.bc = bc
        self.dt = dt  # time step
        self.grid_size = grid_size  # length of grid
        self.save_nth = save_nth
        self.ss_prec = ss_prec
        self.ssIteration = -1
        self.T = T  # wall time

        self.dx = sys_size/self.grid_size  # space step
        self.n = int(self.T/self.dt)

        self.alpha = alpha
        self.beta = beta
        self.dA = dA
        self.dP = dP
        self.kAP = kAP
        self.kPA = kPA
        self.konA = konA
        self.konP = konP
        self.koffA = koffA
        self.koffP = koffP
        self.Ptot = Ptot
        self.ratio = ratio
        self.StoV = StoV

        self.Atot = self.ratio * Ptot
        self.Acy = ones(int(ceil(self.n/self.save_nth)))
        self.Pcy = ones(int(ceil(self.n/self.save_nth)))

    def calc_step_change(self, plotting=False, startPlot=0):
        e = argmax(sum(self.A, 0) == self.grid_size) - 1
        self.stepDiff = [sum(self.A[:, i+1])-sum(self.A[:, i])
                         for i in range(0, e)]
        self.relStepChange = [abs(mean(self.A[:, i+1] / self.A[:, i]))
                              for i in range(0, e)]

        if plotting:
            plot(self.stepDiff[startPlot:])
            show()
            plot(self.relStepChange[startPlot:])
            show()

    def plot_steady_state(self):
        x = linspace(0, self.grid_size*self.dx, self.grid_size)
        if self.ssIteration == -1:
            plot(x, self.A[:, -1])
            plot(x, self.P[:, -1])            
        else:
            plot(x, self.A[:, int(self.ssIteration/self.save_nth)])
            plot(x, self.P[:, int(self.ssIteration/self.save_nth)])
        show()

    def set_init_profile(self):
        self.A = 2*ones((self.grid_size, int(ceil(self.n/self.save_nth))))*1.0
        self.P = ones((self.grid_size, int(ceil(self.n/self.save_nth))))*1.0

        if self.bc == 'PER':
            quarter = int(round(self.grid_size/4))
            self.A[quarter:int(round(self.grid_size/2)+quarter), 0] = 0
            self.P[0:quarter, 0] = 0
            self.P[-quarter:, 0] = 0
        elif self.bc == 'NEU':
            self.A[int(round(self.grid_size/2)):, 0] = 0
            self.P[0:int(round(self.grid_size/2)), 0] = 0
            # Straight gradient initial conditions:
            # self.P[:, 0] = linspace(0, 1, self.grid_size)
            # self.A[:, 0] = flipud(self.A[:, 0])
        if self.applyNoise:
            self.A[:, 0] = (random.normal(0, 1, 100) /
                            10*self.A[:, 0] + self.A[:, 0])
            self.P[:, 0] = (random.normal(0, 1, 100) /
                            10*self.P[:, 0] + self.P[:, 0])

    def show_movie(self, fname=None):
        # fig, ax = subplots()
        size_factor = 134.6/2/(self.grid_size*self.dx)
        fig = figure(figsize=(5, 4.0), facecolor='black')
        ax = gca()
        # Add a subplot with no frame
        # ax = subplot(111, frameon=False)
        # fig = figure(facecolor='black')
        # ax = gca()
        # fig.figsize = (1, 6)
        ax.set_xlim((0, self.grid_size*self.dx))
        ax.set_ylim((0, 5.2))
        line, = ax.plot([], [], lw=2)
        lineA, = ax.plot([], [], lw=2, color='red')
        lineP, = ax.plot([], [], lw=2, color='dodgerblue')
        t_text = ax.text(0.1*self.grid_size*self.dx,
                         0.2, r't = ' + str(0) + ' s',
                         fontsize=12, color='darkorange')
        ax.tick_params(axis='x', colors='darkorange', which='both')
        ax.tick_params(axis='y', colors='darkorange', which='both')
        xlabel(r'Position $[\mu m]$', color='darkorange', fontsize=12)
        ylabel(r'Concentration $[a.u.]$', color='darkorange', fontsize=12)
        ax.spines['left'].set_color('darkorange')
        ax.spines['bottom'].set_color('darkorange')
        ax.set_facecolor((0, 0, 0))

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
            t_text.set_text(r't = ' + str(round(i*self.dt*self.save_nth, 0)) +
                            ' s')
            return (line,)

        # call the animator.
        if self.ssIteration is -1:
            i_e = self.A.shape[1]
        else:
            i_e = argmax(sum(self.A, 0) == 2*self.grid_size)
        a = animation.FuncAnimation(fig, animate, init_func=init,
                                    frames=range(0, int(i_e), int(i_e/20)),
                                    interval=100, blit=False)
        if fname is None:
            a.save('lines.mp4', savefig_kwargs={'facecolor': 'black'}, dpi=200)
        else:
            a.save(fname + '.mp4', savefig_kwargs={'facecolor': 'black'},
                   dpi=500)

    def simulate(self):
        def check_steady_state():
            # Save every save_nth frame, check whether steady state reached
            if (i is not 0) and (i % self.save_nth) == 0:
                j = int(i/self.save_nth)
                self.A[:, j] = An
                self.P[:, j] = Pn
                self.Acy[j] = Acy
                self.Acy[j] = Pcy
                # Strict criteria
                if sum(self.A[:, j-1]-An) == 0:
                    # Append the next frame (An), to show that steady
                    # state was reached
                    self.A[:, j+1] = An
                    self.P[:, j+1] = Pn
                    self.Acy[j+1] = self.Atot-self.StoV*sum(An)/self.grid_size
                    self.Acy[j+1] = self.Ptot-self.StoV*sum(Pn)/self.grid_size
                    print('steady state reached')
                    self.ssIteration = i
                    return True
                # Not so strict criteria
                ind = int(i/self.save_nth)
                if ((i*self.dt > 600) and
                    (abs(mean(self.A[:, ind] /
                     self.A[:, ind-int(600/self.dt/self.save_nth)]-1)) <
                     self.ss_prec)):
                    # Append the next frame (An), to show that steady
                    # state was reached
                    self.A[:, j+1] = An
                    self.P[:, j+1] = Pn
                    self.Acy[j+1] = self.Atot-self.StoV*sum(An)/self.grid_size
                    self.Acy[j+1] = self.Ptot-self.StoV*sum(Pn)/self.grid_size
                    print('steady state not so strict reached.')
                    self.ssIteration = i
                    return True
            return False

        if self.bc == 'PER':
            def laplacian(Z):
                l = len(Z)
                Z = r_[Z, Z, Z]
                Zleft = Z[0:-2]
                Zright = Z[2:]
                Zcenter = Z[1:-1]
                LAP = (Zleft + Zright - 2*Zcenter) / self.dx**2
                return LAP[l-1:2*l-1]
            self.set_init_profile()
            An = copy(self.A[:, 0])
            Pn = copy(self.P[:, 0])
            for i in range(self.n-1):
                deltaA = laplacian(An)
                deltaP = laplacian(Pn)
                Acy = self.Atot - self.StoV*sum(An)/self.grid_size
                Pcy = self.Ptot - self.StoV*sum(Pn)/self.grid_size
                An = An+self.dt*(self.dA*deltaA - self.koffA*An +
                                 self.konA*Acy - self.kAP*An*Pn**self.alpha)
                Pn = Pn+self.dt*(self.dP*deltaP - self.koffP*Pn +
                                 self.konP*Pcy - self.kPA*An**self.beta*Pn)
                if check_steady_state() is True:
                    break

        elif self.bc == 'NEU':
            def laplacian(Z):
                Zleft = Z[0:-2]
                Zright = Z[2:]
                Zcenter = Z[1:-1]
                return (Zleft + Zright - 2*Zcenter) / self.dx**2
            self.set_init_profile()
            An = copy(self.A[:, 0])
            Pn = copy(self.P[:, 0])
            for i in range(self.n-1):
                deltaA = laplacian(An)
                deltaP = laplacian(Pn)
                Acy = self.Atot - self.StoV*sum(flipud(An), 0)/self.grid_size
                Pcy = self.Ptot - self.StoV*sum(Pn, 0)/self.grid_size
                # Defining Ra and Rp separately is necessary in order to not
                # update An before Pn (which would then have an effect on Pn
                # in the same iteration, making the system asymmetric)
                Rp = self.dt*(self.dP*deltaP-self.koffP*Pn[1:-1]+self.konP*Pcy -
                              self.kPA*An[1:-1]**self.beta*Pn[1:-1])
                Ra = self.dt*(self.dA*deltaA-self.koffA*An[1:-1]+self.konA*Acy -
                              self.kAP*Pn[1:-1]**self.alpha*An[1:-1])
                Pn[1:-1] = Pn[1:-1] + Rp
                An[1:-1] = An[1:-1] + Ra
                # Neumann conditions
                for Z in (An, Pn):
                    Z[0] = Z[1]
                    Z[-1] = Z[-2]
                if check_steady_state() is True:
                    break
        if self.ssIteration == -1:
            print('Steady state not reached within ' + str(self.ssIteration) +
                  ' iterations.')


class Sim_Container:

    def __init__(self, param_dict, no_workers=8, sys='symmetric'):
        self.param_dict = param_dict
        # self.r = r

        self.no_workers = no_workers
        self.sys = sys

    def init_simus(self):
        self.simList = []
        for k in self.param_dict['S']:
            if self.sys == 'symmetric':
                self.simList.append(ParSim(alpha=2, beta=2, dA=0.15, dP=0.15,
                                           dt=0.05, kAP=1, kPA=1, koffA=0.005,
                                           koffP=0.005, konA=0.006, konP=0.006,
                                           ratio=1.01, grid_size=200,
                                           ss_prec=0.000001, sys_size=k,
                                           T=120000, save_nth=400))
            elif self.sys == 'asymmetric':
                self.simList.append(ParSim(alpha=2, beta=2, dA=0.015, dP=0.15,
                                           dt=0.05, kAP=1, kPA=1, koffA=0.005,
                                           koffP=0.005, konA=0.06, konP=0.06,
                                           ratio=1.01, grid_size=200,
                                           ss_prec=0.000001, sys_size=k,
                                           T=120000, save_nth=400))

            elif self.sys == 'PARsys':
                self.simList.append(ParSim())

    def pickle_data(self, fname, postfix=''):
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
            # print(str(i.D))
            # print(str(i.koff))
            # print(str(i.dx*i.grid_size))
        show()

    def run_simus(self):
        # Create pool, use starmap to pass more than one parameter, do work
        # pool = Pool(processes=self.no_workers)
        # res = pool.map(sim_indiv, self.simList)
        # self.simus = res
        # pool.close()
        # pool.join()
        self.simus = [sim_indiv(self.simList[0])]


def s_to_v(mode, describers):
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


def msum(iterable):
    "Full precision summation using multiple floats for intermediate values"
    # Rounded x+y stored in hi with the round-off stored in lo.  Together
    # hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
    # to each partial so that the list of partial sums remains exact.
    # Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
    # www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps

    partials = []               # sorted, non-overlapping partial sums
    for x in iterable:
        i = 0
        for y in partials:
            if abs(x) < abs(y):
                x, y = y, x
            hi = x + y
            lo = y - (hi - x)
            if lo:
                partials[i] = lo
                i += 1
            x = hi
        partials[i:] = [x]
    return sum(partials, 0.0)
