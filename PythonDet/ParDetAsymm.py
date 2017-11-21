from IPython.core.debugger import Tracer
from IPython.display import HTML
from itertools import repeat
from numpy import (abs, argmax, array_equal, ceil, copy, flipud, float64, linspace,
                   maximum, mean, ones, r_, random, round, sqrt, sum)
from matplotlib.pyplot import (cla, figure, gca, gcf, Normalize, plot, show,
                               subplots, subplot, xlabel, ylabel)
from matplotlib import animation, rc
from multiprocessing import Pool
from pickle import dump, HIGHEST_PROTOCOL
from scipy import arcsin
from scipy.special import ellipe


class ParSim(object):

    def __init__(self, applyNoise=False, alpha=1, bc='NEU', beta=2,
                 dA=0.28, dP=0.15, dt=0.01, grid_size=100, kAP=0.19, kPA=2,
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
        self.A = ones((self.grid_size, int(ceil(self.n/self.save_nth))))*1.0
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
            i_e = int(self.ssIteration/self.save_nth)
        # Tracer()()
        a = animation.FuncAnimation(fig, animate, init_func=init,
                                    frames=range(0, int(i_e), 10),
                                    interval=100, blit=False)
        if fname is None:
            a.save('lines.mp4', savefig_kwargs={'facecolor': 'black'}, dpi=200)
        else:
            a.save(fname + '.mp4', savefig_kwargs={'facecolor': 'black'},
                   dpi=500)

    def simulate(self):

        # Adaptive step size parameters
        atol = 0.000001
        rtol = 0.000001
        # 5TH ORDER RK COEFFICIENTS for Dormand-Prince
        a21 = 1/5
        a31 = 3/40
        a32 = 9/40
        a41 = 44/45
        a42 = -56/15
        a43 = 32/9
        a51 = 19372/6561
        a52 = -25360/2187
        a53 = 64448/6561
        a54 = -212/729
        a61 = 9017/3168
        a62 = -355/33
        a63 = 46732/5247
        a64 = 49/176
        a65 = -5103/18656
        a71 = 35/384
        a72 = 0
        a73 = 500/1113
        a74 = 125/192
        a75 = -2187/6784
        a76 = 11/84

        b1 = 35/384
        b2 = 0
        b3 = 500/1113
        b4 = 125/192
        b5 = -2187/6784
        b6 = 11/84
        b7 = 0

        bs1 = 5179/57600
        bs2 = 0
        bs3 = 7571/16695
        bs4 = 393/640
        bs5 = -92097/339200
        bs6 = 187/2100
        bs7 = 1/40

        c2 = 1/5
        c3 = 3/10
        c4 = 4/5
        c5 = 8/9
        c6 = 1
        c7 = 1

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
            self.set_init_profile()
            An = copy(self.A[:, 0])
            Pn = copy(self.P[:, 0])
            for i in range(self.n-1):
                delA = laplacianPER(An, self.dx)
                delP = laplacianPER(Pn, self.dx)
                Acy = self.Atot - self.StoV*sum(An)/self.grid_size
                Pcy = self.Ptot - self.StoV*sum(Pn)/self.grid_size
                An = An+self.dt*(self.dA*delA - self.koffA*An +
                                 self.konA*Acy - self.kAP*An*Pn**self.alpha)
                Pn = Pn+self.dt*(self.dP*delP - self.koffP*Pn +
                                 self.konP*Pcy - self.kPA*An**self.beta*Pn)
                if check_steady_state() is True:
                    break

        elif self.bc == 'NEU':
            def neu(An, Pn):
                delA = laplacianNEU(An, self.dx)
                delP = laplacianNEU(Pn, self.dx)
                Acy = self.Atot - self.StoV*sum(flipud(An), 0)/self.grid_size
                Pcy = self.Ptot - self.StoV*sum(Pn, 0)/self.grid_size
                # Defining Ra and Rp separately is necessary in order to not
                # update An before Pn (which would then have an effect on Pn
                # in the same iteration, making the system asymmetric)
                Rp = self.dt*(self.dP*delP-self.koffP*Pn[1:-1]+self.konP*Pcy -
                              self.kPA*An[1:-1]**self.beta*Pn[1:-1])
                # print(sum(Rp))
                Ra = self.dt*(self.dA*delA-self.koffA*An[1:-1]+self.konA*Acy -
                              self.kAP*Pn[1:-1]**self.alpha*An[1:-1])
                # Return arrays with first and last two elements equal
                # respecitvely, to impose zero derivatives
                return r_[Ra[0], Ra, Ra[-1]], r_[Rp[0], Rp, Rp[-1]]

            self.set_init_profile()
            A0 = self.A[:, 0]
            P0 = self.P[:, 0]
            self.t = [0]
            self.totalerror = [0]
            self.errRatio = []

            while self.t[-1] < self.T:
                # Calculate increments for RK45
                A1, P1 = neu(A0, P0)
                A2, P2 = neu(A0+a21*A1, P0+a21*P1)
                A3, P3 = neu(A0+a31*A1+a32*A2, P0+a31*P1+a32*P2)
                A4, P4 = neu(A0+a41*A1+a42*A2+a43*A3, P0+a41*P1+a42*P2+a43*P3)
                A5, P5 = neu(A0+a51*A1+a52*A2+a53*A3+a54*A4,
                             P0+a51*P1+a52*P2+a53*P3+a54*P4)
                A6, P6 = neu(A0+a61*A1+a62*A2+a63*A3+a64*A4+a65*A5,
                             P0+a61*P1+a62*P2+a63*P3+a64*P4+a65*P5)
                A7, P7 = neu(A0+a71*A1+a73*A3+a74*A4+a75*A5+a76*A6,  # a72=0
                             P0+a71*P1+a73*P3+a74*P4+a75*P5+a76*P6)
                # Update concentrations using A1-A6 and P1-P6, coefficient for
                # A7 and P7 is 0.
                An_new = A0+b1*A1+b3*A3+b4*A4+b5*A5+b6*A6  # b2/7=0 => no A2/7
                Pn_new = P0+b1*P1+b3*P3+b4*P4+b5*P5+b6*P6  # b2/7=0 => no P2/7
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
                dtnew = 0.99*self.dt*abs(1/totalerror)**(1/5)
                # Upper and lower bound for timestep to avoid changing too fast
                if dtnew > 5*self.dt:
                    dtnew = 5*self.dt
                elif dtnew < self.dt/5:
                    dtnew = self.dt/5
                # Set timestep for next round
                self.dt = dtnew
                # Accept step if error is on the order of error scale or below
                if totalerror < 1:
                    self.t.append(self.t[-1]+self.dt)
                    self.errRatio.append(max(A0/An_new))
                    self.totalerror.append(totalerror)
                    A0 = An_new
                    P0 = Pn_new

                # Break if things change by less than 0.1% over
                # the course of 1 min.
                if self.errRatio[-1]**(60/dtnew) < 1.001:
                    break

            self.A = An_new
            self.P = Pn_new
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
                                           ratio=1.01, grid_size=100,
                                           ss_prec=0.000001, sys_size=k,
                                           T=1200, save_nth=400))
            elif self.sys == 'asymmetric':
                self.simList.append(ParSim(alpha=2, beta=2, dA=k, dP=0.35,
                                           dt=0.01, kAP=1, kPA=1, koffA=0.005,
                                           koffP=0.005, konA=0.006, konP=0.006,
                                           ratio=1.00, grid_size=100,
                                           ss_prec=0.0001, sys_size=40.0,
                                           T=100000, save_nth=1000))            
            elif self.sys == 'Tom':
                self.simList.append(ParSim(alpha=2, beta=2, dA=1, dP=1,
                                           dt=0.001, kAP=1.3, kPA=1.3,
                                           koffA=0.3, koffP=0.3, konA=1,
                                           konP=1, ratio=1.00, grid_size=100,
                                           ss_prec=0.001, StoV=0.3,
                                           sys_size=50.0, T=100000,
                                           save_nth=1000))
            elif self.sys == 'TomHill':
                self.simList.append(ParSim(alpha=2, beta=2, dA=1, dP=1,
                                           dt=0.005, kAP=100*0.25,
                                           kPA=100*0.25, koffA=0.3,
                                           koffP=0.3, konA=1, konP=1,
                                           ratio=1.0, grid_size=100,
                                           ss_prec=0.001, sys_size=100,
                                           T=10000, save_nth=1000))   
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


def laplacianNEU(Z, dx):
    Zleft = Z[0:-2]
    Zright = Z[2:]
    Zcenter = Z[1:-1]
    return (Zleft + Zright - 2*Zcenter) / dx**2


def laplacianPER(Z, dx):
    l = len(Z)
    Z = r_[Z, Z, Z]
    Zleft = Z[0:-2]
    Zright = Z[2:]
    Zcenter = Z[1:-1]
    LAP = (Zleft + Zright - 2*Zcenter) / dx**2
    return LAP[l-1:2*l-1]


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
