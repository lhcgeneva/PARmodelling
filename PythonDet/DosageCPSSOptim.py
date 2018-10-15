import itertools
import sys
import numpy as np
from ParDetAsymm import ParSim
import time

sz = np.linspace(3, 30, 8)

t0 = time.time()
symmetric = {'alpha': 2, 'beta': 2, 'dA': 0.1, 'dP': 0.1,
             'kAP': 1, 'kPA': 1, 'koffA': 0.005, 'koffP': 0.005,
             'konA': 0.006, 'konP': 0.006, 'Ptot': 1, 'ratio': 1.0,
             'sys_size': 30}
sym = [copy.deepcopy(symmetric) for _ in sz]
for i, dic in enumerate(sym):
    sym[i]['sys_size'] = sz[i]
ratio_below = [0.01 for _ in sym]
ratio_above = [5 for _ in sym]
ratio = [0.01 for _ in sym]
outcome = [3 for _ in sym]
for j in range(10):
    for i in range(len(sym)):
        sym[i]['ratio'] = ratio[i]
    s = Sim_Container(sym, no_workers=8)
    s.init_simus()
    s.run_simus()
#     p.pickle_data(fname='ratio_' + str(p.ratio) + 'size_' + str(p.sys_size))
    ratio_temp = copy.deepcopy(ratio)
    for i in range(len(sym)):
        outcome[i] = s.simus[i].finished_in_time
        # 0 - unpolarized P won
        # 1 - unpolarized A won
        # 2 - polarized
        if (outcome[i] == 2) or (outcome[i] == 0):
            ratio[i] = ratio[i] + (ratio_above[i]-ratio[i])/2
            ratio_below[i] = ratio_temp[i]
        elif outcome[i] == 1:
            ratio[i] = ratio[i] - (ratio[i]-ratio_below[i])/2
            ratio_above[i] = ratio_temp[i]
    print(ratio)
#         if abs(ratio_above[i]-ratio_below[i])/ratio_below[i] < 0.002:
#             break