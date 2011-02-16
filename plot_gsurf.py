#!/usr/bin/env python

import numpy as np
import datetime
import time
import tables as tb
import matplotlib.cm as cm
import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#file_path = raw_input('Enter path to HDF file: ')
#file_path = '/python/PIXIE_runs/PreRevision/debug-.h5'
#file_path = '/python/PIXIE_runs/Standards/Co60_LOW-5min.h5'
#file_path = '/python/PIXIE_runs/Se_Test/Se_test.h5'
file_path = '/python/PIXIE_runs/WFP_1222/WFP2_1222-3d.h5'

f = tb.openFile(file_path, 'r')
evtstab = f.root.bin_data_parse.readout

emax = 700
gg_list_1 = np.array([row['energy_1'] for row in evtstab.where("""(energy_1 <= emax) & (energy_2 <= emax) & (energy_1 > 0) & (energy_2 > 0)""")])
gg_list_2 = np.array([row['energy_2'] for row in evtstab.where("""(energy_1 <= emax) & (energy_2 <= emax) & (energy_1 > 0) & (energy_2 > 0)""")])

xmin = 0
xmax = gg_list_1.max()
ymin = 0
ymax = gg_list_2.max()

print('Generating cust_list...')
cust_hist = np.ones(shape=(emax,emax), dtype=np.int32)
for i in range(len(gg_list_1)):
    if ((gg_list_1[i]<emax) and (gg_list_2[i]<emax)):
        cust_hist[gg_list_1[i],gg_list_2[i]] += 1

xi = np.arange(emax)
yi = np.arange(emax)
X,Y = np.meshgrid(xi,yi)

print('Initializing plot...')
fig = plt.figure()
ax = fig.gca(projection='3d')

print('Starting plot...')
ax.plot_surface(X,Y,cust_hist, linewidth=0, rstride=2, cstride=2, antialiased=False, cmap=cm.jet)
plt.show()
