#!/usr/bin/env python

import numpy as np
import datetime
import time
import tables as tb
import matplotlib.cm as cm
import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from function_lib import marker2energy

#file_path = raw_input('Enter path to HDF file: ')
#file_path = '/python/PIXIE_runs/PreRevision/debug-.h5'
#file_path = '/python/PIXIE_runs/Standards/Co60_LOW-5min.h5'
#file_path = '/python/PIXIE_runs/Se_Test/Se_test.h5'
file_path = '/python/PIXIE_runs/WFP_1222/WFP2_1222-3d.h5'

en_coeff1 = [0.118379, 0.41] #3-day coeff
en_coeff2 = [-1.646097, 0.47]

#en_coeff1 = [-0.931608, 0.391917] #28-day coeff
#en_coeff2 = [-1.646097, 0.414058]

f = tb.openFile(file_path, 'r')
evtstab = f.root.bin_data_parse.readout

emax = 800
gg_list_1 = np.array([row['energy_1'] for row in evtstab.where("""(energy_1 <= emax) & (energy_2 <= emax) & (energy_1 > 0) & (energy_2 > 0)""")])
gg_list_2 = np.array([row['energy_2'] for row in evtstab.where("""(energy_1 <= emax) & (energy_2 <= emax) & (energy_1 > 0) & (energy_2 > 0)""")])

gg_en_1 = np.zeros(len(gg_list_1))
gg_en_2 = np.zeros(len(gg_list_2))

for i in range(len(gg_list_1)):
    gg_en_1[i] = marker2energy(gg_list_1[i], en_coeff1)
    
for i in range(len(gg_list_2)):
    gg_en_2[i] = marker2energy(gg_list_2[i], en_coeff2)

cust_hist = np.ones(shape=(emax,emax), dtype=np.int32)
for i in range(len(gg_list_1)):
    if ((gg_list_1[i]<emax) and (gg_list_2[i]<emax)):
        cust_hist[gg_list_1[i],gg_list_2[i]] += 1

xi = np.arange(emax)
yi = np.arange(emax)

xcoord = en_coeff1[0] + en_coeff1[1]*xi
ycoord = en_coeff2[0] + en_coeff2[1]*yi

X,Y = np.meshgrid(xcoord,ycoord)

fig = plt.figure(figsize=(15,12))
ax = fig.gca(projection='3d')

ax.plot_surface(X,Y,cust_hist, linewidth=0, rstride=5, cstride=5, antialiased=False, cmap=cm.jet)
ax.set_xlabel('Lower Detector Energy (keV)')
ax.set_ylabel('Upper Detector Energy (keV)')
ax.set_zlabel('Counts')
plt.show()
