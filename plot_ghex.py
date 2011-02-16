#!/usr/bin/env python

import numpy as np
import datetime
import time
import tables as tb
import matplotlib.cm as cm
import  matplotlib.pyplot as plt

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

hexplot = plt.hexbin(gg_list_1, gg_list_2, cmap=cm.jet, gridsize=200)
plt.axis([xmin, xmax, ymin, ymax])
plt.title(file_path)
cb = plt.colorbar()
cb.set_label('Counts')

plt.show()