#!/usr/bin/env python

# Data parser for PIXIE List Mode binary data *.bin
# The data contained in the .bin file is reformatted into an HDF5 file
# that stores event information is a series of related table entries for quick
# extraction of the necessary events used in spectra construction.

import struct
import time
import datetime
import os
import numpy as np
import math
from tables import *

# pyramds_cfg contains user-supplied info on detector systems and which data
# files are to be used in the analysis
from pyramds_cfg import *

file_counter = 1
file_path = (file_series + '%04d') % file_counter

# ******************************INFO READ-IN************************************
with open(file_path + '.ifm','rU') as finfo:
    # Read in the entire .ifm file as a series of lines. From info_str_list, the
    # necessary sections can be acquired for various data parameters. This works
    # granted the .ifm file never changes its format.
    info_str_list = finfo.readlines()
    
    date_str = info_str_list[1][23:-2]
    run_start_time = datetime.datetime(*time.strptime(date_str, \
                        '%I:%M:%S %p %a, %b %d, %Y')[0:6])
    
    # Total time = real time
    total_time = info_str_list[6].split()[3]
    
    # Live time for each detector channel
    live_time =[0]*4
    for channel in range(4):
        live_time[channel] = info_str_list[9 + channel].split()[2]
    
    bufheadlen = int(info_str_list[33].split()[1])
    eventheadlen = int(info_str_list[34].split()[1])
    chanheadlen = 2 #int(info_str_list[35].split()[1])

# Only initialize the buffer counter before the entire run... not each file.
buffer_no = 0

# ******************************DATA READ-IN************************************

# Setup class variables for each field in the .h5 table
class GammaEvent(IsDescription):
    energy_0    = Int16Col(pos=0)    # energy reading from channels 0, 1, 2
    energy_1    = Int16Col(pos=1)
    energy_2    = Int16Col(pos=2)
    deltaT_01   = Float32Col(pos=3)  # time difference between energy_0 & _1
    deltaT_02   = Float32Col(pos=4)  # and so on...
    deltaT_12   = Float32Col(pos=5)  # and so on...
    timestamp   = Float32Col(pos=6)

h5file = openFile(file_series + '.h5', mode = 'w', title = 'Data - ' + file_series)
group = h5file.createGroup(h5file.root, 'bin_data_parse', 'PIXIE Binary Parse')

# This table is where the data will be placed after unpacking it from binary
table = h5file.createTable(group, 'readout', GammaEvent, "Data readout")

# Pointer object to place values on in the row for each event
event = table.row

while os.path.exists(file_path + '.bin'):
    
    file_pos = 0
    
    with open(file_path + '.bin', 'rb') as fin:
        
        print('Working on ' + file_path + '.bin ...')
        starttime = time.time()
        
        # Pointer object to place values on in the row for each event. Must
        # create new instance each time .flush is called (every .bin file)
        event = table.row
        
        word = fin.read(2)
        
        while word:
        
            # Determine the length of the buffer to be used in recognizing new buffer
            buf_ndata = struct.unpack('<H', word)[0]
            fin.seek(-2,1)
            file_pos += (buf_ndata * 2) # Increase limits of the file pointer
        
            # Read in buffer header data
            header = fin.read(bufheadlen * 2)
            buf_ndata, buf_modnum, buf_format, buf_timehi, buf_timemi, buf_timelo = \
                                    struct.unpack('<' + str(bufheadlen) + 'H', header)
            
            while fin.tell() < file_pos:
                # Read in event header data
                header = fin.read(eventheadlen * 2)
                evt_pattern, evt_timehi, evt_timelo = \
                                struct.unpack('<' + str(eventheadlen) + 'H', header)
                
                read_pattern = list(map(int, bin(evt_pattern)[-4:]))
                read_pattern.reverse()
                hit_pattern = ''.join(map(str, read_pattern))
                
                trigger_vals = [float('nan')]*3
                for channel in range(3):
                    if read_pattern[channel] == 1:
                        words = fin.read(2 * 2)
                        chan_trigtime, energy = \
                                    struct.unpack('<' +str(chanheadlen) + 'H', words)
                        trigger_vals[channel] = (float((evt_timehi * 64000 + chan_trigtime) * tunits))
                        
                        # Store the data read in from the binary file into the
                        # HDF5 table.
                        if energy <= energy_max:
                            event['energy_' + str(channel)] = energy
                        else: event['energy_' + str(channel)] = -1
                    elif read_pattern[channel] == 0:
                        event['energy_' + str(channel)] = -1
                
                event['deltaT_01'] = abs(trigger_vals[0] - trigger_vals[1])
                event['deltaT_02'] = abs(trigger_vals[0] - trigger_vals[2])
                event['deltaT_12'] = abs(trigger_vals[1] - trigger_vals[2])
                
                event['timestamp'] = (buf_timehi * 64000 * 64000 + evt_timehi * 64000 + evt_timelo) * tunits * 1e-9
                
                event.append()
        
            # Read word, buf_ndata, to continue loop or break
            word = fin.read(2)    
        
            buffer_no += 1
            if buffer_no%100 == 0:
                print('\rBuffer No. {0}'.format(buffer_no))
            
            # Flush data to the HFD5 table and start new buffer    
            table.flush()
            
    file_counter += 1
    file_path = (file_series + '%04d') % file_counter
    
# Define new groups within HDF5 for storing each spectra type for each detector
spectra = h5file.createGroup(h5file.root, "spectra", "Spectra Arrays")
gNormal = h5file.createGroup(h5file.root.spectra, "normal", "Normal Data")
gCompton = h5file.createGroup(h5file.root.spectra, "compton", "Compton-Supp Data")
gGGcoinc = h5file.createGroup(h5file.root.spectra, "ggcoinc", "Gamma-Gamma Data")

gN1 = h5file.createGroup(h5file.root.spectra.normal, "det1", "Det 1 Normal Data")
gN1_T = h5file.createGroup(h5file.root.spectra.normal.det1, "t_arrays", "Time array - Det 1")
gN2 = h5file.createGroup(h5file.root.spectra.normal, "det2", "Det 2 Normal Data")
gN2_T = h5file.createGroup(h5file.root.spectra.normal.det2, "t_arrays", "Time array - Det 2")

gC1 = h5file.createGroup(h5file.root.spectra.compton, "det1", "Det 1 Compton-Supp Data")
gC1_T = h5file.createGroup(h5file.root.spectra.compton.det1, "t_arrays", "Time array - Det 1")
gC2 = h5file.createGroup(h5file.root.spectra.compton, "det2", "Det 2 Compton-Supp Data")
gC2_T = h5file.createGroup(h5file.root.spectra.compton.det2, "t_arrays", "Time array - Det 2")

gGG1 = h5file.createGroup(h5file.root.spectra.ggcoinc, "det1", "Det 1 Gamma-Gamma Data")
gGG1_T = h5file.createGroup(h5file.root.spectra.ggcoinc.det1, "t_arrays", "Time array - Det 1")
gGG2 = h5file.createGroup(h5file.root.spectra.ggcoinc, "det2", "Det 2 Gamma-Gamma Data")
gGG2_T = h5file.createGroup(h5file.root.spectra.ggcoinc.det2, "t_arrays", "Time array - Det 2")

################################################################################
# Store Det-1 & 2 arrays containing normal counts###############################
################################################################################

norm_evts = [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_1 != -1) | (energy_2 != -1)""") ]

h5file.createArray(gNormal, 'norm_evts12', np.array(norm_evts), "Normal Timestamped Events - Det 1")

# Det 1

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[0] for row in norm_evts if row[0] != -1]:
    spec_temp[x[0]] +=1
    
h5file.createArray(gN1, 'norm1_spec', spec_temp, "Normal Spec Array - Det 1")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[2] for i in norm_evts if (int(i[0]) == chn)])
    h5file.createArray(gN1_T, 'chan' + chn, t_array_dump)

# Det 2

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[0] for row in norm_evts if row[0] != -1]:
    spec_temp[x[0]] +=1
    
h5file.createArray(gN2, 'norm2_spec', spec_temp, "Normal Spec Array - Det 2")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[1] for i in norm_evts if (int(i[1]) == chn)])
    h5file.createArray(gN2_T, 'chan' + chn, t_array_dump)
    
################################################################################
# Store Det-1 & 2 arrays containing Compton counts##############################
################################################################################

compt_evts = [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 == -1)""") ]
compt_evts += [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 == -1) & (energy_2 != -1)""") ]

h5file.createArray(gCompton, 'compt_evts12', np.array(compt_evts), "Compton-Supp Timestamped Events - Det 1/2")

# Det 1

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[0] for row in compt_evts if row[0] != -1]:
    spec_temp[x] += 1
   
h5file.createArray(gC1, 'compt1_spec', spec_temp, "Compton-Supp Spec Array - Det 1")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[2] for i in compt_evts if (int(i[0]) == chn)])
    h5file.createArray(gC1_T, 'chan' + chn, t_array_dump)

# Det 2

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[1] for row in compt_evts if row[1] != -1]:
    spec_temp[x] += 1

h5file.createArray(gC2, 'compt2_spec', spec_temp, "Compton-Supp Spec Array - Det 2")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[2] for i in compt_evts if (int(i[1]) == chn)])
    h5file.createArray(gC2_T, 'chan' + chn, t_array_dump)

################################################################################
# Store Det 1 & 2 Gamma-Gamma arrays containg counts############################
################################################################################

gg_evts = [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 != -1) & (deltaT_12 < 90)""") ]

h5file.createArray(gGGcoinc, 'gg_evts', np.array(gg_evts), "Gamma-Gamma Timestamped Events - Det 1/2")

# Det 1

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in gg_evts:
    spec_temp[x[0]] += 1

h5file.createArray(gGG1, 'gg1_spec', np.array(spec_temp), "Gamma-Gamma Spec Array - Det 1")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[2] for i in gg_evts if (int(i[0]) == chn)])
    h5file.createArray(gGG1_T, 'chan' + chn, t_array_dump)

# Det 2

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in gg_evts:
    spec_temp[x[1]] += 1
    
h5file.createArray(gGG2, 'gg2_spec', np.array(spec_temp), "Gamma-Gamma Spec Array - Det 2")

for chn in range(energy_max + 1):
    t_array_dump = np.array([i[2] for i in gg_evts if (int(i[1]) == chn)])
    h5file.createArray(gGG2_T, 'chan' + chn, t_array_dump)

# Specgram coincidence array for Det 1/2
#x = [row[0] for row in f.root.spectra.gg_evts]
#y = [row[0] for row in f.root.spectra.gg_evts]
#plt.hexbin(x,y)    
    
h5file.close()