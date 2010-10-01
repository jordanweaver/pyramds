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
            
            # Flush data to the HFD5 table and start new buffer    
            table.flush()
            
    file_counter += 1
    file_path = (file_series + '%04d') % file_counter
    
# Define new group within HDF5 for storing each spectra type
spectra = h5file.createGroup(h5file.root, "spectra", "Spectra Arrays")

################################################################################
# Store Det-1 arrays containing normal counts###################################
################################################################################
norm_evts_1 = [[row['energy_1'],row['timestamp']] for row in table.where("""(energy_1 != -1)""") ]

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in norm_evts_1:
    spec_temp[x[0]] +=1

h5file.createArray(spectra, 'norm1_evts', np.array(norm_evts_1), "Normal Timestamped Events - Det 1")
h5file.createArray(spectra, 'norm1_spec', spec_temp, "Normal Spec Array - Det 1")

################################################################################
# Store Det-2 arrays containing normal counts###################################
################################################################################
norm_evts_2 = [[row['energy_2'],row['timestamp']] for row in table.where("""(energy_2 != -1)""") ]

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in norm_evts_2:
    spec_temp[x[0]] +=1
    
h5file.createArray(spectra, 'norm2_evts', np.array(norm_evts_2), "Normal Timestamped Events - Det 2")
h5file.createArray(spectra, 'norm2_spec', spec_temp, "Normal Spec Array - Det 2")

################################################################################
# Store Det-1 & 2 arrays containing Compton counts##############################
################################################################################
compt_evts = [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 == -1)""") ]
compt_evts += [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 == -1) & (energy_2 != -1)""") ]

h5file.createArray(spectra, 'compt12_evts', np.array(compt_evts), "Compton Timestamped Events - Det 1/2")

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[0] for row in compt_evts if row[0] != -1]:
    spec_temp[x] += 1
    
h5file.createArray(spectra, 'compt1_spec', spec_temp, "Compton Spec Array - Det 1")

spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in [row[1] for row in compt_evts if row[1] != -1]:
    spec_temp[x] += 1

h5file.createArray(spectra, 'compt2_spec', spec_temp, "Compton Spec Array - Det 2")

################################################################################
# Store Det 1 & 2 Gamma-Gamma arrays containg counts############################
################################################################################
gg_evts = [[row['energy_1'],row['energy_2'],row['timestamp']]
    for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 != -1) & (deltaT_12 < 90)""") ]

h5file.createArray(spectra, 'gg_evts', np.array(gg_evts), "Gamma-Gamma Timestamped Events - Det 1/2")

# Spectrum array for Det 1
spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in gg_evts:
    spec_temp[x[0]] +=1

h5file.createArray(spectra, 'gg1_spec', np.array(spec_temp), "Gamma-Gamma Spec Array - Det 1")

# Spectrum array for Det 2
spec_temp = np.zeros(energy_max + 1, dtype=np.int32)
for x in gg_evts:
    spec_temp[x[1]] +=1
    
h5file.createArray(spectra, 'gg2_spec', np.array(spec_temp), "Gamma-Gamma Spec Array - Det 2")

# Specgram coincidence array for Det 1/2
#x = [row[0] for row in f.root.spectra.gg_evts]
#y = [row[0] for row in f.root.spectra.gg_evts]
#plt.hexbin(x,y)    
    
h5file.close()