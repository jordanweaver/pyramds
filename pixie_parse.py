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
    
class AggEvent1(IsDescription):
    energy    = Int16Col(pos=1)
    timestamp   = Float32Col(pos=6)
    
class AggEvent2(IsDescription):
    energy_1    = Int16Col(pos=1)
    energy_2    = Int16Col(pos=2)
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
                
            # Remember the time of the first buffer of entire run. Use this for comparing time stops.
            if buffer_no == 0:
                t_start_hi = buf_timehi * 64000 * 64000
                t_start_mi = buf_timemi * 64000
                t_start_lo = buf_timelo
                t_start = (t_start_hi + t_start_mi + t_start_lo) * tunits * 1e-9 # in seconds
            
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
    
t_final = (buf_timehi * 64000 * 64000 + evt_timehi * 64000 + evt_timelo) * tunits * 1e-9 # in seconds
t_duration = t_final - t_start
t_array_dim = int(np.ceil(t_duration/t_steps))
    
# Define new groups within HDF5 for storing each spectra type for each detector
spectra = h5file.createGroup(h5file.root, "spectra", "Spectra Arrays")
gNormal = h5file.createGroup(h5file.root.spectra, "normal", "Normal Data")
gCompton = h5file.createGroup(h5file.root.spectra, "compton", "Compton-Supp Data")
gGGcoinc = h5file.createGroup(h5file.root.spectra, "ggcoinc", "Gamma-Gamma Data")

gN1 = h5file.createGroup(h5file.root.spectra.normal, "det1", "Det 1 Normal Data")
gN2 = h5file.createGroup(h5file.root.spectra.normal, "det2", "Det 2 Normal Data")

gC1 = h5file.createGroup(h5file.root.spectra.compton, "det1", "Det 1 Compton-Supp Data")
gC2 = h5file.createGroup(h5file.root.spectra.compton, "det2", "Det 2 Compton-Supp Data")

gGG1 = h5file.createGroup(h5file.root.spectra.ggcoinc, "det1", "Det 1 Gamma-Gamma Data")
gGG2 = h5file.createGroup(h5file.root.spectra.ggcoinc, "det2", "Det 2 Gamma-Gamma Data")

################################################################################
# Store Det-1 & 2 arrays containing normal counts###############################
################################################################################

norm12table = h5file.createTable(gNormal, 'norm_evts12', AggEvent2, "Normal Time-Stamped Data")
norm_event = norm12table.row

for row in table.where("""(energy_1 != -1) | (energy_2 != -1)"""):
    norm_event['energy_1'] = row['energy_1']
    norm_event['energy_2'] = row['energy_2']
    norm_event['timestamp'] = row['timestamp']
    norm_event.append()
    norm12table.flush()

# Det 1
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in norm12table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy_1'] >= 0:
        dt_temp[row['energy_1']] += 1
dt_array[-1] = dt_temp.copy()

    
h5file.createArray(gN1, 'norm1_spec', dt_array, "Normal Time-Chunked Spec Array - Det 1")

# Det 2
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in norm12table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy_2'] >= 0:
        dt_temp[row['energy_2']] += 1
dt_array[-1] = dt_temp.copy()
        
h5file.createArray(gN2, 'norm2_spec', dt_array, "Normal Time-Chunked Spec Array - Det 2")

################################################################################
# Store Det-1 & 2 arrays containing Compton counts##############################
################################################################################

compt1table = h5file.createTable(gCompton, 'compt_evts1', AggEvent1, "Compton Time-Stamped Data - Det 1")
compt_event = compt1table.row

for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 == -1)"""):
    compt_event['energy'] = row['energy_1']
    compt_event['timestamp'] = row['timestamp']
    compt_event.append()
    compt1table.flush()
    
compt2table = h5file.createTable(gCompton, 'compt_evts2', AggEvent1, "Compton Time-Stamped Data - Det 2")
compt_event = compt2table.row

for row in table.where("""(energy_0 == -1) & (energy_1 == -1) & (energy_2 != -1)"""):
    compt_event['energy'] = row['energy_2']
    compt_event['timestamp'] = row['timestamp']
    compt_event.append()
    compt2table.flush()
    
# Det 1
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in compt1table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy'] >= 0:
        dt_temp[row['energy']] += 1
dt_array[-1] = dt_temp.copy()
   
h5file.createArray(gC1, 'compt1_spec', dt_array, "Compton-Supp Time-Chunked Spec Array - Det 1")

# Det 2
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in compt2table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy'] >= 0:
        dt_temp[row['energy']] += 1
dt_array[-1] = dt_temp.copy()
   
h5file.createArray(gC2, 'compt2_spec', dt_array, "Compton-Supp Time-Chunked Spec Array - Det 2")

################################################################################
# Store Det 1 & 2 Gamma-Gamma arrays containg counts############################
################################################################################

gg12table = h5file.createTable(gGGcoinc, 'gg_evts12', AggEvent2, "G-G Time-Stamped Data")
gg_event = gg12table.row

for row in table.where("""(energy_0 == -1) & (energy_1 != -1) & (energy_2 != -1) & (deltaT_12 < %f)"""%short_window):
    gg_event['energy_1'] = row['energy_1']
    gg_event['energy_2'] = row['energy_2']
    gg_event['timestamp'] = row['timestamp']
    gg_event.append()
    gg12table.flush()

# Det 1
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in gg12table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy_1'] >= 0:
        dt_temp[row['energy_1']] += 1
dt_array[-1] = dt_temp.copy()

h5file.createArray(gGG1, 'gg1_spec', dt_array, "G-G Time-Chunked Spec Array - Det 1")

# Det 2
dt_temp = np.zeros(energy_max + 1, dtype=np.int32)
dt_array = np.zeros((t_array_dim + 1, energy_max + 1), dtype=np.int32)
ti = 0
for row in gg12table:
    if (row['timestamp'] - t_start) >= (ti + 1) * t_steps:
        ti += 1
        dt_array[ti] = dt_temp.copy()
    if row['energy_2'] >= 0:
        dt_temp[row['energy_2']] += 1
dt_array[-1] = dt_temp.copy()

h5file.createArray(gGG1, 'gg2_spec', dt_array, "G-G Time-Chunked Spec Array - Det 2")
    
h5file.close()