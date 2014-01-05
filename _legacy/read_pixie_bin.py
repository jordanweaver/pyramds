#!/usr/bin/env python

# Data parser for PIXIE List Mode binary data *.bin
# Used to read in the run info from .ifm file corresponding to the .bin file
import struct
import time
import datetime
import os
import numpy as np

import pyramds_config
from funlib import build_gamma, write_spec

with open('Input.txt', 'r') as finput:
    
    en_coeff = {}
    fwhm_coeff = {}
    for c,channel in enumerate(detector_choices):
        en_coeff[channel] = map(float, enerfit[channel].split())
        fwhm_coeff[channel] = map(float, fwhmfit[channel].split())
    
    # Run build_gamma to get gamma_lib, sig_lookup
    build_gamma(lib_name, fwhm_coeff, en_coeff)

    gate_choices = {}
    gate_option = raw_input('Would you like gamma gating? (y/n) ')
    
    if gate_option == 'y':
        for i,dc in enumerate(detector_choices):
            gate_choices[dc] = []
            for j,ht in enumerate(spe_choices[dc]):
                gates = []
                gate_add = raw_input('Enter desired gamma energy gates (use ID from lib file) for Detector-' \
                                     + dc + ' Spectra-' + str(j+1) + ': ').split()
                for g,sig in enumerate(gate_add):
                    gates.append(sig_lookup[gate_add[g]])
                gate_choices[dc].append(list(gates))
        

roi_list = dict([(x, []) for x in detector_choices])
roi_user_input = raw_input('Enter signature IDs for desired ROIs (separated by spaces): ').split()
for i,dc in enumerate(detector_choices):
    for z,zaid in enumerate(roi_user_input):
        roi_list[dc].append(sig_lookup[int(dc)][zaid])

detector_choices = map(int, detector_choices)
    
time_stops = [300, 600, 900] # units of seconds
stop_counter = 0
    
file_counter = 1
file_path = (file_series + '%04d') % file_counter

tunits = 1000/75 # Nanoseconds
energy_max = 8192

# This array is for the storage of eventual output by the code. The code eventually
# writes .Spe files for viewing in a third party software and this array is where
# the information is stored until it is written. The dimension of '4' is for each
# of the four channels recorded by the PIXIE. 'max_spec_num' keeps tabs on the
# user input for the maximum number of eventual spectra that will be required by
# the user. Finally, 'energy_max' is the largest value of the "energy" (channel
# number in this sense) that we wish to store.
spec_markers = np.zeros( (4, int(max_spec_num), energy_max + 1), dtype=np.int32)
spec_stop_copy = np.zeros( (len(time_stops), 4, int(max_spec_num), energy_max + 1), dtype=np.int32)

# ******************************INFO READ-IN***********************************
with open(file_path + '.ifm','rU') as finfo:
    # Read in the entire .ifm file as a series of lines. From info_str_list, the
    # necessary sections can be acquired for various data parameters. This works
    # granted the .ifm file never changes its format.
    info_str_list = finfo.readlines()
    
    date_str = info_str_list[1][23:-2]
    run_start_time = datetime.datetime(*time.strptime(date_str, '%I:%M:%S %p %a, %b %d, %Y')[0:6])
    
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

while os.path.exists(file_path + '.bin'):
    
    file_pos = 0
    
    # ******************************DATA READ-IN***********************************
    with open(file_path + '.bin', 'rb') as fin:
        
        print('Working on ' + file_path + '.bin ...')
        starttime = time.time()
        
        word = fin.read(2)
        
        while word:
        
            # Determine the length of the buffer to be used in recognizing new buffer
            buf_ndata = struct.unpack('<H', word)[0]
            fin.seek(-2,1)
            file_pos = file_pos + (buf_ndata * 2) # Increase limits of the file pointer
        
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
                t_check = [x + t_start for x in time_stops]
            
            while fin.tell() < file_pos:
                # Read in event header data
                header = fin.read(eventheadlen * 2)
                evt_pattern, evt_timehi, evt_timelo = \
                                struct.unpack('<' + str(eventheadlen) + 'H', header)
                
                read_pattern = list(map(int, bin(evt_pattern)[-4:]))
                read_pattern.reverse()
                hit_pattern = ''.join(map(str, read_pattern))
                
                # Check if event time has past a time stop. If yes, copy spectra data and write spectra.
                evt_time = (buf_timehi * 64000 * 64000 + evt_timehi * 64000 + evt_timelo) * tunits * 1e-9
                time_lapse = evt_time - t_start
                if stop_counter < len(t_check) and evt_time > t_check[stop_counter]:
                    spec_stop_copy[stop_counter] = spec_markers.copy()
                    file_title = file_series + '_' + str(time_stops[stop_counter]) + 'sec'
                    write_spec(spe_choices, file_title, run_start_time, live_time, time_lapse, detector_choices, energy_max, spec_markers, enerfit, mca_cal, shape_cal)
                    stop_counter += 1
                
                for channel in range(4):
                    if read_pattern[channel] == 1:
                        words = fin.read(2 * 2)
                        chan_trigtime, energy = \
                                    struct.unpack('<' +str(chanheadlen) + 'H', words)
                        trigger = (evt_timehi * 64000 + chan_trigtime) * tunits
                        #data_line = map(str, [trigger, channel, hit_pattern, energy]) 
                        
                        # Store the data read in from the binary file into the ndarray,
                        # spec_markers, for eventual writing to .Spe file. A series of
                        # checks need to be met (allowed energy? channel is of interest
                        # to the user? hit pattern is of interest given the channel?).
                        # If all checks are met, increment that marker (i.e bin).
                        if energy <= energy_max:
                            if str(channel) in spe_choices.keys():
                                for i, ht in enumerate(spe_choices[str(channel)]):
                                    if hit_pattern in ht:
                                        if str(channel) in gate_choices.keys() and gate_choices[str(channel)][i] != []:
                                            for g,sig in enumerate(gate_choices[str(channel)][i]):
                                                if sig[0] <= energy <= sig[1]:
                                                    spec_markers[channel][i][energy] += 1
                                        else:
                                            spec_markers[channel][i][energy] += 1
        
            # Read word, buf_ndata, to continue loop or break
            word = fin.read(2)    
        
            buffer_no += 1
            
        elapstime = time.time() - starttime
        print('Time: %f' % elapstime)
            
    file_counter += 1
    file_path = (file_series + '%04d') % file_counter
  
# WRITE OUT CURRENT VALUES FOR 'spec_markers'
write_spec(spe_choices, file_series, run_start_time, live_time, total_time, detector_choices, energy_max, spec_markers, enerfit, mca_cal, shape_cal)

t_final = (buf_timehi * 64000 * 64000 + evt_timehi * 64000 + evt_timelo) * tunits * 1e-9 # in seconds
total_elapsed_time = t_final - t_start
print('Total elapsed: %f', total_elapsed_time)