#!/usr/bin/env python
import numpy as np

# path to data series (leave off run number and extension)
#file_series = '/python/PIXIE_runs/debug-' #EuCal/Eu152_92710-'
file_series = 'debug-' #EuCal/Eu152_92710-'

#lib_name = '/python/PYRAMDS/sig_library.txt' 
lib_name = 'sig_library.txt' 

spec_user_input = {
    'type'      : [],   # 0-Normal, 1-Compton, 2-GG, 3-GG/Compton, 4-Specgram
    'det'       : [],   # Declare detectors requiring spectra
    'tstart'    : 0.,   # Enter float value for detector start (in seconds)
    'tend'      : 0.,   # Enter value for period end (zero for end of file)
    'roi'       : [],   # Enter signature IDs for ROI reporting
    'e_gates'   : [],   # Enter signature IDs for gating
}

# DETECTOR CHARACTERISTICS (i.e. energy calibrations, resolution curves, etc.)

# Edit these values to correspond to current system being used (see ORTEC pdf).
enerfit = {
    '0' :   '1.568643 0.382061',
    '1' :   '1.568643 0.382061',
    '2' :   '2.745898 0.417368',
    '3' :   ''
    }

fwhmfit = {
    '0' :   '1.8468E+000 7E-005 3E-008',
    '1' :   '1.8468E+000 7E-005 3E-008',
    '2' :   '1.8468E+000 7E-005 3E-008',
    '3' :   ''
    }

mca_cal = {
    '0' :   '3\r\n1.568643E+000 3.820606E-001 4.662797E-007',
    '1' :   '3\r\n1.568643E+000 3.820606E-001 4.662797E-007',
    '2' :   '3\r\n2.745898E+000 4.173675E-001 1.642237E-006',
    '3' :   ''
    }

shape_cal = {
    '0' :   '',
    '1' :   '3\r\n1.8468E+000 7E-005 3E-008',
    '2' :   '',
    '3' :   ''
    }

en_coeff = {}
fwhm_coeff = {}
for channel in range(3):
    en_coeff[str(channel)] = map(float, enerfit[str(channel)].split())
    fwhm_coeff[str(channel)] = map(float, fwhmfit[str(channel)].split())    

# Supply the path to the custom signature library. Builds necessary lookup arrays

with open(lib_name, 'r') as ginput:
    gamma_lib = {}
    
    line_data = ginput.readline().split()
    while line_data:
        ZAID = line_data[0]
        name_str = line_data[1]
        energy = line_data[2]
        
        gamma_lib[ZAID] = [name_str, float(energy)]
        
        line_data = ginput.readline().split()
    
    sig_lookup = [{}]*4
    for c,chn in enumerate(en_coeff.keys()):
        for z,zaid in enumerate(gamma_lib.keys()):
            top_val = gamma_lib[zaid][1] - en_coeff[chn][0]
            bottom_val = en_coeff[chn][1]
            cent_chn = top_val/bottom_val
            width = fwhm_coeff[chn][0] + fwhm_coeff[chn][1] * cent_chn + fwhm_coeff[chn][2] * cent_chn * cent_chn
            
            left_marker = int(round(cent_chn - 1.75*width))
            right_marker = int(round(cent_chn + 1.75*width))
            
            sig_lookup[int(chn)][zaid] = [left_marker, right_marker]

# Detector system variables
energy_max = 8192 # maximum number of bins (energies) to be stored for detectors
tunits = 1000./75. # 13.3333... nanoseconds (for PIXIE timing)
short_window = 90. # timing window for gamma-gamma condition (nanoseconds)

# Timing chunks for storing spectrum states
t_steps = 60. # seconds

# Store the bin-energy mapping for plotting spectra
en_plot_axis = {}
for c, chn in enumerate(en_coeff.keys()):
    en_plot_axis[chn] = np.zeros( (energy_max + 1), dtype=np.float32)
    for marker in range(energy_max + 1):
        en_plot_axis[chn][marker] = en_coeff[chn][0] + en_coeff[chn][1] * marker
