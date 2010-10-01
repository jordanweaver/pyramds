#!/usr/bin/env python

# path to data series (leave off run number and extension)
file_series = 'Eu152_cal/Eu152_12hr'

detector_choices = ['1', '2'] # channels requiring spectra
spec_num = [1, 1] # number of spectra for above channels
max_spec_num = max(spec_num)

pattern_choices

for i,dc in enumerate(detector_choices):
    pattern_choices = []
    for j in range(spec_num[i]):
        patt_add = pattern_choices
            #raw_input('Enter desired hit patterns (separated by spaces) ' + \
                      #'for Detector-' + dc + ' Spectra-' + str(j+1) + ': ').split()
        pattern_choices.append(patt_add)
    spe_choices[dc] = pattern_choices

# 'spe_choices' is a dictionary containing all the choices the user defines. The
# form is to be that of keys consisting of the elements stored in 'detector_choices'
# and the values are nested lists. The lists first level is that of the different
# spectra that the user requires. Within each element of that first level is a
# list of each hit pattern that will be included in constructing that spectra.
#
# Visually: spe_choices = {'1': [['0100']], '2': [['0010'], ['1110', '1010']]}
spe_choices = {'1':[['1110', '0110', '0100', '1100']], '2':[['1110', '0110', '1010', '0010']]}

# Edit these values to correspond to current system being used (see ORTEC pdf).
enerfit = {'0':'1.568643 0.382061', '1':'1.568643 0.382061', '2':'2.745898 0.417368', '3':''}
fwhmfit = {'0':'', '1':'1.8468E+000 7E-005 3E-008', '2':'1.8468E+000 7E-005 3E-008', '3':''}
mca_cal = {'0':'3\r\n1.568643E+000 3.820606E-001 4.662797E-007', '1':'3\r\n1.568643E+000 3.820606E-001 4.662797E-007', '2':'3\r\n2.745898E+000 4.173675E-001 1.642237E-006', '3':''}
shape_cal = {'0':'', '1':'3\r\n1.8468E+000 7E-005 3E-008', '2':'', '3':''}

en_coeff = {}
fwhm_coeff = {}
for c,channel in enumerate(detector_choices):
    en_coeff[channel] = map(float, enerfit[channel].split())
    fwhm_coeff[channel] = map(float, fwhmfit[channel].split())