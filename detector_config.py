# PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)
#
# Author: Jordan Weaver

# Detector Configuration File
#
# Edit these values to correspond to current system being used (see ORTEC pdf).
# The dictionary keys 0, 1, 2 correspond to the Pixie Input Channel

enerfit = {
    '0': '1.568643 0.382061',
    '1': '0.118379 0.506684',
    '2': '1.024204 0.413356',
    '3': ''}

fwhmfit = {
    '0': '1.8468E+000 7E-005 3E-008',
    '1': '1.828179E+000 5.470951E-005 3.872378E-008',
    '2': '1.703389E+000 -3.886789E-005 6.447196E-008',
    '3': ''}

mca_cal = {
    '0': '3\r\n1.568643E+000 3.820606E-001 4.662797E-007',
    '1': '2\r\n1.18379E-001 5.06684E-001',
    '2': '2\r\n1.024204E+000 4.13356E-001',
    '3': ''}

shape_cal = {
    '0': '3\r\n1.8468E+000 7E-005 3E-008',
    '1': '3\r\n1.828179E+000 5.470951E-005 3.872378E-008',
    '2': '3\r\n1.703389E+000 -3.886789E-005 6.447196E-008',
    '3': ''}
