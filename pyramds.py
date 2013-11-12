"""
PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)

Data parser for PIXIE List Mode binary data *.bin
The data contained in the .bin file is reformatted into an HDF5 file
that stores event information is a series of related table entries for quick
extraction of the necessary events used in spectra construction.

"""

import os
import sys
import struct
import time

import numpy as np
import tables as tb

from traits.api import HasTraits, Instance
from enaml.widgets.api import Window, RadioButton


class PyramdsView(HasTraits):

    # Filename for .bin files
    bin_filename = File

    # PLACEHOLDER FOR PARSING OPTIONS

    # Save information
    save_direcotry = Directory(".")
    parse_run_button = Button(label="Parse Data")

