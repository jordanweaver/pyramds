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
import datetime
import time

import numpy as np
import tables as tb

from traits.api import HasTraits, Instance, File, Str, Dict, Int
from traitsui.api import View, Item


class PyramdsModel(HasTraits):

    # Filecounter that tracks progress through file series
    file_counter = Int(1)

    def get_info(self):

        # Chop off run number and extension
        base_filepath = bin_filename[:-8]

        with open(self.base_filepath + '0001.ifm','rU') as finfo:
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

            times = {
                'start' :   run_start_time,
                'total' :   total_time,
                'live'  :   live_time
            }

            bufheadlen = int(info_str_list[33].split()[1])
            eventheadlen = int(info_str_list[34].split()[1])
            # Due to bug in PIXIE IGOR Software, need to declare head length
            # chanheadlen = #int(info_str_list[35].split()[1])
            chanheadlen = 2

        return self.times, self.bufheadlen, self.eventheadlen, self.chanheadlen





class PyramdsView(HasTraits):
    pyramds = Instance(PyramdsModel)

    bin_filename = File
    
    # Add button to submit filename and get info
    pyramds.get_info(bin_filename)

    view = View(
                Item('pyramds', style='custom', show_label=False, ),
            )

if __name__ == '__main__':
    MainWindow = PyramdsView(pyramds=PyramdsModel())
    MainWindow.configure_traits()

