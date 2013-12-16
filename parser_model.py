# PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)
#
# Author: Jordan Weaver

import struct
import os
import numpy as np

from tables import Float32Col, Int32Col, IsDescription

# Internal Imports
from parser_setup import PyramdsBase

# Setup PyTables metaclasses for use in Table constructor
class GammaEvent(IsDescription):

    # Energy reading from Pixie detector channels 0, 1, 2
    energy_0 = Int32Col(pos=0)
    energy_1 = Int32Col(pos=1)
    energy_2 = Int32Col(pos=2)

    # Time differences between detector events
    # (e.g. "deltaT_01" is between energy_0 and energy_1)
    deltaT_01 = Float32Col(pos=3)
    deltaT_02 = Float32Col(pos=4)
    deltaT_12 = Float32Col(pos=5)

    timestamp = Float32Col(pos=6)

class AggEvent1(IsDescription):
    energy = Int32Col(pos=0)
    timestamp = Float32Col(pos=1)

class AggEvent2(IsDescription):
    energy_1 = Int32Col(pos=0)
    energy_2 = Int32Col(pos=1)
    timestamp = Float32Col(pos=2)

class PyramdsParser(PyramdsBase):

    def start_parse(self):

        self.create_h5()
        self.record_time_stats()

        # This table is where the data will be placed after unpacking it
        # from binary
        self.table = self.h5file.createTable(self.h5_group, 'readout',
                                             GammaEvent, "Data readout")

        # Only start the buffer count before the entire run, not each file
        buffer_no = 0

        for data_file in self.get_file_series('bin'):

            data_path = os.path.join(self.data_cwd, data_file)

            file_pos = 0

            with open(data_path, 'rb') as fin:
                print('Working on ' + data_file)

                # Pointer object to place values on in the row for
                # each event. Must create new instance each time .flush
                # is called (every .bin file)
                event = self.table.row

                word = fin.read(2)

                while word:

                    # Determine the length of the buffer to be used in
                    # recognizing new buffer
                    buf_ndata = struct.unpack('<H', word)[0]
                    fin.seek(-2, 1)

                    # Increase limits of the file pointer
                    file_pos += (buf_ndata * 2)

                    # Read in buffer header data
                    head = fin.read(self.bufheadlen * 2)
                    head_fmt = '<' + str(self.bufheadlen) + 'H'

                    (buf_ndata, buf_modnum,
                     buf_format, buf_timehi,
                     buf_timemi, buf_timelo) = struct.unpack(head_fmt, head)

                    # Remember the time of the first buffer of entire run.
                    # Use this for comparing time stops.
                    if buffer_no == 0:
                        t_start_hi = buf_timehi * 64000 * 64000
                        t_start_mi = buf_timemi * 64000
                        t_start_lo = buf_timelo
                        self.t_start = \
                            (t_start_hi + t_start_mi + t_start_lo) * \
                            self.tunits * 1e-9  # in seconds

                    while fin.tell() < file_pos:
                        # Read in event header data
                        head = fin.read(self.eventheadlen * 2)
                        head_fmt = '<' + str(self.eventheadlen) + 'H'

                        (evt_pattern,
                         evt_timehi,
                         evt_timelo) = struct.unpack(head_fmt, head)

                        read_pattern = list(map(int, bin(evt_pattern)[-4:]))
                        read_pattern.reverse()

                        trigger_vals = [float('nan')] * 3
                        for chan in range(3):
                            if read_pattern[chan] == 1:
                                words = fin.read(2 * 2)

                                head_fmt = '<' + str(self.chanheadlen) + 'H'
                                chan_trigtime, energy = struct.unpack(head_fmt,
                                                                      words)

                                trigger_vals[chan] = \
                                    float((evt_timehi * 64000 +
                                          chan_trigtime) * self.tunits)

                                # Store the data read in from the binary file
                                # into the HDF5 table.
                                if ((energy > self.energy_max) & (chan > 0)):
                                    event['energy_' + str(chan)] = -1
                                else:
                                    event['energy_' + str(chan)] = energy
                            elif read_pattern[chan] == 0:
                                event['energy_' + str(chan)] = -1

                        #ftest.write('\n')

                        event['deltaT_01'] = abs(trigger_vals[0] -
                                                 trigger_vals[1])
                        event['deltaT_02'] = abs(trigger_vals[0] -
                                                 trigger_vals[2])
                        event['deltaT_12'] = abs(trigger_vals[1] -
                                                 trigger_vals[2])

                        event['timestamp'] = (buf_timehi * 64000 * 64000 +
                                              evt_timehi * 64000 +
                                              evt_timelo) * self.tunits * 1e-9

                        event.append()
                        self.table.flush()

                    # Read word, buf_ndata, to continue loop or break
                    word = fin.read(2)

                    buffer_no += 1
                    if buffer_no % 100 == 0:
                        print('Buffer No. %d' % buffer_no)

                    # Flush data to the HFD5 table and start new buffer
                    self.table.flush()

        # in seconds
        self.t_final = (buf_timehi * 64000 * 64000 +
                        evt_timehi * 64000 + evt_timelo) * self.tunits * 1e-9
        self.t_duration = self.t_final - self.t_start
        self.t_array_dim = int(np.ceil(self.t_duration / self.t_steps))

    def store_spectra_h5(self):

        ######################################################################
        # Store Det-1 & 2 arrays containing Normal counts ####################
        ######################################################################

        print('Started creating normal data aggregates...')

        norm12table = self.h5file.createTable(self.h5_gNormal, 'norm_evts12',
                                              AggEvent2,
                                              "Normal Log Data")
        norm_event = norm12table.row

        for row in self.table.where("""(energy_1 != -1) | (energy_2 != -1)"""):
            norm_event['energy_1'] = row['energy_1']
            norm_event['energy_2'] = row['energy_2']
            norm_event['timestamp'] = row['timestamp']
            norm_event.append()
            norm12table.flush()

        # Det 1
        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)
        ti = 0
        for row in norm12table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy_1'] >= 0:
                dt_temp[row['energy_1']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gNormal, 'norm1_spec',
                                dt_array,
                                "Normal Time-Chunked Spec Array - Det 1")

        # Det 2
        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)
        ti = 0
        for row in norm12table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy_2'] >= 0:
                dt_temp[row['energy_2']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gNormal, 'norm2_spec',
                                dt_array,
                                "Normal Time-Chunked Spec Array - Det 2")

        ######################################################################
        # Store Det-1 & 2 arrays containing Compton counts ###################
        ######################################################################

        print('Started creating anti-coincidence data aggregates...')

        compt1table = self.h5file.createTable(self.h5_gCompton, 'compt_evts1',
                                              AggEvent1,
                                              "Compton Log Data - Det 1")
        compt_event = compt1table.row

        query = "(energy_0 == -1) & (energy_1 != -1) & (energy_2 == -1)"
        for row in self.table.where(query):
            compt_event['energy'] = row['energy_1']
            compt_event['timestamp'] = row['timestamp']
            compt_event.append()
            compt1table.flush()

        compt2table = self.h5file.createTable(self.h5_gCompton, 'compt_evts2',
                                              AggEvent1,
                                              "Compton Log Data - Det 2")
        compt_event = compt2table.row

        query = "(energy_0 == -1) & (energy_1 == -1) & (energy_2 != -1)"
        for row in self.table.where(query):
            compt_event['energy'] = row['energy_2']
            compt_event['timestamp'] = row['timestamp']
            compt_event.append()
            compt2table.flush()

        # Det 1
        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)

        ti = 0
        for row in compt1table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy'] >= 0:
                dt_temp[row['energy']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gCompton, 'compt1_spec', dt_array,
                                "Compton-Supp Time-Chunked Spec Array - Det 1")

        # Det 2
        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)

        ti = 0
        for row in compt2table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy'] >= 0:
                dt_temp[row['energy']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gCompton, 'compt2_spec', dt_array,
                                "Compton-Supp Time-Chunked Spec Array - Det 2")

        ######################################################################
        # Store Det-1 & 2 arrays containing Gamma-Gamma counts ###############
        ######################################################################

        print('Started creating coincidence data aggregates...')

        gg1table = self.h5file.createTable(self.h5_gGGcoinc, 'gg_evts1',
                                           AggEvent1, "G-G Log Data - Det 1")
        gg_event = gg1table.row

        # Det 1
        teststr1 = ("((energy_0 == -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_12 < {0}))")

        teststr2 = ("((energy_0 != -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 == -1) & "
                    "(deltaT_01 < {0}))")

        teststr3 = ("((energy_0 != -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_01 < {0}))")

        teststr4 = ("((energy_0 != -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_12 < {0}))")

        teststr = ' | '.join([teststr1, teststr2, teststr3, teststr4])
        query = teststr.format(self.short_window)

        for row in self.table.where(query):
            gg_event['energy'] = row['energy_1']
            gg_event['timestamp'] = row['timestamp']
            gg_event.append()
            gg1table.flush()

        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)

        ti = 0
        for row in gg1table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy'] >= 0:
                dt_temp[row['energy']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gGGcoinc, 'gg1_spec', dt_array,
                                "G-G Time-Chunked Spec Array - Det 1")

        # Det 2
        gg2table = self.h5file.createTable(self.h5_gGGcoinc, 'gg_evts2',
                                           AggEvent1, "G-G Log Data - Det 2")
        gg_event = gg2table.row

        teststr1 = ("((energy_0 == -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_12 < {0}))")

        teststr2 = ("((energy_0 != -1) & "
                    "(energy_1 == -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_02 < {0}))")

        teststr3 = ("((energy_0 != -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_02 < {0}))")

        teststr4 = ("((energy_0 != -1) & "
                    "(energy_1 != -1) & "
                    "(energy_2 != -1) & "
                    "(deltaT_12 < {0}))")

        teststr = ' | '.join([teststr1, teststr2, teststr3, teststr4])
        query = teststr.format(self.short_window)

        for row in self.table.where(query):
            gg_event['energy'] = row['energy_2']
            gg_event['timestamp'] = row['timestamp']
            gg_event.append()
            gg2table.flush()

        dt_temp = np.zeros(self.energy_max + 1, dtype=np.int32)
        dt_array = np.zeros((self.t_array_dim + 1, self.energy_max + 1),
                            dtype=np.int32)
        ti = 0
        for row in gg2table:
            if (row['timestamp'] - self.t_start) >= (ti + 1) * self.t_steps:
                ti += 1
                dt_array[ti] = dt_temp.copy()
            if row['energy'] >= 0:
                dt_temp[row['energy']] += 1
        dt_array[-1] = dt_temp.copy()

        self.h5file.createArray(self.h5_gGGcoinc, 'gg2_spec', dt_array,
                                "G-G Time-Chunked Spec Array - Det 2")

class SpectrumExporter(PyramdsParser):
    pass
