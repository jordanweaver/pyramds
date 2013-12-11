"""
PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)

Data parser for PIXIE List Mode binary data *.bin
The data contained in the .bin file is reformatted into an HDF5 file
that stores event information is a series of related table entries for quick
extraction of the necessary events used in spectra construction.

"""
import os
from datetime import datetime
from fnmatch import fnmatch
from os.path import dirname, join

import numpy as np
from tables import Float32Col, Int32Col, IsDescription, openFile
from traits.api import (Button, Dict, File, HasTraits, Instance, Int, List,
                        Property)
from traitsui.api import (FileEditor, Group, HSplit, Item, ListStrEditor,
                          VGroup, View)

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
    energy    = Int32Col(pos=0)
    timestamp = Float32Col(pos=1)

class AggEvent2(IsDescription):
    energy_1  = Int32Col(pos=0)
    energy_2  = Int32Col(pos=1)
    timestamp = Float32Col(pos=2)

class PyramdsParser(HasTraits):

    # Path to data file selected in UI
    selected_data_file = File()

    # Path to current working directory
    data_cwd = Property

    # File series base name (series number and extension removed)
    series_basename = Property

    # Filecounter that tracks progress through file series
    file_counter = Int(1)

    # File path currently being processed (no extension)
    active_file_path = Property

    # Only initialize buffer counter before the entire run
    buffer_no = 0

    def get_file_series(self, ext):

        file_series = []

        # Populate file_series list
        file_wildcard = self.series_basename + '*.' + ext

        for data_file in os.listdir(self.data_cwd):
            if fnmatch(join(self.data_cwd, data_file), file_wildcard):
                file_series.append(data_file)
        return file_series

    def get_bin_info(self):

        ifm_file_series = self.get_file_series('ifm')
        ifm_file_series.sort()

        total_time = 0.0
        live_time = np.zeros(4)

        for ifm_file in ifm_file_series:
            with open(ifm_file, 'rU') as f:

                info_str_list = f.readlines()

                # On first pass/file, store the series run start time
                if ifm_file is ifm_file_series[0]:
                    date_str = info_str_list[1][23:-2]
                    date_format = '%I:%M:%S %p %a, %b %d, %Y'
                    run_start_time = datetime.strptime(date_str, date_format)

                # Total time = real time; increment value for each file read
                total_time_str = info_str_list[6].split()[3]
                total_time += float(total_time_str)

                # Live time for each detector channel
                live_time_str = [info_str_list[9 + channel].split()[2]
                                 for channel in range(4)]
                live_time += np.array(map(float, live_time_str))

                self.bufheadlen = int(info_str_list[33].split()[1])
                self.eventheadlen = int(info_str_list[34].split()[1])

                # Bug in PIXIE IGOR Software, explicity set chanheadlen
                self.chanheadlen = 2  # int(info_str_list[35].split()[1])

        self.stats = {'start': run_start_time,
                      'total': str(total_time),
                      'live': map(str, live_time)}

        return self.stats

    def create_h5(self):
        """
        Open a file pointer to the HDF5 file and set up each group that will
        make up the framework of the overall file structure.
        """

        h5_filename = self.series_basename + '.h5'
        h5_title = 'Data - ' + self.series_basename

        self.h5file = openFile(h5_filename, mode='w', title=h5_title)

        self.h5_group = self.h5file.createGroup(
            self.h5file.root,
            'bin_data_parse', 'PIXIE Binary Parse')

        self.h5file.createGroup(
            self.h5file.root, 'stats', 'File time information')

        self.h5file.createGroup(
            self.h5file.root, "spectra", "Spectra Arrays")

        self.h5_gNormal = self.h5file.createGroup(
            self.h5file.root.spectra, "normal", "Normal Data")

        self.h5_gCompton = self.h5file.createGroup(
            self.h5file.root.spectra,
            "compton", "Compton-Supp Data")

        self.h5_gGGcoinc = self.h5file.createGroup(
            self.h5file.root.spectra,
            "ggcoinc", "Gamma-Gamma Data")

    def record_time_stats(self):

        self.h5file.createArray(
            self.h5file.root.stats, 'live',
            self.stats['live'], "Live stats of run")

        startt = self.stats['start'].timetuple()

        self.h5file.createArray(
            self.h5file.root.stats, 'start', [x for x in startt][:-3],
            "Start time list of run")

        self.h5file.createArray(
            self.h5file.root.stats, 'total', [self.stats['total']])

    def _get_data_cwd(self):
        return dirname(parser.series_basename)

    def _get_series_basename(self):
        return self.selected_data_file[:-8]

    def _get_active_file_path(self):
        active_file_path = \
            "{self.series_basename}{self.file_counter:0>4}".format(self=self)
        return active_file_path

class PyramdsView(HasTraits):
    Parser = Instance(PyramdsParser)

    bin_file_editor = FileEditor(filter=['*.bin'])
    hdf_file_editor = FileEditor(filter=['*.h5'])
    series_editor = ListStrEditor(editable=False)

    bin_filename = File()
    hdf_filename = File()

    bin_file_series = List()
    ifm_file_series = List()

    stats = Dict()

    parse_button = Button(label="Parse Series")

    traits_view = View(
        Group(
            VGroup(Item('bin_filename', editor=bin_file_editor, label='BIN File'),
                   HSplit(Item('bin_file_series', editor=series_editor, label='Series',width=0.4),
                          Item('bin_file_series', editor=series_editor, label='Stats'),
                          springy=True),
            show_border=True),
            VGroup(Item('parse_button', show_label=False)),
            label="PIXE PARSER"),
        Group(
            VGroup(Item('hdf_filename', editor=bin_file_editor, label='HDF File'),
            show_border=True,),label="SPECTRUM EXPORTER"),
        resizable=True,
        title="PYRAMDS",
        )

    # On filename change, update parser model and pull new .ifm stats
    def _bin_filename_changed(self, new):

        self.Parser.selected_data_file = new
        self.bin_file_series = self.Parser.get_file_series('bin')

        self.ifm_file_series = self.Parser.get_file_series('ifm')
        self.stats = self.Parser.get_bin_info()

if __name__ == '__main__':
    parser = PyramdsParser()

    pv = PyramdsView(Parser=parser)
    pv.configure_traits()
