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
from os.path import basename, dirname, join

from tables import Float32Col, Int32Col, IsDescription, openFile
from traits.api import Button, File, HasTraits, Instance, Int, Property, Str, List
from traitsui.api import Group, VGroup, HSplit, Item, View, ListStrEditor, FileEditor

# Setup PyTables metaclasses for use in Table constructor
class GammaEvent(IsDescription):

    # Energy reading from Pixie detector channels 0, 1, 2
    energy_0  = Int32Col(pos=0)
    energy_1  = Int32Col(pos=1)
    energy_2  = Int32Col(pos=2)

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

        for file in os.listdir(self.data_cwd):
            if fnmatch(join(self.data_cwd, file), file_wildcard):
                file_series.append(file)
        return file_series

    def get_bin_info(self):

        active_file_ifm = self.active_file_path + ".ifm"

        with open(active_file_ifm, 'rU') as finfo:

            info_str_list = finfo.readlines()

            date_str = info_str_list[1][23:-2]
            date_format = '%I:%M:%S %p %a, %b %d, %Y'
            run_start_time = datetime.strptime(date_str, date_format)

            # Total time = real time
            total_time = info_str_list[6].split()[3]

            # Live time for each detector channel
            live_time = [info_str_list[9 + channel].split()[2] for channel in range(4)]

            self.times = {
                'start' :   run_start_time,
                'total' :   total_time,
                'live'  :   live_time
            }

            self.bufheadlen = int(info_str_list[33].split()[1])
            self.eventheadlen = int(info_str_list[34].split()[1])

            # Due to bug in PIXIE IGOR Software, explicity set head length
            self.chanheadlen = 2 #int(info_str_list[35].split()[1])

    def create_h5(self):
        """
        Open a file pointer to the HDF5 file and set up each group that will
        make up the framework of the overall file structure.
        """

        h5_filename = self.series_basename + '.h5'
        h5_title = 'Data - ' + self.series_basename

        self.h5file = openFile(h5_filename, mode = 'w', title = h5_title)

        self.h5_group = self.h5file.createGroup(self.h5file.root,
            'bin_data_parse', 'PIXIE Binary Parse')

        self.h5file.createGroup(self.h5file.root,
            'times', 'File time information')

        self.h5file.createGroup(self.h5file.root,
            "spectra", "Spectra Arrays")

        self.h5_gNormal = self.h5file.createGroup(self.h5file.root.spectra,
            "normal", "Normal Data")

        self.h5_gCompton = self.h5file.createGroup(self.h5file.root.spectra,
            "compton", "Compton-Supp Data")

        self.h5_gGGcoinc = self.h5file.createGroup(self.h5file.root.spectra,
            "ggcoinc", "Gamma-Gamma Data")

    def record_time_stats(self):

        self.h5file.createArray(self.h5file.root.times, 'live', times['live'], "Live times of run")
        startt = times['start'].timetuple()
        self.h5file.createArray(self.h5file.root.times, 'start', [x for x in startt][:-3], "Start time list of run")
        self.h5file.createArray(self.h5file.root.times, 'total', [times['total']])


    def _get_data_cwd(self):
        return dirname(parser.series_basename)

    def _get_series_basename(self):
        return self.selected_data_file[:-8]

    def _get_active_file_path(self):
        return "{self.series_basename}{self.file_counter:0>4}".format(self=self)

class PyramdsView(HasTraits):
    Parser = Instance(PyramdsParser)

    bin_file_editor = FileEditor(filter=['*.bin'])
    hdf_file_editor = FileEditor(filter=['*.h5'])
    series_editor = ListStrEditor(editable=False)

    bin_filename = File()
    hdf_filename = File()

    bin_file_series = List()
    ifm_file_series = List()

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
        resizable = True,
        title = "PYRAMDS",
        )

    # On filename change, update parser model and pull new .ifm stats
    def _bin_filename_changed(self, new):

       self.Parser.selected_data_file = new
       self.bin_file_series = self.Parser.get_file_series('bin')

       self.ifm_file_series = self.Parser.get_file_series('ifm')
       for ifm in self.ifm_file_series:
           pass

if __name__ == '__main__':
    parser = PyramdsParser()

    PyramdsWindow = PyramdsView(Parser=parser)
    PyramdsWindow.configure_traits()
