# PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)
#
# Data parser for PIXIE List Mode binary data *.bin
# The data contained in the .bin file is reformatted into an HDF5 file
# that stores event information is a series of related table entries for quick
# extraction of the necessary events used in spectra construction.
#
# Author: Jordan Weaver

# Standard Library Imports
import textwrap

# External Imports
from traits.api import Button, Dict, File, HasTraits, Instance, List, Str
from traitsui.api import (FileEditor, Group, HSplit, Item, ListStrEditor,
                          TextEditor, VGroup, View)

# Internal Imports
from parser_model import PyramdsParser, SpectrumExporter


class SeriesView(HasTraits):
    bin_file_series = List(label="Files", desc="Select file from series")
    series_editor = ListStrEditor(editable=False)

    view = View(Group(Item('bin_file_series',
                           editor=series_editor, style='readonly',
                           show_label=False),
                show_border=True, label="DATA SERIES"))

class StatsView(HasTraits):
    ifm_file_series = List()
    stats_editor = TextEditor()
    stats = Dict()

    disp_str = Str()

    view = View(Group(Item('disp_str',
                           editor=stats_editor, style='readonly',
                           show_label=False),
                show_border=True, label="SERIES INFORMATION"))

    def _disp_str_default(self):
        disp_str = "No series currently selected"
        return disp_str

    def _stats_changed(self):
        disp_str = textwrap.dedent("""\

            Series was started: {start}

            Total run time: {total}

            Live time for each detector:
            DET0         DET1        DET2
            {live[0]}    {live[1]}   {live[2]}
            """)

        self.disp_str = disp_str.format(**self.stats)

class ExporterView(HasTraits):
    pass

class PyramdsView(HasTraits):
    parser = Instance(PyramdsParser, ())
    exporter = Instance(SpectrumExporter, ())
    series_view = Instance(SeriesView, ())
    stats_view = Instance(StatsView, ())
    exporter_view = Instance(ExporterView, ())

    bin_file_editor = FileEditor(filter=['*.bin'])
    hdf_file_editor = FileEditor(filter=['*.h5'])

    bin_filename = File()
    hdf_filename = File()

    parse_button = Button(label="Parse Series")
    export_button = Button(label="Export Spectra")

    traits_view = View(
        Group(
            VGroup(Item('bin_filename', editor=bin_file_editor,
                        label='BIN File'),
                   HSplit(Item('series_view', style='custom',
                               show_label=False),
                          Item('stats_view', style='custom', show_label=False),
                          springy=True),
                   show_border=True),
            VGroup(Item('parse_button', show_label=False)),
            label="PIXE PARSER"),
        Group(
            VGroup(Item('hdf_filename', editor=hdf_file_editor,
                        label='HDF File'),
                   Item('exporter_view', style='custom', show_label=False),
                   show_border=True,),
            VGroup(Item('export_button', show_label=False)),
            label="SPECTRUM EXPORTER"),
        resizable=True,
        title="PYRAMDS")

    # On bin filename change, update parser model and pull new .ifm stats
    def _bin_filename_changed(self, new):

        self.parser.data_file = new
        self.series_view.bin_file_series = self.parser.get_file_series('bin')

        self.stats_view.stats = self.parser.get_bin_info()

    # On hdf filename change, update exporter model and pull/display info
    def _hdf_filename_changed(self, new):

        self.exporter.data_file = self.hdf_filename

    def _parse_button_fired(self):
        # Check if old HDF5 file is still around
        if self.parser.h5file is not None:
            self.parser.h5file.close()

        # Open new HDF5 file, parse data, store spectra structurs, and close
        self.parser.start_parse()
        self.parser.store_spectra_h5()

        self.hdf_filename = self.parser.h5_filename

    def _export_button_fired(self):
        # Check if old HDF5 file is still around
        if self.parser.h5file is not None:
            self.parser.h5file.close()

        self.exporter.write_spec()

if __name__ == '__main__':

    pv = PyramdsView()
    pv.configure_traits()
