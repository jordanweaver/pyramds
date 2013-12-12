"""
PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)

Data parser for PIXIE List Mode binary data *.bin
The data contained in the .bin file is reformatted into an HDF5 file
that stores event information is a series of related table entries for quick
extraction of the necessary events used in spectra construction.

"""
from traits.api import HasTraits, Instance, File, List, Dict, Button
from traitsui.api import (FileEditor, Group, HSplit, Item,
                          ListStrEditor, TextEditor, VGroup, View)

from pyramds_model import PyramdsParser

class SeriesView(HasTraits):
    bin_file_series = List(label="BIN File", desc="Select file from series")
    series_editor = ListStrEditor(editable=False)

    view = View(Item('bin_file_series',
                     editor=series_editor, style='readonly'))

class PyramdsView(HasTraits):
    parser = Instance(PyramdsParser)
    series_view = Instance(SeriesView)

    bin_file_editor = FileEditor(filter=['*.bin'])
    hdf_file_editor = FileEditor(filter=['*.h5'])
    stats_editor = TextEditor(multi_line=True)

    bin_filename = File()
    hdf_filename = File()

    ifm_file_series = List()

    stats = Dict()

    parse_button = Button(label="Parse Series")

    traits_view = View(
        Group(
            VGroup(Item('bin_filename', editor=bin_file_editor, label='BIN File'),
                   HSplit(Item('series_view', style='custom'),
                          Item('stats', editor=stats_editor, label='Stats'),
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

        self.parser.selected_data_file = new
        self.bin_file_series = self.parser.get_file_series('bin')

        self.stats = self.parser.get_bin_info()

if __name__ == '__main__':
    pp = PyramdsParser()
    sv = SeriesView()

    pv = PyramdsView(parser=pp, series_view=)
    pv.configure_traits()
