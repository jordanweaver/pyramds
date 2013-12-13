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
    bin_file_series = List(label="Files", desc="Select file from series")
    series_editor = ListStrEditor(editable=False)

    view = View(Item('bin_file_series',
                     editor=series_editor, style='readonly'))

class StatsView(HasTraits):
    ifm_file_series = List()
    stats_editor = TextEditor(multi_line=True)
    stats = Dict()

    view = View(Item('stats',
                     editor=stats_editor, style='readonly'))

    # ADD METHODS FOR DISPLAYING THE FULL TEXT OF THE MESSAGE

class PyramdsView(HasTraits):
    parser = Instance(PyramdsParser)
    series_view = Instance(SeriesView)
    stats_view = Instance(StatsView)

    bin_file_editor = FileEditor(filter=['*.bin'])
    hdf_file_editor = FileEditor(filter=['*.h5'])

    bin_filename = File()
    hdf_filename = File()

    parse_button = Button(label="Parse Series")

    traits_view = View(
        Group(
            VGroup(Item('bin_filename', editor=bin_file_editor,
                        label='BIN File'),
                   HSplit(Item('series_view', style='custom',
                               show_label=False),
                          Item('stats', style='custom', show_label=False),
                          springy=True),
                   show_border=True),
            VGroup(Item('parse_button', show_label=False)),
            label="PIXE PARSER"),
        Group(
            VGroup(Item('hdf_filename', editor=bin_file_editor,
                        label='HDF File'),
                   show_border=True,),
            label="SPECTRUM EXPORTER"),
        resizable=True,
        title="PYRAMDS")

    # On filename change, update parser model and pull new .ifm stats
    def _bin_filename_changed(self, new):

        self.parser.selected_data_file = new
        self.series_view.bin_file_series = self.parser.get_file_series('bin')

        self.stats_view.stats = self.parser.get_bin_info()

if __name__ == '__main__':
    pp = PyramdsParser()
    serv = SeriesView()
    stav = StatsView()

    pv = PyramdsView(parser=pp, series_view=serv, stats_view=stav)
    pv.configure_traits()
