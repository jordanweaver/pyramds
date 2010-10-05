#! /usr/bin/env python
import os
import time
import subprocess

import numpy 
import tables
import scipy

# ETS imports
from enthought.traits.api import HasTraits, Instance, Str, Enum, Any, List, Array, Dict, Button, NO_COMPARE, Float, \
    File, Directory, on_trait_change, HTML, Range

from enthought.traits.ui.api import View, Item, Group, VGroup, HGroup, EnumEditor, TableEditor, TabularEditor, spring, Spring, \
    HTMLEditor, HSplit, VSplit, VGrid, RangeEditor
from enthought.traits.ui.ui_editors.array_view_editor import ArrayViewEditor
from enthought.traits.ui.menu import Action

from enthought.chaco.api import Plot, ArrayPlotData
from enthought.enable.api import Component, ComponentEditor, Window

from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom

from spectrum_plot_view import make_spectrum_plot, save_plot
from detection_limits_helpers import detection_limits_to_html

np = numpy
tb = tables

size = (500, 500)

class PyramdsView(HasTraits):
    plot = Instance(Component)
    index_selections = List
    peak_indices = Dict

    # Database info
    filename = File
    datafile = Any()
    dfr      = Any()

    # Save Information
    save_directory = Directory(".")
    save_figure = Button(label="Save Figure")


    # plot data
    chan = Array
    hist = Array
    pchn = Array
    peak = Array

    title = Str("")

    detector = Enum(["1", "2"])
    spectrum = Enum(["Normal", "Compton", "Gamma-Gamma"])
    spectrum_set_names = {
        "Normal":      "norm",
        "Compton":     "compt",
        "Gamma-Gamma": "gg",
        }
    spectrum_group_names = {
        "Normal":      "normal",
        "Compton":     "compton",
        "Gamma-Gamma": "ggcoinc",
        }

    # UI Stuff
    start_time_low = Float(0.0)
    start_time_high = Float(1.0)
    end_time_low = Float(0.0)
    end_time_high = Float(1.0)

    start_time = Range(low=0.0, high=np.inf, value=9999.99)
    end_time = Range(low=0.0, high=np.inf, value=99999.99)

    # Detection Limits 
    isotope_enum = List
    peaknum_enum = List

    isotope = Str
    peaknum = Str

    detection_limits_html = HTML
    save_detection_limits = Button(label="Save Detection Limits")

    detection_limits = Dict({
        'lookup_selected': {},
        'clicked_selected': {},
        })

    detection_limits_title = Str("Detection Limits")

    traits_view = View(
        VGroup(
            Group(Item('filename', label="Data File")),
            HSplit(
                Item('plot', editor=ComponentEditor(size=size), show_label=False),
                VGroup(
                    HGroup(spring, Item('detection_limits_title', style='readonly', show_label=False), spring),
                    VGroup(
                        Item('isotope', editor=EnumEditor(name='isotope_enum'), label="Isotope"),
                        Item('peaknum', editor=EnumEditor(name='peaknum_enum'), label="Peak No.")
                    ),
                    Item('detection_limits_html', editor=HTMLEditor(), show_label=False),
                    Item('save_detection_limits', show_label=False),
                    label="Detection Limits",
                ),
            ),
            HGroup(
                Item('detector'), 
                spring,
                Item('spectrum'), 
            ),
            HGroup(
                Item('save_directory', width=400, label="Save Dir"),
                spring,
                Item('save_figure', show_label=False),
            ),
            VGroup(
                Item('start_time', editor=RangeEditor(low=9999.99, high=99999.99, low_name='start_time_low', high_name='start_time_high', format='%.2F', mode='slider'), 
                    format_str='%.2F',  
                    label="Start Time"),
                Item('end_time', editor=RangeEditor(low=9999.99, high=99999.99, low_name='end_time_low', high_name='end_time_high', format='%.2F', mode='slider'), 
                    format_str='%.2F', 
                    label="End Time"), 
            ),
        ),
        width=500, 
        height=500, 
        resizable=True, 
        title="Pyramds Visualizer",
        )

    def get_histogram_data_set(self):
        hist_grp = getattr(self.dfr.spectra, self.spectrum_group_names[self.spectrum])
        hist_set = getattr(hist_grp, "{0}{1}_spec".format(self.spectrum_set_names[self.spectrum], self.detector))
        return hist_set

    def load_histogram_data(self):
        hist_set = self.get_histogram_data_set()
        hist = hist_set[-1]
        chan = np.arange(len(hist))
        pchn = np.array([])
        peak = np.array([])

        self.peak_indices = {}

        self.chan = chan
        self.hist = hist
        self.pchn = pchn
        self.peak = peak

    def get_total_time(self):
        if self.dfr == None:
            dt = 99999.99
        else:
            dt = self.dfr.bin_data_parse.readout[-1]['timestamp'] - self.dfr.bin_data_parse.readout[0]['timestamp']
        return float(dt)

    def draw_plot(self):
        label = "Detector {0} {1} Spectrum".format(self.detector, self.spectrum)
        plot = make_spectrum_plot(self.chan, self.hist, self.pchn, self.peak, label)
        self.plot = plot

    def draw_peak_plot(self):
        self.plot.plots['plot1'][0].index.set_data(self.pchn)
        self.plot.plots['plot1'][0].value.set_data(self.peak)

    def draw_detection_limits(self):
        detection_limits_html = detection_limits_to_html(self.detection_limits)
        self.detection_limits_html = detection_limits_html

    #
    # Set Trait Defaults
    #

    def _filename_default(self):
        return ""

    def _datafile_default(self):
        if os.path.exists(self.filename):
            return tb.openFile(self.filename, "r")
        else:
            return None

    def _dfr_default(self):
        if self.datafile == None:
            return None
        else:
            return self.datafile.root

    def _chan_default(self):
        if self.dfr == None:
            return np.array([0])
        else:
            return np.arange(len(self.hist))

    def _hist_default(self):
        if self.dfr == None:
            return np.array([1000])
        else:
            hist_set = self.get_histogram_data_set()
            hist = hist_set[-1]
            return hist

    def _pchn_default(self):
        return np.array([])

    def _peak_default(self):
        return np.array([])

    def _plot_default(self):
        plot = make_spectrum_plot(self.chan, self.hist, self.pchn, self.peak)
        return plot

    def _start_time_low_default(self):
         return 9999.99

    def _start_time_high_default(self):
         return self.get_total_time()

    def _end_time_low_default(self):
         return 9999.99

    def _end_time_high_default(self):
         return self.get_total_time()

    def _detector_default(self):
        return "1"

    def _spectrum_default(self):
        return "Normal"

    def _isotope_enum_default(self):
        return ["None"]

    def _peaknum_enum_default(self):
        return ["0"]

    def _isotope_default(self):
        return "None"

    def _peaknum_default(self):
        return "0"

    def _detection_limits_html_default(self):
        return "<b>No peaks selected.</b>"

    #
    # Define Traits Changed
    #

    def _filename_changed(self, new):
        # Close old hdf5 file
        if self.datafile != None:
            self.datafile.close()

        # Load new hdf5 file
        self.datafile = tb.openFile(self.filename, "r")
        self.dfr = self.datafile.root

        self.load_histogram_data()

        # reset times
        self.start_time_low = 0.0
        self.end_time_low = 0.0
        self.end_time_high = self.get_total_time()
        self.start_time_high = self.get_total_time()
        self.start_time = 0.0
        self.end_time = self.get_total_time()

        # Redraw plot
        self.draw_plot()

    def _start_time_changed(self, old, new):
        self.end_time_low = new

    def _end_time_changed(self, old, new):
        self.start_time_high = new

    def _detector_changed(self):
        self.load_histogram_data()
        self.draw_plot()

    def _spectrum_changed(self):
        self.load_histogram_data()
        self.draw_plot()

    def _plot_changed(self):
        self.plot.plots['plot0'][0].index.on_trait_change(self._index_metadata_handler, 'metadata_changed')
        self.index_selections = self.plot.plots['plot0'][0].index.metadata['selections']

    def _index_metadata_handler(self):
        self.index_selections = self.plot.plots['plot0'][0].index.metadata['selections']

    def _index_selections_changed(self, old, new):
        # Some fake way of selecting the peak
        pi = []
        for si in self.index_selections:
            # Prevent the peak-finding algorithm from running 
            # unnecessarily often by caching the results
            if si not in self.peak_indices:
                self.peak_indices[si] = range(si-10, si+11) 
            pi.extend( self.peak_indices[si] )

        # Remove old selections
        for key in self.peak_indices.keys():
            if key not in self.index_selections:
                del self.peak_indices[key]

        # Set the points to highlight
        pchn = np.array(pi)
        if pi == []:
            peak = np.array([])
        else:
            peak = self.hist[pchn]

        # Set the traits
        self.pchn = pchn
        self.peak = peak

        # Redraw the plot
        self.draw_peak_plot()

        # Calculate detection limits
        for si in self.index_selections:
            if si not in self.detection_limits['clicked_selected']:
                self.detection_limits['clicked_selected'][si] = 0.5

        # Remove old detection limits
        for key in self.detection_limits['clicked_selected'].keys():
            if key not in self.index_selections:
                del self.detection_limits['clicked_selected'][key]

        # Redraw the detection limits display
        self.draw_detection_limits()
        
    #
    # Buttons Fired methods
    #

    def _save_figure_fired(self):
        w = 800
        h = 600
        savefile = "{0}/{1}{2}_Spectrum.png".format(self.save_directory, self.spectrum, self.detector)
        save_plot(self.plot, savefile, w, h)

if __name__ == "__main__":
    pv = PyramdsView()
    pv.configure_traits()

    # close the hdf5 file that is still open
    if pv.datafile != None:
        pv.datafile.close()
