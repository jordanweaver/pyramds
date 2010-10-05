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
from detection_limits_helpers import detection_limits_to_html, detection_limits_to_tsv

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

    start_time = Range(low=0.0, high=np.inf, value=0.0)
    end_time = Range(low=0.0, high=np.inf, value=1.0)


    # Detection Limits 
    detection_limits_title = Str("Detection Limits")

    isotope_enum = List
    peaknum_enum = List

    isotope = Str
    peaknum = Str

    add_peak = Button(label="Add Peak")
    del_peak = Button(label="Remove Peak")

    detection_limits = Dict({
        'lookup_selected': {},
        'clicked_selected': {},
        })


    detection_limits_html = HTML

    save_detection_limits = Button(label="Save Detection Limits")

    # Traits view
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
                    HGroup(
                        Item('add_peak', show_label=False, width=0.5),
                        Item('del_peak', show_label=False, width=0.5),
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
                Item('start_time', 
                    editor=RangeEditor(low=0.0, high=1.0, 
                        low_name='start_time_low', high_name='start_time_high', 
                        format='%.1F', label_width=90, mode='slider',
                        ),
                    format_str='%.1F',  
                    label="Start Time"),
                Item('end_time', 
                    editor=RangeEditor(low=0.0, high=1.0, 
                        low_name='end_time_low', high_name='end_time_high', 
                        format='%.1F', label_width=90, mode='slider'
                    ), 
                    format_str='%.1F', 
                    label="End Time"), 
            ),
        ),
        width=500, 
        height=500, 
        resizable=True, 
        title="Pyramds Visualizer",
        )

    def get_total_time(self):
        if self.dfr == None:
            dt = 1.0
        else:
            dt = self.dfr.bin_data_parse.readout[-1]['timestamp'] - self.dfr.bin_data_parse.readout[0]['timestamp']
        return float(dt)

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

    def calc_histogram_data(self):
        hist_set = self.get_histogram_data_set()
        len_set = len(hist_set)

        # Calculate lower index
        lower_index = int( np.floor(len_set * self.start_time / self.end_time_high) )
        if lower_index == len_set:
            lower_index = lower_index - 1

        # Calculate upper index
        upper_index = int( np.floor(len_set * self.end_time   / self.end_time_high) )
        if upper_index == len_set:
            upper_index = upper_index - 1

        # Calculate histogram
        if lower_index == upper_index:
            hist = hist_set[upper_index]
        else:
            hist = hist_set[upper_index] - hist_set[lower_index]
        return hist

    def draw_plot(self):
        label = "Detector {0} {1} Spectrum".format(self.detector, self.spectrum)
        plot = make_spectrum_plot(self.chan, self.hist, self.pchn, self.peak, label)
        self.plot = plot

    def redraw_hist_plot(self):
        self.plot.plots['plot0'][0].index.set_data(self.chan)
        self.plot.plots['plot0'][0].value.set_data(self.hist)

    def redraw_peak_plot(self):
        self.plot.plots['plot1'][0].index.set_data(self.pchn)
        self.plot.plots['plot1'][0].value.set_data(self.peak)

    def redraw_plot(self):
        self.redraw_hist_plot()
        self.redraw_peak_plot()

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
         return 0.0

    def _start_time_high_default(self):
         return self.get_total_time()

    def _end_time_low_default(self):
         return 0.0

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

        # Calculate new hisogram
        hist = self.calc_histogram_data()
        self.hist = hist

        # Recalculate peak
        if 0 < len(self.pchn):
            peak = self.hist[self.pchn]
            self.peak = peak

        # Redraw points on plot
        self.redraw_plot()

    def _end_time_changed(self, old, new):
        self.start_time_high = new

        # Calculate new hisogram
        hist = self.calc_histogram_data()
        self.hist = hist

        # Recalculate peak
        if 0 < len(self.pchn):
            peak = self.hist[self.pchn]
            self.peak = peak

        # Redraw points on plot
        self.redraw_plot()

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
        self.redraw_peak_plot()

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

    def _save_detection_limits_fired(self):
        tsv = detection_limits_to_tsv(self.detection_limits)
        savefile = "{0}/{1}{2}_Detection_Limits.txt".format(self.save_directory, self.spectrum, self.detector)
        with open(savefile, 'w') as f:
            f.write(tsv)

if __name__ == "__main__":
    pv = PyramdsView()
    pv.configure_traits()

    # close the hdf5 file that is still open
    if pv.datafile != None:
        pv.datafile.close()
