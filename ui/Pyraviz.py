#! /usr/bin/env python
import os
import time
import subprocess

import numpy 
import tables
import scipy

# ETS imports
from enthought.traits.api import HasTraits, Instance, Str, Enum, Any, List, Array, Dict, Button, NO_COMPARE, Float, \
    File, Directory

from enthought.traits.ui.api import View, Item, Group, VGroup, HGroup, EnumEditor, TableEditor, TabularEditor, spring
from enthought.traits.ui.ui_editors.array_view_editor import ArrayViewEditor
from enthought.traits.ui.menu import Action

from enthought.chaco.api import Plot, ArrayPlotData
from enthought.enable.api import Component, ComponentEditor, Window

from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom

from spectrum_plot_view import make_spectrum_plot, save_plot

np = numpy
tb = tables

size = (500, 500)

class PyramdsView(HasTraits):
    plot = Instance(Component)

    # Database info
    filename = File
    datafile = Any()
    dfr      = Any()

    # Save Information
    save_directory = Directory(".")
    save_figure = Button(label="Save Figure")
#    save_table = Button(label="Save Table")
    

    # plot data
    chan = Array
    hist = Array
    peak = Array

    title = Str("")

    detector = Enum(["1", "2"])
    spectrum = Enum(["Normal", "Compton", "Gamma-Gamma"])
    spectrum_set_names = {
        "Normal":      "norm",
        "Compton":     "compt",
        "Gamma-Gamma": "gg",
        }

    # UI Stuff
    start_time = Float
    end_time = Float

    traits_view = View(
        VGroup(
            Item('filename', label="Data File"),
            Item('plot', editor=ComponentEditor(size=size), show_label=False), 
            HGroup(
                Item('detector'), 
                spring,
                Item('spectrum'), 
                spring,
                Item('start_time', format_str='%.2F',  label="Start Time"),
                spring,
                Item('end_time', format_str='%.2F', label="End Time"), 
            ),
            HGroup(
                Item('save_directory', width=400, label="Save Dir"),
                spring,
                Item('save_figure', show_label=False),
            ),
        ),
        width=500, 
        height=500, 
        resizable=True, 
        title="Pyramds Visualizer",
        )

    def load_histogram_data(self):
        hist_set = getattr(self.dfr.spectra, "{0}{1}_spec".format(self.spectrum_set_names[self.spectrum], self.detector))
        hist = hist_set.read()
        chan = np.arange(len(hist))
        peak = np.array([])

        self.chan = chan
        self.hist = hist
        self.peak = peak

    def draw_plot(self):
        pass

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
            return np.array([])
        else:
            return np.arange(len(self.hist))

    def _hist_default(self):
        if self.dfr == None:
            return np.array([])
        else:
            hist_set = getattr(self.dfr.spectra, "{0}{1}_spec".format(self.spectrum_set_names[self.spectrum], self.detector))
            hist = hist_set.read()
            return hist

    def _peak_default(self):
        return np.array([])

    def _plot_default(self):
        plot = make_spectrum_plot(self.chan, self.hist, self.peak)
        return plot

    def _start_time_default(self):
        return 0.0

    def _end_time_default(self):
        if self.dfr == None:
            return 0.0
        else:
            return self.dfr.bin_data_parse.readout[-1]['timestamp'] - self.dfr.bin_data_parse.readout[0]['timestamp']

    def _detector_default(self):
        return "1"

    def _spectrum_default(self):
        return "Normal"

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

        self.start_time = 0.0
        self.end_time = self.dfr.bin_data_parse.readout[-1]['timestamp'] - self.dfr.bin_data_parse.readout[0]['timestamp']

        plot = make_spectrum_plot(self.chan, self.hist, self.peak)
        self.plot = plot

    def _detector_changed(self):
        self.load_histogram_data()
        plot = make_spectrum_plot(self.chan, self.hist, self.peak)
        self.plot = plot

    def _spectrum_changed(self):
        self.load_histogram_data()
        plot = make_spectrum_plot(self.chan, self.hist, self.peak)
        self.plot = plot

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
