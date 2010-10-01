#! /usr/bin/env python
import os
import time
import subprocess

import numpy 
import tables
import scipy

# ETS imports
from enthought.traits.api import HasTraits, Instance, Str, Enum, Any, List, Array, Dict, Button, NO_COMPARE, Float
from enthought.traits.ui.api import View, Item, VGroup, HGroup, EnumEditor, TableEditor, TabularEditor
from enthought.traits.ui.ui_editors.array_view_editor import ArrayViewEditor
from enthought.traits.ui.menu import Action
from enthought.chaco.api import Plot, ArrayPlotData
from enthought.enable.api import Component, ComponentEditor, Window

from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, \
         TraitsTool, DragZoom

np = numpy
tb = tables

class PyramdsView(HasTraits):
    plot = Instance(Plot)

    # Database info
    filename = Str("Co60_92710-.h5")
    datafile = Any()
    dfr      = Any()

    # plot data
    chan = Array
    hist = Array
    peak = Array

    title = Str("")

#    save_directory = Str("vizfigs/")
#    save_figure = Button(label="Save Figure")
#    save_table = Button(label="Save Table")

    # UI Stuff
    start_time = Float
    end_time = Float

    traits_view = View(
        Item('plot', editor=ComponentEditor(), show_label=False),
        HGroup(
            Item('start_time', format_str='%.2F',  label="Start Time"), 
            Item('end_time', format_str='%.2F', label="End Time"), 
        ),
        width=500, 
        height=500, 
        resizable=True, 
        title="Pyramds Visualizer",
        )

    def draw_plot(self):
        pass

    # Set default Values
    def _datafile_default(self):
        return tb.openFile(self.filename, "r")

    def _dfr_default(self):
        return self.datafile.root

    def _chan_default(self):
        return np.arange(len(self.dfr.spectra.norm1_spec))

    def _hist_default(self):
        return self.dfr.spectra.norm1_spec.read()

    def _peak_default(self):
        return np.array([])

    def _plot_default(self):
        # Add Data
        plotdata = ArrayPlotData(x=self.chan, y=self.hist, y2=self.peak)
        plot = Plot(plotdata)
        plot.plot(("x", "y"), type="scatter", color="black")
        plot.plot(("x", "y2"), type="scatter", color="red")

        # Add nice zooming and Panning
        plot.tools.append(PanTool(plot))
        zoom = ZoomTool(plot, tool_mode="box", always_on=False)
        plot.overlays.append(zoom)
        dragzoom = DragZoom(plot, drag_button="right")
        plot.tools.append(dragzoom)

        return plot

    def _start_time_default(self):
        return 0.0

    def _end_time_default(self):
        return self.dfr.bin_data_parse.readout[-1]['timestamp'] - self.dfr.bin_data_parse.readout[0]['timestamp']

if __name__ == "__main__":
    pv = PyramdsView()
#    pv.draw_plot()
    pv.configure_traits()
#    pv.edit_traits()
