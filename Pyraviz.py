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

from function_lib import calc_det_limit, calc_det_limit_sel

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

    # Signal Lookup
    sig_lookup = Dict({
        "1": {}, 
        "2": {}, 
        })
    en_coeff_1 = Array(dtype=float)
    en_coeff_2 = Array(dtype=float)
    fwhm_coeff_1 = Array(dtype=float)
    fwhm_coeff_2 = Array(dtype=float)

    # Save Information
    save_directory = Directory(".")
    save_figure = Button(label="Save Figure")

    # plot data
    chan = Array(dtype=int)
    hist = Array(dtype=int)
    pchn = Array(dtype=int)
    peak = Array(dtype=int)

    title = Str("")
    linlog_toggle = Enum(["Linear", "Log"])

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
        'lookup_selected': [],
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
                spring,
                Item('linlog_toggle', show_label=False),
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

    def load_sig_lookup(self):
        sig_lookup = {"1": {}, "2": {}}

        for n in ["1", "2"]:
            sig_table = getattr(self.dfr.sig_lookup, "det{0}_sig".format(n))
            # init isolists
            for row in sig_table:
                sig_lookup[n][row['name']] = []

            # Add data
            for row in sig_table:
                sig_lookup[n][row['name']].append( (row['LM'], row['RM']) )

        self.sig_lookup = sig_lookup

        self.en_coeff_1 = self.dfr.sig_lookup.en_coeff_1.read()
        self.en_coeff_2 = self.dfr.sig_lookup.en_coeff_2.read()
        self.fwhm_coeff_1 = self.dfr.sig_lookup.fwhm_coeff_1.read()
        self.fwhm_coeff_2 = self.dfr.sig_lookup.fwhm_coeff_2.read()


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
        en_coeff = getattr(self, "en_coeff_{0}".format(self.detector))
        plot = make_spectrum_plot(self.chan, self.hist, self.pchn, self.peak, en_coeff, label, self.linlog_toggle)
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

    def reset_isotope(self):
        isotope_enum = sorted(self.sig_lookup[self.detector].keys())
        isotope = isotope_enum[0]
        self.isotope_enum = isotope_enum
        self.isotope = isotope

    def reset_peaknum(self):
        peaknum_enum = [str(i) for i in range(1, len(self.sig_lookup[self.detector][self.isotope]) + 1)]
        peaknum = peaknum_enum[0]
        self.peaknum_enum = peaknum_enum
        self.peaknum = peaknum

    def time_changed(self):
        # Calculate new hisogram
        hist = self.calc_histogram_data()
        self.hist = hist

        # Recalculate peak
        if 0 < len(self.pchn):
            peak = self.hist[self.pchn]
            self.peak = peak

        # Redraw points on plot
        self.redraw_plot()

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
        en_coeff = getattr(self, "en_coeff_{0}".format(self.detector))
        plot = make_spectrum_plot(self.chan, self.hist, self.pchn, self.peak, en_coeff, scale=self.linlog_toggle)
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

    def _en_coeff_1_default(self):
        return np.array([0.0, 1.0], dtype=float)

    def _en_coeff_2_default(self):
        return np.array([0.0, 1.0], dtype=float)

    def _fwhm_coeff_1_default(self):
        return np.array([1.0, 1.0, 1.0], dtype=float)

    def _fwhm_coeff_2_default(self):
        return np.array([1.0, 1.0, 1.0], dtype=float)

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

        # Load signal lookup and add to drop-downs
        self.load_sig_lookup()
        self.reset_isotope()
        self.reset_peaknum()

        # Redraw plot
        self.draw_plot()


    def _start_time_changed(self, old, new):
        self.end_time_low = new
        self.time_changed()

    def _end_time_changed(self, old, new):
        self.start_time_high = new
        self.time_changed()

    def _isotope_changed(self):
        self.reset_peaknum()

    def _detector_changed(self):
        self.load_histogram_data()
        self.draw_plot()

        self.reset_isotope()
        self.reset_peaknum()

    def _spectrum_changed(self):
        self.load_histogram_data()
        self.draw_plot()

    def _linlog_toggle_changed(self, old, new):
        if new == "Log":
            self.plot.value_scale = 'log'

        if new == "Linear":
            self.plot.value_scale = 'linear'

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
                en_coeff = getattr(self, "en_coeff_{0}".format(self.detector))
                fwhm_coeff = getattr(self, "fwhm_coeff_{0}".format(self.detector))
                ld, lc = calc_det_limit_sel(self.detector, si, en_coeff, fwhm_coeff, self.hist)
                self.detection_limits['clicked_selected'][si] = [ld, lc]

        # Remove old detection limits
        for key in self.detection_limits['clicked_selected'].keys():
            if key not in self.index_selections:
                del self.detection_limits['clicked_selected'][key]

        # Redraw the detection limits display
        self.draw_detection_limits()
        
    #
    # Buttons Fired methods
    #

    def _add_peak_fired(self):
        markers = self.sig_lookup[self.detector][self.isotope][int(self.peaknum)-1]
        ld, lc = calc_det_limit(markers[0], markers[1], self.hist)

        # Set detection limits
        row = [self.isotope, self.peaknum, ld, lc]
        self.detection_limits['lookup_selected'].append(row)

        # Redraw the detection limits display
        self.draw_detection_limits()

        # Add peak to plot
        new_pchn = np.arange(markers[0], markers[1] + 1, dtype=int)
        pchn = np.append(self.pchn, new_pchn)
        peak = self.hist[pchn]

        self.pchn = pchn
        self.peak = peak

        # Redraw the plot
        self.redraw_peak_plot()        

    def _del_peak_fired(self):
        del_index = None
        for n in range(len(self.detection_limits['lookup_selected'])):
            row = self.detection_limits['lookup_selected'][n]
            if (row[0] == self.isotope) and (row[1] == self.peaknum):
                del_index = n
                break

        if del_index != None:
            del self.detection_limits['lookup_selected'][del_index]

            # Redraw the detection limits display
            self.draw_detection_limits()

            # remove peak from plot
            markers = self.sig_lookup[self.detector][self.isotope][int(self.peaknum)-1]
            rm_chn = np.arange(markers[0], markers[1] + 1, dtype=int)
            new_pchn = np.array([chn for chn in self.pchn if chn not in rm_chn], dtype=int)
            if 0 < len(new_pchn):
                peak = self.hist[new_pchn]
            else:
                peak = np.array([])

            self.pchn = new_pchn
            self.peak = peak

            # Redraw the plot
            self.redraw_peak_plot()
    
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
