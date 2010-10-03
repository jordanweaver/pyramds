from enthought.chaco.api import Plot, ArrayPlotData
from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom

from enthought.chaco.plot_graphics_context import PlotGraphicsContext

def make_spectrum_plot(chan, hist, peak):
    plotdata = ArrayPlotData(x=chan, y=hist, y2=peak)

    container = Plot(plotdata)
    container.plot(("x", "y"), type="scatter", color="black")
    container.plot(("x", "y2"), type="scatter", color="red")
              
    # Add nice zooming and Panning
    container.tools.append(PanTool(container))

    zoom = ZoomTool(container, tool_mode="box", always_on=False)
    container.overlays.append(zoom)

    dragzoom = DragZoom(container, drag_button="right")
    container.tools.append(dragzoom)

    return container

def save_plot(plot, filename, width, height):
    orig_bounds = plot.outer_bounds
    plot.outer_bounds = [width, height]
    plot.do_layout(force=True)
    gc = PlotGraphicsContext((width, height), dpi=72)
    gc.render_component(plot)
    gc.save(filename)
    plot.outer_bounds = orig_bounds
