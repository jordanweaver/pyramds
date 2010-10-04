from enthought.chaco.api import Plot, ArrayPlotData, ScatterInspectorOverlay
from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom, ScatterInspector

from enthought.chaco.plot_graphics_context import PlotGraphicsContext

def make_spectrum_plot(chan, hist, pchn, peak):
    plotdata = ArrayPlotData(x=chan, y=hist, x2=pchn, y2=peak)

    container = Plot(plotdata, border_visible=True, overlay_border=True)
    scatter_plot = container.plot(("x", "y"), type="scatter", color="black")[0]
    container.plot(("x2", "y2"), type="scatter", color="red")

    # Add nice zooming and Panning
    container.tools.append(PanTool(container))

    zoom = ZoomTool(container, tool_mode="box", always_on=False)
    container.overlays.append(zoom)

    dragzoom = DragZoom(container, drag_button="right")
    container.tools.append(dragzoom)

    scatter_plot.tools.append(ScatterInspector(scatter_plot))
    scatter_overlay = ScatterInspectorOverlay(scatter_plot, 
        hover_color="yellow", 
        hover_marker_size=6,
        selection_marker_size=6,
        selection_color="yellow",
        selection_outline_color="black",
        selection_line_width=1
        )
    scatter_plot.overlays.append(scatter_overlay)

    return container

def save_plot(plot, filename, width, height):
    orig_bounds = plot.outer_bounds
    plot.outer_bounds = [width, height]
    plot.do_layout(force=True)
    gc = PlotGraphicsContext((width, height), dpi=72)
    gc.render_component(plot)
    gc.save(filename)
    plot.outer_bounds = orig_bounds
