from enthought.chaco.api import Plot, ArrayPlotData, ScatterInspectorOverlay, PlotLabel, PlotAxis, LabelAxis
from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom, ScatterInspector

from enthought.chaco.plot_graphics_context import PlotGraphicsContext

def make_spectrum_plot(chan, hist, pchn, peak, label="Spectrum Plot"):
    plotdata = ArrayPlotData(x=chan, y=hist, x2=pchn, y2=peak)

    container = Plot(
        plotdata,
        border_visible=True, 
        overlay_border=True,
        )

    scatter_plot = container.plot(("x", "y"), type="scatter", color="black")[0]
    container.plot(("x2", "y2"), type="scatter", color="red")

    # Add nice zooming and Panning
    container.tools.append(PanTool(container))

    zoom = ZoomTool(container, tool_mode="box", always_on=False)
    container.overlays.append(zoom)

    dragzoom = DragZoom(container, drag_button="right")
    container.tools.append(dragzoom)

    # Add scatter plot selector
    scatter_plot.tools.append(ScatterInspector(scatter_plot))
    scatter_overlay = ScatterInspectorOverlay(scatter_plot, 
        hover_color="yellow", 
        hover_marker_size=6,
        selection_marker_size=6,
        selection_color="limegreen",
        selection_outline_color="black",
        selection_line_width=1
        )
    scatter_plot.overlays.append(scatter_overlay)

    # Add the title at the top
    container.overlays.append(PlotLabel(label,
        component=container,
        font = "Times 20",
        overlay_position="top")
        )

    # Set x-axis
    container.x_axis.title = "Energy (MeV)"
    container.x_axis.title_font = "Times 16"
    container.x_axis.tick_label_font = "Times 12"

    # Set y-axis
    container.y_axis.title = "Counts"
    container.y_axis.title_font = "Times 16"
    container.y_axis.tick_label_font = "Times 12"

    # Get Title Spacing right
    container.padding_left = 65

    return container

def save_plot(plot, filename, width, height):
    orig_bounds = plot.outer_bounds
    plot.outer_bounds = [width, height]
    plot.do_layout(force=True)
    gc = PlotGraphicsContext((width, height), dpi=72)
    gc.render_component(plot)
    gc.save(filename)
    plot.outer_bounds = orig_bounds
