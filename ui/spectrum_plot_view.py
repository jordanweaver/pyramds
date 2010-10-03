from enthought.chaco.api import Plot, ArrayPlotData
from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, TraitsTool, DragZoom

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
