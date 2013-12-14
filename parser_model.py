# PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)
#
# Author: Jordan Weaver

from tables import Float32Col, Int32Col, IsDescription

# Internal Imports
from parser_setup import PyramdsBase

# Setup PyTables metaclasses for use in Table constructor
class GammaEvent(IsDescription):

    # Energy reading from Pixie detector channels 0, 1, 2
    energy_0 = Int32Col(pos=0)
    energy_1 = Int32Col(pos=1)
    energy_2 = Int32Col(pos=2)

    # Time differences between detector events
    # (e.g. "deltaT_01" is between energy_0 and energy_1)
    deltaT_01 = Float32Col(pos=3)
    deltaT_02 = Float32Col(pos=4)
    deltaT_12 = Float32Col(pos=5)

    timestamp = Float32Col(pos=6)

class AggEvent1(IsDescription):
    energy = Int32Col(pos=0)
    timestamp = Float32Col(pos=1)

class AggEvent2(IsDescription):
    energy_1 = Int32Col(pos=0)
    energy_2 = Int32Col(pos=1)
    timestamp = Float32Col(pos=2)

class PyramdsParser(PyramdsBase):

