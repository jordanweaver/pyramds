"""
PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)

Data parser for PIXIE List Mode binary data *.bin
The data contained in the .bin file is reformatted into an HDF5 file
that stores event information is a series of related table entries for quick
extraction of the necessary events used in spectra construction.

"""

from traits.api import HasTraits, Instance, Property, File, Int
from traitsui.api import View, Item, Button


class PyramdsModel(HasTraits):

    # Filecounter that tracks progress through file series
    file_counter = Int(1)
    
    # Filename for .bin file in series that is currently being processed
    bin_file_name = File()
        
    # Filename base (series number and extension removed)
    base_filepath = Property
    
    def _get_base_filepath(self):
        return self.bin_file_name[:8]

class PyramdsView(HasTraits):
    pmds = Instance(PyramdsModel)
    
    parse_button = Button(label="Parse Series")

    traits_view = View(Item(pmds.bin_file_name))
    
    # On filename change, update "Parse Series" button

if __name__ == '__main__':
    model = PyramdsModel()
    
    PyramdsWindow = PyramdsView(pmds=model)
    PyramdsWindow.configure_traits()

