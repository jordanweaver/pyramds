"""
PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression)

Data parser for PIXIE List Mode binary data *.bin
The data contained in the .bin file is reformatted into an HDF5 file
that stores event information is a series of related table entries for quick
extraction of the necessary events used in spectra construction.

"""

from datetime import datetime

from traits.api import (HasTraits, Instance, Property, File, Int, Button)
from traitsui.api import View, Item
from tables import IsDescription, Int32Col, Float32Col

# Setup PyTables metaclasses for use in Table constructor
class GammaEvent(IsDescription):
    
    # Energy reading from Pixie detector channels 0, 1, 2
    energy_0  = Int32Col(pos=0)
    energy_1  = Int32Col(pos=1)
    energy_2  = Int32Col(pos=2)
    
    # Time differences between detector events 
    # (e.g. "deltaT_01" is between energy_0 and energy_1)
    deltaT_01 = Float32Col(pos=3)
    deltaT_02 = Float32Col(pos=4)
    deltaT_12 = Float32Col(pos=5)
    
    timestamp = Float32Col(pos=6)

class AggEvent1(IsDescription):
    energy    = Int32Col(pos=0)
    timestamp = Float32Col(pos=1)

class AggEvent2(IsDescription):
    energy_1  = Int32Col(pos=0)
    energy_2  = Int32Col(pos=1)
    timestamp = Float32Col(pos=2)

class PyramdsParser(HasTraits):

    # Filename for .bin file selected in UI
    selected_bin_file = File()
        
    # Filename base (series number and extension removed)
    file_series = Property
    
    # Filecounter that tracks progress through file series
    file_counter = Int(1)
    # File path currently being processed (no extension)
    active_file_path = file_series + '%04d'.format(file_counter)
    
    # Only initialize the buffer counter before the entire run... not each file.
    buffer_no = 0

    def get_bin_info(self):
        
        active_file_ifm = (self.file_series + '%04d' + '.ifm') % self.file_counter
        
        with open(active_file_ifm, 'rU') as finfo:
    
            info_str_list = finfo.readlines()
    
            date_str = info_str_list[1][23:-2]
            date_format = '%I:%M:%S %p %a, %b %d, %Y'
            run_start_time = datetime.strptime(date_str, date_format)
    
            # Total time = real time
            total_time = info_str_list[6].split()[3]
    
            # Live time for each detector channel
            live_time = [info_str_list[9 + channel].split()[2] for channel in range(4)]

            self.times = {
                'start' :   run_start_time,
                'total' :   total_time,
                'live'  :   live_time
            }
    
            self.bufheadlen = int(info_str_list[33].split()[1])
            self.eventheadlen = int(info_str_list[34].split()[1])
            
            # Due to bug in PIXIE IGOR Software, explicity set head length
            self.chanheadlen = 2 #int(info_str_list[35].split()[1])
        
    
    def _get_file_series(self):
        return self.selected_bin_file[:8]

class PyramdsView(HasTraits):
    Parser = Instance(PyramdsParser)
    
    filename = File()  
    parse_button = Button(label="Parse Series")

    traits_view = View(Item('filename', label = 'BIN File'),
        width = 700,
        height = 700,
        resizable = True,
        title = "PYRAMDS",
        )
    
    # On filename change, update parser model and pull new .ifm stats
    def _filename_changed(self, new):
       self.Parser.selected_bin_file = new
       # self.parser_view.get_bin_info()

if __name__ == '__main__':
    parser = PyramdsParser()
    
    PyramdsWindow = PyramdsView(Parser=parser)
    PyramdsWindow.configure_traits()

