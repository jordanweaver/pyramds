# PYRAMDS #

PYRAMDS (Python for Radioisotope Analysis & Multi-Detector Suppression) is at most a data parser for PIXIE Gamma Detector List Mode binary data. The data contained in PIXIE's .bin files can be reformatted into an HDF5 file. This file stores event information in a series of related table entries for quick extraction of the necessary events to be used in spectra construction.

Focus is placed on the ability to post-process gamma detector data for producing various flavors of spectra, such as Compton suppressed or Gamma-Gamma coincidence.

This is most definitely a work in progress. And by "this" I mean both the software and my ability to write it.

------------------------------------

## Requirements ##
Regarding older versions of the following packages, I have not investigated how "required" the versions I've listed are. Instead, think of these as the tools I used to build it. I used the [Enthought Canopy distribution](https://www.enthought.com/products/canopy/), supplemented with more recent versions of certain packages.

* Python 2.7
* Numpy
* PyTables
* Traits
* [Enaml](https://github.com/nucleic/enaml)

## Installation & Usage ##
After installing dependencies, clone the repo and run: `python pyramds.py`

/Include screenshot of GUI

## Configuration for Detector Systems ##

KNOWN LIMITATIONS/BUGS

CREDITS/ACKNOWLEDGEMENTS

CONTACT INFORMATION FOR QUESTIONS, BUGS REPORTS
