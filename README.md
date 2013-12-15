# PYRAMDS #

PYRAMDS (Python for Radioisotope Analysis & Multi-Detector Suppression) is at most a data parser for PIXIE Gamma Detector List Mode binary data. The data contained in PIXIE's .bin files can be reformatted into an HDF5 file. This file stores event information in a series of related table entries for quick extraction of the necessary events to be used in spectra construction.

Focus is placed on the ability to post-process gamma detector data for producing various flavors of spectra, such as Compton suppressed or Gamma-Gamma coincidence.

This is most definitely a work in progress. And by "this" I mean both the code and my ability to write it.

------------------------------------

## Requirements ##
Regarding older versions of the following packages, I have not investigated how "required" the versions I've listed are. Instead, think of these as the tools I used to build it at the time. I used the [Enthought Canopy distribution](https://www.enthought.com/products/canopy/).

* Python 2.7
* Numpy
* PyTables
* Traits
* TraitsUI

## Installation & Usage ##
After installing dependencies, clone the repo and run: `python pyramds.py`

/Include screenshot of GUI

## Configuration for Detector Systems ##
The application still requires a significant amount of configuration on the part of the researcher. This configuration has mostly to do with collecting the necessary calibration/resolution spectra to characterize the detector performance, and then saving these parameters in the PYRAMDS application for use in the spectrum exporter.

## Known Limitations/Bugs ##
The application still lacks most (or any) types of exception handling and is probably still rather fragile in how one operates the parser. I have not tested the limits of this and still recommend following a fairly routine recipe/workflow when operating the GUI parser.

## Credits & Acknowledgements ##
The author would like to thank Anthony Scopatz for his initial push many years ago for using Python to do my doctoral research, of which the legacy branch is the majority of that work. His help with quickly hacking out the first GUI, Pyraviz, over the course of a few days has helped immensely in having a working example to learn Traits/TraitsUI in the context of this work.

Second, thanks also goes to Jonathan March for his discussions, guidance, and motivation for trying out new things and using the right tool for the job, even if I felt it was overkill, this has been very educational.

CONTACT INFORMATION FOR QUESTIONS, BUGS REPORTS
