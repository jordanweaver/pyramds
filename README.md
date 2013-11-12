# PYRAMDS (Python for Radioisotope Analysis & Multidetector Suppression) #

PYRAMDS is basically a data parser for PIXIE Gamma Detector List Mode binary data. The data contained in PIXIE's .bin files can be reformatted into an HDF5 file that stores event information in a series of related table entries for quick extraction of the necessary events used in spectra construction.

Focus is placed on the ability to post-process gamma detector data for producing various flavors of spectra, such as Compton suppressed or Gamma-Gamma coincidence.


INCLUDE DEPENDENCIES (ENAML, TRAITS, TABLES)

