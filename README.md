# Werts-veusz-ProcSpec-reader-plugin
Veusz Import Plugin for reading Ocean Optics ProcSpec files (written by Martinus Werts)

# Description
- Veusz Import Plugin for reading Ocean Optics ProcSpec files
- Note: the ProcSpec files should contain a 'Processed' spectrum
- (e.g. an absorbance spectrum, or a background-corrected emission spectrum)
- Files containing only a raw spectrum will not be read, resulting
- in an error message.
- 
- by Martinus Werts, CNRS, ENS Rennes, France
- v140416   first release, 

# about mwave class
- mwave class, enables storing abscissa and ordinate data simultaneously in a single object (such as wavelength and intensity for a spectrum) allowing for certain operations to be more easily achieved this module has been copied from the original source file in order to have this plugin as a monolithic file

# What's New
- new mwave object working!
- to do: thorough test baseline correct
