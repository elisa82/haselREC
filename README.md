# OpenInsel
Open-source tool for INtegrated hazard analysis and ground motion record SELection (OpenInsel)

It perfomes a ground motion record selection following a target conditional spectrum, using the OpenQuake libraries and disaggregation results.

Main contributors:
* Elisa Zuccolo - EUCENTRE Foundation, Italy
* Gerard J. O'Reilly - Scuola Universitaria Superiore IUSS Pavia, Italy

# Dependencies
OpenInsel requires the following dependencies:
NumPy
Pandas
Scipy
OS/Sys/Glob
ObsPy
Matplotlib
Openquake.hazardlib

# Usage
* cd OI-selection
* python select_accelerograms.py $filename.ini

# License
Copyright (c) 2020, OpenInsel Developers
OpenInsel is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3.0 of the License, or (at your option) any later version.
You should have received a copy of the GNU General Public License with this download. If not, see http://www.gnu.org/licenses/

# Disclaimer
OpenInsel is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
The authors of the software assume no liability for use of the software.

# Known Issues
A list of issues that need to be fixed:
* Need to download complete NGA West 2 database and keep original filenames to make adding/removing easier. Should try to find a way to download the whole database

# Potential Improvements
A list of things that could be improved:
* The need to directly state the IM and period info could be removed and taken from the chosen disaggregation outputs - ???
* Possibility to directly specify the strong-motion database to be used - DONE
* To increase the database of GMPE -DONE
* Adding the possibility to directly execute a disaggregation analysis - NO

# Main References

* Baker JW, Lee C. (2018) An Improved Algorithm for Selecting Ground Motions to Match a Conditional Spectrum. Journal of Earthquake Engineering 22(4): 708–723. 
* Kaklamanos, J., L. G. Baise, and D. M. Boore (2011). Estimating unknown input parameters when implementing the NGA ground-motion prediction equations in engineering practice, Earthquake Spectra 27(4): 1219-1235.
* Kohrangi M, Bazzurro P, Vamvatsikos D, Spillatura A. (2017) Conditional spectrum-based ground motion record selection using average spectral acceleration. Earthquake Engng Struct. Dyn. 46: 1667–1685.
