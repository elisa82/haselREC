# OpenInsel
Open-source tool for INtegrated hazard analysis and ground
motion record SELection (OpenInsel)

It utilises the OpenQuake hazard computation engine to perform a
probabilistic seismic hazard analysis (PSHA) for a single or number of sites and
ground motion intensity measures, followed by a selection of ground motion
records following a target conditional spectrum.

Main contributors:
* Elisa Zuccolo - EUCENTRE Foundation, Italy
* Gerard J. O'Reilly - Scuola Universitaria Superiore IUSS Pavia, Italy


# Known Issues
A list of issues that need to be fixed:
* Need to download complete NGA West 2 database and keep original filenames to make adding/removing easier. Should try to find a way to download the whole database

# Potential Improvements
A list of things that could be improved:
* The need to directly state the IM and period info could be removed and taken from the chosen disaggregation outputs
* Possibility to directly specify the strong-motion database to be used
* To increase the database of GMPE
* Adding the possibility to directly execute a disaggregation analysis


# References
