![haselREC Logo](docs/logo.png)

# haselREC - HAzard-based SELection of RECords

The code is described in:
Zuccolo E, Poggi V, O'Reilly GJ, Monteiro R (2021). 
haselREC: an automated open-source ground motion record selection tool
compatible with OpenQuake. Submitted to SDEE

Contributors:
* Elisa Zuccolo - EUCENTRE Foundation, Italy, elisa.zuccolo@eucentre.it
* Gerard J. O'Reilly - Scuola Universitaria Superiore IUSS Pavia, Italy, 
gerard.oreilly@iusspavia.it
* Andrea Francia - andrea@andreafrancia.it 

# Dependencies
haselREC has the following dependencies:

 * numpy
 * pandas
 * scipy
 * obspy
 * matplotlib
 * openquake.hazardlib
 * shakelib.conversions.imc.boore_kishida_2017 (Optional. Only of the input 
 intensity measure component is the larger between the two horizontal 
 components).

# Installation
Add the lib folder to PYTHONPATH. For Linux open the file `~/.bashrc` in your 
text editor and add the following line at the end:
```
export PYTHONPATH=/path/to/haselREC
```
# How to make haselREC to automatically download ESM recordings
ESM recordings can be stored in advance or automatically downloaded from 
internet using a token file (token.txt). To obtain the token file you need 
at first to register at: https://esm-db.eu/ and then you can run the command:

```
curl -X POST -F 'message={"user_email": "email","user_password": "password"}' "https://esm-db.eu/esmws/generate-signed-message/1/query" > token.txt
```

# Documentation
The documentation of the code can be generated using Sphinx (v1.8.5).
To run the documentation type:
```
make html
```
The documentation will be generated in the *build* folder

# Demos
Some demos can be found in the *demo* folder. See the README file in the *demo*
folder



# License
Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
haselREC is free software: you can redistribute it and/or modify it under the 
terms of the GNU Affero General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later 
version. See the GNU Affero General Public License for more details. You should
have received a copy of the GNU Affero General Public  License along with 
haselREC. If not, see <http://www.gnu.org/licenses/>.

# Disclaimer
haselREC is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
The authors of the software assume no liability for use of the software.

# Potential Improvements
A list of things that could be improved:
* Use of original names for recordings from the NGA-West2 database 
* Computation of exact CS 
