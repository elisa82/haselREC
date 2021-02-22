# HaselREC 
HAzard-based SELection of RECords

It perfoms a ground motion record selection following a target conditional spectrum, using the OpenQuake libraries and hazard results.
The code is described in:
*Zuccolo E, Poggi V, O'Reilly G, Monteiro R (2021). HaselREC: an open-source ground motion record selection tool bridging seismic hazard and structural analyses. Submitted to SDEE

Main contributors:
* Elisa Zuccolo - EUCENTRE Foundation, Italy, elisa.zuccolo@eucentre.it
* Gerard J. O'Reilly - Scuola Universitaria Superiore IUSS Pavia, Italy, gerard.oreilly@iusspavia.it

# Dependencies
HaselREC has the following dependencies:

 * NumPy
 * Pandas
 * Scipy
 * ObsPy
 * Matplotlib
 * Openquake.hazardlib
 * shakelib.conversions.imc.boore_kishida_2017 (Optional. Only of the input intensity measure component is the larger between the two horizontal components).

# Installation
Add the lib folder to PYTHONPATH. For Linux open the file ~/.bashrc in your text editor and add the following line at the end:
```
export PYTHONPATH=/path/to/HaselREC/lib
```

##Database
- All NGA-West2 recordings have to be stored in a folder and renamed as:
*RSN#NUM_1.AT2* (1st horiz comp)
*RSN#NUM_2.AT2* (2nd horiz comp)
*RSN#NUM_3.AT2* (vertical component)
#NUM is the record sequence number of NGA recordings

- ESM recordings can be stored in advance or automatically downloaded from internet using a token file. To obtain the token file you need at first to register at: https://esm-db.eu/ and then run the command:
```
curl -X POST -F 'message={"user_email": "email","user_password": "password"}' "https://esm-db.eu/esmws/generate-signed-message/1/query" > token.txt
```

# Usage
```
haselrec.py $filename.ini [option]
```
Possible [option] are: 
- --run-complete 
- --run-selection
- --run-scaling
- --check-NGArec

# License
Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
HaselREC is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
HaselREC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with HaselREC. If not, see <http://www.gnu.org/licenses/>.

# Disclaimer
HaselREC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
The authors of the software assume no liability for use of the software.

# Potential Improvements
A list of things that could be improved:
* Use of original names for recordings from the NGA-West2 database 
* Computation of exact CS 
