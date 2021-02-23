# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# HaselREC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HaselREC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HaselREC. If not, see <http://www.gnu.org/licenses/>.

def check_module(output_folder, site_code, probability_of_exceedance_num,
                 intensity_measures, n_gm, path_nga_folder):

    """

    This module is called when mode :code:`--check-NGArec` is specified.

    It identifies NGA-West2 records not already stored on the computer
    from the list of selected IDs.

    It requires to have run mode :code:`--run-selection` in advance since it reads
    in input the summary file created by mode :code:`--run-selection`

    It generates the file::

        missing_NGArec.txt

    which contains the list of missing NGA-West2 records.
    Example::

        Missing NGA-West2 IDrecords
        65
        68
        158
        159
        165
        170
        171
        179
        180
        181
        183
        184
        191
        266

    """

    import os
    import numpy as np
    import pandas as pd

    missing_file = output_folder + '/missing_NGArec.txt'
    missing_ID=[]
    for ii in np.arange(len(site_code)):
        site = site_code[ii]
        for jj in np.arange(len(probability_of_exceedance_num)):
            poe = probability_of_exceedance_num[jj]
            for im in np.arange(len(intensity_measures)):
                name = intensity_measures[im] + '-site_' + str(
                    site) + '-poe-' + str(poe)
                name_summary = (output_folder + '/' + name + '/' + name +
                                "_summary_selection.txt")
                summary = pd.read_csv(name_summary, sep=' ', skiprows=3)
                for i in np.arange(n_gm):
                    if summary.source[i] == 'NGA-West2':
                        start_string = 'RSN' + str(summary.recID_NGA[i]) +\
                                       '_'
                        if not os.path.isfile(
                                path_nga_folder + '/' + start_string +
                                '1.AT2') or not \
                                os.path.isfile(path_nga_folder + '/' +
                                               start_string + '2.AT2'):
                                    missing_ID.append(summary.recID_NGA[i])
    IDs=np.unique(np.asarray(missing_ID))
    with open(missing_file, "w") as f:
        f.write("Missing NGA-West2 IDrecords\n")
        for i in np.arange(len(IDs)):
            f.write("{}\n".format(IDs[i]))
    return
