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
    import os
    import numpy as np
    import pandas as pd

    missing_file = output_folder + '/missing_NGArec.txt'
    with open(missing_file, "w") as f:
        f.write("Missing NGA-West2 IDrecords\n")
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
                            # comp1 = ['1', 'EW', '-W', '-E']
                            # comp2 = ['2', 'NS', '-N']
                            # import glob
                            # exist1 = 0
                            # exist2 = 0
                            # for k in np.arange(len(comp1)):
                            #    end_string=comp1[k]+'.AT2'
                            #    if glob.glob(path_NGA_folder+'/'+start_string+
                            #    '*'+end_string):
                            #        exist1=1
                            # for k in np.arange(len(comp2)):
                            #    end_string=comp2[k]+'.AT2'
                            #    if glob.glob(path_NGA_folder+'/'+start_string+
                            #    '*'+end_string):
                            #        exist2=1
                            # if exist1==0 or exist2==0:
                            #    f.write("{}\n".format(summary.recID_NGA[i]))
                            if not os.path.isfile(
                                    path_nga_folder + '/' + start_string +
                                    '1.AT2') or not \
                                    os.path.isfile(path_nga_folder + '/' +
                                                   start_string + '2.AT2'):
                                f.write("{}\n".format(summary.recID_NGA[i]))
    return
