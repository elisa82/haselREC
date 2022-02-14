# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# haselREC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# haselREC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with haselREC. If not, see <http://www.gnu.org/licenses/>.

def compute_conditioning_value(rlz, intensity_measures, site, poe, num_disagg,
                               probability_of_exceedance, num_classical,
                               path_results_disagg, investigation_time,
                               path_results_classical, meanMag_disagg, 
                               meanDist_disagg, hazard_value, hazard_mode):
    """
    Reads 2 output files ('.csv') from OpenQuake: the file with disaggregation
    results and the map with hazard values and computes the IM value at which to
    condition the CS, along with the mean magnitude and distance from the
    disaggregation analysis.
    """
    import pandas as pd
    import numpy as np

    if(hazard_mode==0):

        # Get the name of the disaggregation file to look in
        disagg_results = 'rlz-' + str(rlz) + '-' + intensity_measures + '-sid-' + \
                         str(site) + '-poe-' + str(poe) + '_Mag_Dist_' + \
                         str(num_disagg) + '.csv'

        probability_of_exceedance=np.float64(probability_of_exceedance)

        selected_column = intensity_measures + '-' + str(probability_of_exceedance)
        file_with_oq_acc_value = 'hazard_map-mean_' + str(num_classical) + '.csv'

        # Retrieve disaggregation results
        df = pd.read_csv(''.join([path_results_disagg, '/', disagg_results]),
            skiprows=1)
        df['rate'] = -np.log(1 - df['poe']) / investigation_time
        df['rate_norm'] = df['rate'] / df['rate'].sum()
        # mode = df.sort_values(by='rate_norm', ascending=False)[0:1]
        mean_mag = np.sum(df['mag'] * df['rate_norm'])
        mean_dist = np.sum(df['dist'] * df['rate_norm'])

        # Retrieve conditioning value
        df = pd.read_csv(''.join(
            [path_results_classical, '/', file_with_oq_acc_value]), skiprows=1)
        output_oq = df[selected_column]
        im_star = output_oq[site]

        dist = np.array([mean_dist])
        mag = mean_mag

    else:
        mag = meanMag_disagg
        dist = np.array([meanDist_disagg])
        im_star = hazard_value

    return im_star, dist, mag
