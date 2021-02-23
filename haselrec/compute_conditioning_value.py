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

def compute_conditioning_value(rlz, intensity_measures, site, poe, num_disagg,
                               probability_of_exceedance, num_classical,
                               path_results_disagg, investigation_time,
                               radius_dist, radius_mag, path_results_classical):

    import pandas as pd
    import numpy as np

    # Get the name of the disaggregation file to look in
    disagg_results = 'rlz-' + str(rlz) + '-' + intensity_measures + '-sid-' + \
                     str(site) + '-poe-' + str(poe) + '_Mag_Dist_' + \
                     str(num_disagg) + '.csv'

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
    allowed_recs_d = [mean_dist - radius_dist, mean_dist + radius_dist]
    allowed_recs_mag = [mean_mag - radius_mag, mean_mag + radius_mag]

    # Retrieve conditioning value
    df = pd.read_csv(''.join(
        [path_results_classical, '/', file_with_oq_acc_value]), skiprows=1)
    output_oq = df[selected_column]
    im_star = output_oq[site]

    return im_star, allowed_recs_d, allowed_recs_mag
