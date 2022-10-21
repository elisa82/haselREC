#!/usr/bin/env python3

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

# %% Import libraries
# Standard built-in libraries

"""
haselREC (HAzard-based SELection of RECords)
It is a useful open-source tool for OpenQuake users, able to select and scale
recorded accelerograms to be used for dynamic analyses.

It is described in:
Zuccolo E, Poggi V, O'Reilly G, Monteiro R (2021).
haselREC: an automated open-source ground motion record selection tool
compatible with OpenQuake. Under review

haselREC can be launched with the following command::

    python3 -m haselrec <input_file> <mode>

Four modes are permitted:

    - :code:`--run-selection`: it performs record selection only
    - :code:`--run-scaling`: it performs record scaling only (requires to have run
       mode :code:`--run-selection` in advance)
    - :code:`--run-complete`: it performs both record selection and scaling
    - :code:`--check-rec`: it identifies NGA-West2 and KiK-net record IDs not already stored
       on the computer (it requires to have run mode :code:`--run-selection` in
       advance)

The output files are stored in a folder, which has the following name structure::

    <IM>-site_<num_site>-poe-<num_poe>

where:
    - `<IM>` is the required intensity measure
    - `<num_site>` is the site number
    - `<num_poe>` is the probability of exceedance number

"""

import sys
import os
from .read_input_data import read_input_data
from .scaling_module import scaling_module
from .check_module import check_module
from .selection_module import selection_module
from .scaling_nodes_NTC18 import scaling_nodes_NTC18

if __name__ == '__main__':

    # %% Initial setup
    try:
        fileini = sys.argv[1]
        calculation_mode = sys.argv[2]
    except IndexError:
        sys.exit('usage:\n'
                 'python3 -m haselrec #input_file [mode]' + "\n"
                 + '       [--run-complete]' + "\n"
                 + '       [--run-selection]' + "\n"
                 + '       [--run-scaling]' + "\n"
                 + '       [--check-rec]'
                 + '       [--scaling-nodes-NTC18] [TR]')
    if calculation_mode == '--scaling-nodes-NTC18':
        try: 
            TR = sys.argv[3]
            path_to_share = sys.argv[4]
        except IndexError:
            sys.exit('usage:\n'
                 'python3 -m haselrec #input_file [mode]' + "\n"
                 + '       [--scaling-nodes-NTC18] [TR]'
                 + 'path_to_share')

    # Read fileini

    [intensity_measures, site_code, rlz_code, path_results_classical,
     path_results_disagg, num_disagg, num_classical,
     probability_of_exceedance_num, probability_of_exceedance,
     investigation_time, target_periods, tstar, im_type,
     im_type_lbl, avg_periods, corr_type, gmpe_input, rake, vs30, vs30type,
     hypo_depth, dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
     database_path, allowed_database, allowed_recs_vs30, allowed_ec8_code,
     maxsf_input, radius_dist_input, dist_range_input, radius_mag_input, allowed_depth, n_gm,
     random_seed, n_trials, weights, n_loop, penalty, path_nga_folder,
     path_esm_folder, output_folder, meanMag_disagg, meanDist_disagg, 
     hazard_value, hazard_mode, component, correlated_motion,
     code_spectrum_file, period_range_spectrumcompatibility, 
     threshold_up, threshold_low, selection_type, path_kiknet_folder,
     radius_dist_type_input, radius_mag_type_input, vertical_component,
     tstar1,tstar2] = read_input_data(fileini)

    if calculation_mode == '--run-complete' or \
            calculation_mode == '--run-selection':

        selection_module(intensity_measures, site_code, rlz_code,
                         path_results_classical, path_results_disagg,
                         num_disagg, num_classical,
                         probability_of_exceedance_num,
                         probability_of_exceedance, investigation_time,
                         target_periods, tstar, im_type, im_type_lbl,
                         avg_periods, corr_type, gmpe_input, rake, vs30,
                         vs30type, hypo_depth, dip, azimuth, fhw, z2pt5, z1pt0,
                         upper_sd, lower_sd, database_path, allowed_database,
                         allowed_recs_vs30, allowed_ec8_code, maxsf_input,
                         radius_dist_input, dist_range_input, radius_mag_input, allowed_depth,
                         n_gm, random_seed, n_trials, weights, n_loop, penalty,
                         output_folder, meanMag_disagg, meanDist_disagg, 
                         hazard_value, hazard_mode, component, correlated_motion,
                         code_spectrum_file, period_range_spectrumcompatibility,
                         threshold_up, threshold_low, selection_type,
                         radius_dist_type_input, radius_mag_type_input,
                         tstar1,tstar2)
                        

    if calculation_mode == '--check-rec':
        check_module(output_folder, site_code, probability_of_exceedance_num,
                     intensity_measures, n_gm, path_nga_folder, path_kiknet_folder,
                     selection_type)

    if calculation_mode == '--scaling-nodes-NTC18':
        scaling_nodes_NTC18(TR, output_folder, n_gm,
                       path_nga_folder, path_esm_folder, path_kiknet_folder,
                       selection_type, path_to_share, component)

    if calculation_mode == '--run-complete' or \
            calculation_mode == '--run-scaling':
        if not os.path.exists(path_nga_folder):
            os.makedirs(path_nga_folder)
        if not os.path.exists(path_esm_folder):
            os.makedirs(path_esm_folder)
        if not os.path.exists(path_kiknet_folder):
            os.makedirs(path_kiknet_folder)

        scaling_module(site_code, probability_of_exceedance_num,
                       intensity_measures, output_folder, n_gm,
                       path_nga_folder, path_esm_folder, path_kiknet_folder,
                       selection_type, vertical_component)
