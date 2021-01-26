#!/usr/bin/env python

# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# OpenSel is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenSel is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenSel. If not, see <http://www.gnu.org/licenses/>.

# %% Import libraries
# Standard built-in libraries

import sys
import os
from lib.read_input_data import read_input_data
from lib.scaling_module import scaling_module
from lib.check_module import check_module
from lib.selection_module import selection_module

# %% Initial setup
try:
    fileini = sys.argv[1]
    calculation_mode = sys.argv[2]
except IndexError:
    sys.exit('usage: opensel JOB.INI [option]' + "\n"
             + '       [--run-complete]' + "\n"
             + '       [--run-selection]' + "\n"
             + '       [--run-scaling]' + "\n"
             + '       [--check-NGArec]')
# fileini = 'demo/job_selection.ini '

# Read fileini

[intensity_measures, site_code, rlz_code, path_results_classical,
 path_results_disagg, num_disagg, num_classical, probability_of_exceedance_num,
 probability_of_exceedance, investigation_time, target_periods, tstar, im_type,
 im_type_lbl, avg_periods, corr_type, gmpe_input, rake, vs30, vs30type,
 hypo_depth, dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd, database_path,
 allowed_database, allowed_recs_vs30, allowed_ec8_code, maxsf_input,
 radius_dist_input, radius_mag_input, allowed_depth, n_gm, random_seed,
 n_trials, weights, n_loop, penalty, path_nga_folder, path_esm_folder,
 output_folder] = read_input_data(fileini)

if calculation_mode == '--run-complete' or \
        calculation_mode == '--run-selection':

    selection_module(intensity_measures, site_code, rlz_code,
                     path_results_classical, path_results_disagg, num_disagg,
                     num_classical, probability_of_exceedance_num,
                     probability_of_exceedance, investigation_time,
                     target_periods, tstar, im_type, im_type_lbl, avg_periods,
                     corr_type, gmpe_input, rake, vs30, vs30type, hypo_depth,
                     dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
                     database_path, allowed_database, allowed_recs_vs30,
                     allowed_ec8_code, maxsf_input, radius_dist_input,
                     radius_mag_input, allowed_depth, n_gm, random_seed,
                     n_trials, weights, n_loop, penalty, output_folder)

if calculation_mode == '--check-NGArec':
    check_module(output_folder, site_code, probability_of_exceedance_num,
                 intensity_measures, n_gm, path_nga_folder)
    command='tail -n +2 '+output_folder+'/missing_NGArec.txt  | sort | uniq > '+output_folder+'/temp'
    os.system(command)
    command='mv '+output_folder+'/temp '+output_folder+'/missing_NGArec.txt'
    os.system(command)

if calculation_mode == '--run-complete' or calculation_mode == '--run-scaling':
    scaling_module(site_code, probability_of_exceedance_num,
                   intensity_measures, output_folder, n_gm,
                   path_nga_folder, path_esm_folder)

