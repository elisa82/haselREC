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

def selection_module(intensity_measures, site_code, rlz_code,
                     path_results_classical, path_results_disagg, num_disagg,
                     num_classical, probability_of_exceedance_num,
                     probability_of_exceedance, investigation_time,
                     target_periods, tstar, im_type, im_type_lbl, avg_periods,
                     corr_type, gmpe_input, rake, vs30, vs30type, hypo_depth,
                     dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
                     database_path, allowed_database, allowed_recs_vs30,
                     allowed_ec8_code, maxsf_input, radius_dist_input,
                     radius_mag_input, allowed_depth, n_gm, random_seed,
                     n_trials, weights, n_loop, penalty, output_folder):
    import os
    import pandas as pd
    import numpy as np
    from lib.screen_database import screen_database
    from lib.simulate_spectra import simulate_spectra
    from lib.plot_final_selection import plot_final_selection
    from lib.input_GMPE import inizialize_gmm
    from lib.create_output_files import create_output_files
    from lib.compute_cs import compute_cs
    from lib.find_ground_motion import find_ground_motion
    from lib.optimizing_ground_motion import optimizing_ground_motion

    # %% Start the routine
    print('Inputs loaded, starting selection....')
    ind = 1

    # For each site investigated
    for ii in np.arange(len(site_code)):

        # Get the current site and realisation indices
        site = site_code[ii]
        rlz = rlz_code[ii]

        # For each hazard of poe level investigated
        for jj in np.arange(len(probability_of_exceedance_num)):

            poe = probability_of_exceedance_num[jj]

            if hasattr(maxsf_input, '__len__'):
                maxsf = maxsf_input[jj]
            else:
                maxsf = maxsf_input
            if hasattr(radius_dist_input, '__len__'):
                radius_dist = radius_dist_input[jj]
            else:
                radius_dist = radius_dist_input
            if hasattr(radius_mag_input, '__len__'):
                radius_mag = radius_mag_input[jj]
            else:
                radius_mag = radius_mag_input

            # For each intensity measure investigated
            for im in np.arange(len(intensity_measures)):

                # Get the name of the disaggregation file to look in
                disagg_results = 'rlz-' + str(rlz) + '-' + intensity_measures[
                    im] + '-sid-' + str(site) + '-poe-' + str(
                    poe) + '_Mag_Dist_' + str(num_disagg) + '.csv'
                name = intensity_measures[im] + '-site_' + str(
                    site) + '-poe-' + str(poe)
                selected_column = intensity_measures[im] + '-' + str(
                    probability_of_exceedance[jj])
                file_with_oq_acc_value = 'hazard_map-mean_' + str(
                    num_classical) + '.csv'

                # Print some on screen feedback
                print('Processing ' + name + ' Case: ' + str(ind) + '/' + str(
                    len(site_code) * len(probability_of_exceedance_num) * len(
                        intensity_measures)))
                ind += 1

                # Retrieve disaggregation results
                df = pd.read_csv(
                    ''.join([path_results_disagg, '/', disagg_results]),
                    skiprows=1)
                df['rate'] = -np.log(1 - df['poe']) / investigation_time
                df['rate_norm'] = df['rate'] / df['rate'].sum()
                # mode = df.sort_values(by='rate_norm', ascending=False)[0:1]
                mean_mag = np.sum(df['mag'] * df['rate_norm'])
                mean_dist = np.sum(df['dist'] * df['rate_norm'])
                allowed_recs_d = [mean_dist - radius_dist, mean_dist +
                                  radius_dist]
                allowed_recs_mag = [mean_mag - radius_mag,
                                    mean_mag + radius_mag]

                # Retrieve conditioning value
                df = pd.read_csv(''.join(
                    [path_results_classical, '/', file_with_oq_acc_value]),
                    skiprows=1)
                output_oq = df[selected_column]
                im_star = output_oq[site]

                # -----------------------------------------------------------------------------

                rjb = np.array([mean_dist])
                mag = mean_mag

                [bgmpe, sctx, rctx, dctx, vs30, rrup] = \
                    inizialize_gmm(ii, gmpe_input, rjb, mag, hypo_depth, dip,
                                   rake, upper_sd, lower_sd, azimuth, fhw,
                                   vs30type, vs30, z2pt5, z1pt0)

                # Screen the database of available ground motions

                [sa_known, ind_per, tgt_per, n_big, allowed_index, event_id,
                 station_code, source, record_sequence_number_nga, event_mw,
                 event_mag, acc_distance, station_vs30, station_ec8] = \
                    screen_database(database_path, allowed_database,
                                    allowed_recs_vs30, allowed_recs_mag,
                                    allowed_recs_d, allowed_ec8_code,
                                    target_periods, n_gm, allowed_depth, vs30)

                # Compute the target spectrum

                [mean_req, cov_req, stdevs] = \
                    compute_cs(tgt_per, bgmpe, sctx, rctx, dctx, im_type[im],
                               tstar[im], rrup, mag, avg_periods, corr_type,
                               im_star)

                simulated_spectra = simulate_spectra(random_seed,
                                                     n_trials,
                                                     mean_req,
                                                     cov_req,
                                                     stdevs,
                                                     n_gm,
                                                     weights)

                [sample_small, sample_big, id_sel, ln_sa1,
                 rec_id, im_scale_fac] = \
                    find_ground_motion(tgt_per, tstar[im], avg_periods,
                                       intensity_measures[im], n_gm,
                                       sa_known, ind_per, mean_req,
                                       n_big, simulated_spectra, maxsf)

                # Further optimize the ground motion selection

                [final_records, final_scale_factors, sample_small] = \
                    optimizing_ground_motion(n_loop, n_gm, sample_small, n_big,
                                             id_sel, ln_sa1, maxsf, sample_big,
                                             tgt_per, mean_req, stdevs, weights,
                                             penalty, rec_id,
                                             im_scale_fac)

                # Create the outputs folder
                folder = output_folder + '/' + name
                if not os.path.exists(folder):
                    os.makedirs(folder)

                # Plot the figure
                plot_final_selection(name, im_type_lbl[im], n_gm, tgt_per,
                                     sample_small, mean_req, stdevs,
                                     output_folder)

                # Collect information of the final record set
                rec_idx = [allowed_index[i] for i in final_records]
                # Create the summary file along with the file with the CS
                create_output_files(output_folder, name, im_star, mean_mag,
                                    mean_dist, n_gm, rec_idx, source, event_id,
                                    station_code, event_mw, acc_distance,
                                    station_vs30, station_ec8,
                                    final_scale_factors, tgt_per, mean_req,
                                    stdevs, record_sequence_number_nga,
                                    event_mag)

    return
