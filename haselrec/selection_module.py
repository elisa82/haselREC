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

def selection_module(intensity_measures, site_code, rlz_code,
                     path_results_classical, path_results_disagg, num_disagg,
                     num_classical, probability_of_exceedance_num,
                     probability_of_exceedance, investigation_time,
                     target_periods, tstar, im_type, im_type_lbl, avg_periods,
                     corr_type, gmpe_input, rake, vs30_input, vs30type, hypo_depth,
                     dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
                     database_path, allowed_database, allowed_recs_vs30,
                     allowed_ec8_code, maxsf_input, radius_dist_input, dist_range_input,
                     radius_mag_input, allowed_depth, n_gm, random_seed,
                     n_trials, weights, n_loop, penalty, output_folder,
                     meanMag_disagg, meanDist_disagg, 
                     hazard_value, hazard_mode, component, correlated_motion,
                     code_spectrum_file, period_range_spectrumcompatibility, 
                     threshold_up, threshold_low, selection_type, 
                     radius_dist_type_input, radius_mag_type_input,tstar1,tstar2,
                     lon,lat,path_to_AllegatoB):
    """
    This module is called when mode :code:`--run-selection` is specified.

    It performs record selection following these steps:

        1) retrieve conditioning value (:code:`compute_conditioning_value` module)
        2) defines all inputs necessary to apply ground motion prediction equations
           (:code:`inizialize_gmm` module)
        3) screening of the database of candidate ground motion
           (:code:`screen_database` module)
        4) computation of the target response spectrum distribution
           (:code:`compute_cs` module)
        5) statistical simulation of response spectra from the target
           distribution
           (:code:`simulate_spectra` module)
        6) selection of ground motions from the database that individually match
           the statistically simulated spectra
           (:code:`find_ground_motion` module)
        7) execution of incremental changes to the initially selected ground
           motion set to further optimize its fit to the target spectrum
           distribution (:code:`optimize_ground_motion` module)
        8) produce output files (3 figures created by :code:`plot_final_selection`
           module and 2 `.txt` files created by :code:`create_output_files` modules)

    """
    import os
    import numpy as np
    import pandas as pd
    from .compute_conditioning_value import compute_conditioning_value
    from .screen_database import screen_database
    from .simulate_spectra import simulate_spectra
    from .plot_final_selection import plot_final_selection
    from .input_GMPE import inizialize_gmm
    from .create_output_files import create_output_files
    from .compute_cs import compute_cs
    from .find_ground_motion import find_ground_motion
    from .optimize_ground_motion import optimize_ground_motion
    from .create_NTC18 import create_NTC18

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
            if isinstance(radius_dist_type_input, str):
                radius_dist_type = radius_dist_type_input
            else:
                radius_dist_type = radius_dist_type_input[jj]
            if isinstance(radius_mag_type_input, str):
                radius_mag_type = radius_mag_type_input
            else:
                radius_mag_type = radius_mag_type_input[jj]

            # For each intensity measure investigated
            for im in np.arange(len(intensity_measures)):

                if(selection_type=='conditional-spectrum'):

                    name = intensity_measures[im] + '-site_' + str(
                        site) + '-poe-' + str(poe)

                    # Print some on screen feedback
                    print('Processing ' + name + ' Case: ' + str(ind) + '/' + str(
                        len(site_code) * len(probability_of_exceedance_num) * len(
                            intensity_measures)))
                    ind += 1

                    [im_star, rjb, mag] = \
                        compute_conditioning_value(rlz, intensity_measures[im],
                                                   site, poe, num_disagg,
                                                   probability_of_exceedance[jj],
                                                   num_classical,
                                                   path_results_disagg,
                                                   investigation_time,
                                                   path_results_classical,
                                                   meanMag_disagg, meanDist_disagg, 
                                                   hazard_value, hazard_mode, ii )

                    [bgmpe, sctx, rctx, dctx, vs30, rrup] = \
                        inizialize_gmm(ii, gmpe_input, rjb, mag, hypo_depth, dip,
                                       rake, upper_sd, lower_sd, azimuth, fhw,
                                       vs30type, vs30_input, z2pt5, z1pt0, site_code)

                    # Screen the database of available ground motions
                    if(bgmpe()=='[Ambraseys1996]'):
                        mag=np.exp(1.421+0.108*mag)-1.863
                    else:
                        mag=mag

                elif(selection_type=='code-spectrum'):
                    TR=int(probability_of_exceedance[jj])
                    name = 'site_' + str(site) + '-RT-' + str(TR) + '_years'

                    # Print some on screen feedback
                    print('Processing ' + name + ' Case: ' + str(ind) + '/' + str(
                        len(site_code) * len(probability_of_exceedance_num)))
                    ind += 1
                    if code_spectrum_file[ii]=='NTC18':
                        target_periods=np.arange(0.,4.01,0.01)
                        code_spectrum = create_NTC18(lon[ii],lat[ii],TR,target_periods,path_to_AllegatoB)
                    else:                    
                        code = pd.read_csv(code_spectrum_file[ii])
                        target_periods = code['T'].to_numpy()
                        code_spectrum = code['Sa'].to_numpy()

                    mag = meanMag_disagg[ii]
                    rjb = np.array([meanDist_disagg[ii]])

                    vs30=[]

                [sa_known, ind_per, tgt_per, n_big, allowed_index, event_id,
                 station_code, source, record_sequence_number_nga, event_mw,
                 event_mag, acc_distance, station_vs30, station_ec8, comp_allowed,
                 fminNS2, fmaxNS2, fminEW2, fmaxEW2, cluster] = \
                    screen_database(database_path, allowed_database,
                                    allowed_recs_vs30, radius_dist, dist_range_input,
                                    radius_mag, rjb[0], mag, allowed_ec8_code,
                                    target_periods, n_gm, allowed_depth, vs30, 
                                    component, radius_dist_type, radius_mag_type,
                                    selection_type)

                # Compute the target spectrum

                if(selection_type=='conditional-spectrum'):

                    [mean_req, cov_req, stdevs] = \
                        compute_cs(tgt_per, bgmpe, sctx, rctx, dctx, im_type[im],
                                   tstar[im], rrup, mag, avg_periods, corr_type,
                                   im_star, gmpe_input)

                    simulated_spectra = simulate_spectra(random_seed,
                                                         n_trials,
                                                         mean_req,
                                                         cov_req,
                                                         stdevs,
                                                         n_gm,
                                                         weights)
                elif(selection_type=='code-spectrum'):

                    ind_per_code = []
                    for itp in np.arange(len(tgt_per)):
                        ind_per_code.append(np.argmin(np.abs(target_periods - tgt_per[itp])))
                    ind_per_code=np.asarray(ind_per_code)
                    sa_ref=code_spectrum[ind_per_code]
                    if component == 'two-component':
                        sa_ref=code_spectrum[ind_per_code]*1.3
                    mean_req=np.log(sa_ref)
                    simulated_spectra=[]
                    stdevs=np.zeros(len(mean_req))
                    im_star=0

                [sample_small, sample_big, id_sel, ln_sa1,
                 rec_id, im_scale_fac, w, id_spectrum_compatibility] = \
                    find_ground_motion(tgt_per, tstar[im], avg_periods,
                                       intensity_measures[im], n_gm,
                                       sa_known, ind_per, mean_req,
                                       n_big, simulated_spectra, maxsf,
                                       event_id, station_code, allowed_index,
                                       correlated_motion,selection_type,
                                       period_range_spectrumcompatibility,
                                       cluster,tstar1,tstar2)
                                      

                # Further optimize the ground motion selection

                [final_records, final_scale_factors, sample_small] = \
                    optimize_ground_motion(n_loop, n_gm, sample_small, n_big,
                                             id_sel, ln_sa1, maxsf, sample_big,
                                             tgt_per, mean_req, stdevs, weights,
                                             penalty, rec_id, im_scale_fac, 
                                             event_id, station_code, allowed_index,
                                             correlated_motion, selection_type,
                                             period_range_spectrumcompatibility,
                                             threshold_up, threshold_low, w, 
                                             id_spectrum_compatibility, cluster)
                                             
                # Create the outputs folder
                folder = output_folder + '/' + name
                if not os.path.exists(folder):
                    os.makedirs(folder)

                # Plot the figure
                plot_final_selection(name, im_type_lbl[im], n_gm, tgt_per,
                                     sample_small, mean_req, stdevs,
                                     folder, selection_type,
                                     period_range_spectrumcompatibility,
                                     threshold_up, threshold_low)

                # Collect information of the final record set
                rec_idx = [allowed_index[i] for i in final_records]
                comp_idx = [comp_allowed[i] for i in final_records]

                min_misfit=[]
                max_misfit=[]
                average_misfit=[]
                if(selection_type=='code-spectrum'):
                    misfit=((np.mean(np.exp(sample_small[:,id_spectrum_compatibility]),axis=0)-np.exp(mean_req[id_spectrum_compatibility]))/np.exp(mean_req[id_spectrum_compatibility]))*100
                    average_misfit=np.mean(np.abs(misfit))
                    if np.min(misfit)>0:
                        min_misfit=0
                    else:
                        min_misfit=np.max(np.abs(misfit[misfit<0]))
                    if np.max(misfit) > 0:
                        max_misfit=np.max(misfit[misfit>0])
                    else:
                        max_misfit = 0


                # Create the summary file along with the file with the CS
                create_output_files(folder, name, im_star, mag,
                                    rjb[0], n_gm, rec_idx, source, event_id,
                                    station_code, event_mw, acc_distance,
                                    station_vs30, station_ec8,
                                    final_scale_factors, tgt_per, mean_req,
                                    stdevs, record_sequence_number_nga,
                                    event_mag, comp_idx, 
                                    fminNS2, fmaxNS2, fminEW2, fmaxEW2,
                                    selection_type,min_misfit, max_misfit, 
                                    average_misfit, sample_small)
                
    return 
