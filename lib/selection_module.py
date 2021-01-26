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
    from scipy.stats import skew
    from openquake.hazardlib import imt, const
    from lib.im_correlation import akkar_correlation
    from lib.im_correlation import baker_jayaram_correlation
    from lib.screen_database import screen_database
    from lib.compute_avgSA import compute_avgsa
    from lib.compute_avgSA import compute_rho_avgsa
    from lib.simulate_spectra import simulate_spectra
    from lib.plot_final_selection import plot_final_selection
    from lib.input_GMPE import inizialize_gmm

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
                file_with_OQ_acc_value = 'hazard_map-mean_' + str(
                    num_classical) + '.csv'

                # Print some on screen feedback
                print('Processing ' + name + ' Case: ' + str(ind) + '/' + str(
                    len(site_code) * len(probability_of_exceedance_num) * len(
                        intensity_measures)))
                ind += 1

                # Retrieve disaggregation results
                meanLst = [], []
                df = pd.read_csv(
                    ''.join([path_results_disagg, '/', disagg_results]),
                    skiprows=1)
                df['rate'] = -np.log(1 - df['poe']) / investigation_time
                df['rate_norm'] = df['rate'] / df['rate'].sum()
                mode = df.sort_values(by='rate_norm', ascending=False)[0:1]
                meanMag = np.sum(df['mag'] * df['rate_norm'])
                meanDist = np.sum(df['dist'] * df['rate_norm'])
                allowed_recs_d = [meanDist - radius_dist, meanDist +
                                  radius_dist]
                allowed_recs_mag = [meanMag - radius_mag, meanMag + radius_mag]

                # Retrieve conditioning value
                df = pd.read_csv(''.join(
                    [path_results_classical, '/', file_with_OQ_acc_value]),
                    skiprows=1)
                output_oq = df[selected_column]
                im_star = output_oq[site]

                # -----------------------------------------------------------------------------

                rjb = np.array([meanDist])
                mag = meanMag

                [bgmpe, sctx, rctx, dctx, vs30, rrup] = \
                    inizialize_gmm(ii, gmpe_input, rjb, mag, hypo_depth, dip,
                                   rake, upper_sd, lower_sd, azimuth, fhw,
                                   vs30type, vs30, z2pt5, z1pt0)

                # Screen the database of available ground motions

                [SaKnown, indPer, TgtPer, nBig, allowedIndex, event_id,
                 station_code, source,
                 record_sequence_number_NGA, event_mw, event_mag,
                 acc_distance, station_vs30,
                 station_ec8] = screen_database(database_path, allowed_database,
                                                allowed_recs_vs30,
                                                allowed_recs_mag,
                                                allowed_recs_d,
                                                allowed_ec8_code,
                                                target_periods, n_gm,
                                                allowed_depth, vs30)

                # Get the GMPE ouput

                P = None
                if im_type[im] == 'AvgSA':
                    [median_im_cond, sigma_im_cond] = compute_avgsa(avg_periods,
                                                                    sctx,
                                                                    rctx,
                                                                    dctx,
                                                                    bgmpe,
                                                                    corr_type)
                    mu_im_cond = np.log(median_im_cond)
                else:
                    if im_type[im] == 'PGA':
                        P = imt.PGA()
                    elif im_type[im] == 'SA':
                        P = imt.SA(tstar[im])
                    S = [const.StdDev.TOTAL]
                    mu_im_cond, sigma_im_cond = bgmpe().get_mean_and_stddevs(
                        sctx, rctx, dctx, P, S)
                    sigma_im_cond = sigma_im_cond[0]

                # Compute how many standard deviations the PSHA differs from
                # the GMPE value
                epsilon = (np.log(im_star) - mu_im_cond) / sigma_im_cond

                if (bgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT ==
                        'Greater of two horizontal' and
                        (im_type[im] == 'PGA' or im_type[im] == 'SA')):
                    from shakelib.conversions.imc.boore_kishida_2017 import \
                        BooreKishida2017

                    bk17 = BooreKishida2017(const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                            const.IMC.RotD50)
                    mu_im_cond = bk17.convertAmps(P, mu_im_cond, rrup,
                                                  float(mag))
                    sigma_im_cond = bk17.convertSigmas(P, sigma_im_cond[0])

                # Use the same periods as the available spectra to construct the
                # conditional spectrum
                T_CS = TgtPer
                mu_im = np.zeros(len(T_CS))
                sigma_im = np.zeros(len(T_CS))
                rho_T_Tstar = np.zeros(len(T_CS))
                mu_im_im_cond = np.zeros(len(T_CS))

                for i in range(len(T_CS)):
                    # Get the GMPE ouput for a rupture scenario
                    if T_CS[i] == 0.:
                        P = imt.PGA()
                    else:
                        P = imt.SA(T_CS[i])
                    S = [const.StdDev.TOTAL]
                    mu0, sigma0 = bgmpe().get_mean_and_stddevs(sctx, rctx, dctx,
                                                               P, S)

                    if (bgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT ==
                            'Greater of two horizontal' and
                            (im_type[im] == 'PGA' or im_type[im] == 'SA')):
                        from shakelib.conversions.imc.boore_kishida_2017 \
                            import BooreKishida2017

                        bk17 = BooreKishida2017(
                            const.IMC.GREATER_OF_TWO_HORIZONTAL,
                            const.IMC.RotD50)

                        mu0 = bk17.convertAmps(P, mu0, rrup, float(mag))
                        sigma0 = bk17.convertSigmas(P, sigma0[0])

                    mu_im[i] = mu0[0]
                    sigma_im[i] = sigma0[0][0]
                    rho = None
                    if im_type[im] == 'AvgSA':
                        rho = compute_rho_avgsa(T_CS[i], avg_periods, sctx,
                                                rctx, dctx, sigma_im_cond,
                                                bgmpe, corr_type)
                        rho = rho[0]
                    elif im_type[im] == 'SA' or im_type[im] == 'PGA':
                        if corr_type == 'baker_jayaram':
                            rho = baker_jayaram_correlation(T_CS[i], tstar[im])
                        if corr_type == 'akkar':
                            rho = akkar_correlation(T_CS[i], tstar[im])
                    rho_T_Tstar[i] = rho
                    # Get the value of the CMS
                    mu_im_im_cond[i] = \
                        mu_im[i] + rho_T_Tstar[i] * epsilon[0] * sigma_im[i]

                # Compute covariances and correlations at all periods

                # Compute the Covariance
                Cov = np.zeros((len(T_CS), len(T_CS)))
                for i in np.arange(len(T_CS)):
                    for j in np.arange(len(T_CS)):
                        var1 = sigma_im[i] ** 2
                        var2 = sigma_im[j] ** 2
                        varTstar = sigma_im_cond ** 2

                        sigma_Corr = []
                        if corr_type == 'baker_jayaram':
                            sigma_Corr = \
                                baker_jayaram_correlation(T_CS[i], T_CS[j]) * \
                                np.sqrt(var1 * var2)
                        if corr_type == 'akkar':
                            sigma_Corr = akkar_correlation(T_CS[i],
                                                           T_CS[j]) * np.sqrt(
                                var1 * var2)
                        sigma11 = np.matrix(
                            [[var1, sigma_Corr], [sigma_Corr, var2]])
                        sigma22 = np.array(varTstar)
                        sigma12 = np.array(
                            [rho_T_Tstar[i] * np.sqrt(var1 * varTstar),
                             rho_T_Tstar[j] * np.sqrt(varTstar * var2)])
                        sigma_cond = sigma11 - sigma12 * 1. / (
                            sigma22) * sigma12.T
                        Cov[i, j] = sigma_cond[0, 1]

                # find covariance values of zero and set them to a small number
                # so that random number generation can be performed
                Cov[np.absolute(Cov) < 1e-10] = 1e-10

                meanReq = mu_im_im_cond
                covReq = Cov
                stdevs = np.sqrt(np.diagonal(covReq))

                random = np.random.RandomState(random_seed)
                simulated_spectra = simulate_spectra(random,
                                                     n_trials,
                                                     meanReq,
                                                     covReq,
                                                     stdevs,
                                                     n_gm,
                                                     weights)

                sampleBig = np.log(SaKnown[:, indPer])

                id_sel = []
                if intensity_measures[im] == 'AvgSA':
                    id_sel_bool = np.isin(TgtPer, avg_periods)
                    for i in np.arange(len(TgtPer)):
                        if id_sel_bool[i]:
                            id_sel.append(i)
                    id_sel = np.array(id_sel)
                else:
                    id_sel = np.where(TgtPer == tstar[im])
                lnSa1 = np.mean(meanReq[id_sel])

                recID = np.zeros(n_gm, dtype=int)
                sampleSmall = []
                IMScaleFac = np.ones(n_gm)
                # Find database spectra most similar to each simulated spectrum
                for i in np.arange(n_gm):  # for each simulated spectrum
                    err = np.zeros(nBig) * 1000000  # initialize error matrix
                    scaleFac = np.ones(nBig)  # initialize scale factors to 1
                    # compute scale factors and errors for each candidate
                    # ground motion
                    for j in np.arange(nBig):
                        rec_value = np.exp(
                            sum(sampleBig[j, id_sel]) / len(id_sel))
                        # rec_value=rec_value[0]
                        if rec_value == 0:
                            scaleFac[j] = 1000000
                        else:
                            scaleFac[j] = np.exp(lnSa1) / rec_value
                        err[j] = sum(
                            (np.log(
                                np.exp(sampleBig[j, :]) * scaleFac[j]) - np.log(
                                simulated_spectra[i, :])) ** 2)

                    # exclude previously-selected ground motions
                    err[recID[0:i - 1]] = 1000000
                    # exclude ground motions requiring too large SF
                    err[scaleFac > maxsf] = 1000000
                    # exclude ground motions requiring too large SF
                    err[scaleFac < 1. / maxsf] = 1000000

                    # find minimum-error ground motion
                    recID[i] = np.argmin(err)
                    min_err = np.min(err)
                    assert (min_err < 1000), (
                        'Warning: problem with simulated spectrum. '
                        'No good matches found')
                    IMScaleFac[i] = scaleFac[recID[i]]  # store scale factor
                    sampleSmall.append(
                        np.log(np.exp(sampleBig[recID[i], :]) * scaleFac[
                            recID[i]]))  # store scaled log spectrum

                sampleSmall = np.array(sampleSmall)

                # Further optimize the ground motion selection

                print(
                    'Please wait...This algorithm takes a few minutes '
                    'depending on the number of records to be selected')

                for k in np.arange(n_loop):
                    # consider replacing each ground motion in the selected set
                    for i in np.arange(n_gm):
                        minDev = 100000

                        sampleSmall = np.delete(sampleSmall, i, 0)
                        recID = np.delete(recID, i)

                        # Try to add a new spectrum to the subset list
                        for j in np.arange(nBig):
                            rec_value = np.exp(
                                sum(sampleBig[j, id_sel]) / len(id_sel))
                            if rec_value == 0:
                                scaleFac[j] = 1000000
                            else:
                                scaleFac[j] = np.exp(lnSa1) / rec_value
                            added1 = np.reshape(
                                (sampleBig[j, :] + np.log(scaleFac[j])),
                                (1, len(TgtPer)))
                            sampleSmall = np.concatenate(
                                (sampleSmall, added1))  # add candidate to set
                            # Compute deviations from target
                            devMean = np.mean(sampleSmall, axis=0) - meanReq
                            devSkew = skew(sampleSmall, axis=0, bias=True)
                            devSig = np.std(sampleSmall, axis=0) - stdevs
                            devTotal = weights[0] * sum(devMean ** 2) + weights[
                                1] * sum(devSig ** 2)

                            # Penalize bad spectra
                            # (set penalty to zero if this is not required)
                            if penalty != 0:
                                for m in np.arange(len(sampleSmall)):
                                    devTotal = devTotal + sum(
                                        np.absolute(
                                            np.exp(sampleSmall[m, :]) > np.exp(
                                                meanReq + 3 * stdevs))) \
                                               * penalty

                            if scaleFac[j] > maxsf or scaleFac[j] < 1. / maxsf:
                                devTotal = devTotal + 1000000

                            # Should cause improvement and record should not
                            # be repeated
                            if devTotal < minDev and not any(recID == j):
                                minID = np.zeros(1, dtype=int)
                                minID[0] = j
                                minDev = devTotal
                            end = len(sampleSmall)
                            sampleSmall = sampleSmall[0:end - 1, :]

                        # Add new element in the right slot
                        IMScaleFac[i] = scaleFac[minID]
                        end = len(sampleSmall)
                        added2 = np.reshape(
                            (sampleBig[minID, :] + np.log(scaleFac[minID])),
                            (1, len(TgtPer)))
                        if i > 0:
                            sampleSmall = np.concatenate((sampleSmall[0:i, :],
                                                          added2,
                                                          sampleSmall[i:end,
                                                          :]))
                            recID = np.concatenate(
                                (recID[0:i], minID, recID[i:end]))
                        else:
                            sampleSmall = np.concatenate(
                                (added2, sampleSmall[i:end, :]))
                            recID = np.concatenate((minID, recID[i:end]))

                # Collect information of the final record set
                finalRecords = recID
                recIdx = [allowedIndex[i] for i in finalRecords]
                finalScaleFactors = IMScaleFac
                meanrecorded = np.mean(np.exp(sampleSmall), axis=0)
                meanrecorded_p2sigma = np.percentile(np.exp(sampleSmall),
                                                     50 + 34.1 + 13.6, axis=0)
                meanrecorded_n2sigma = np.percentile(np.exp(sampleSmall),
                                                     50 - 34.1 - 13.6, axis=0)
                meanrecorded_eps = (np.log(meanrecorded_p2sigma) - np.log(
                    meanrecorded_n2sigma)) / (2 * 1.96)

                # Create the outputs folder
                folder = output_folder + '/' + name
                if not os.path.exists(folder):
                    os.makedirs(folder)

                # Plot the figure
                plot_final_selection(name, im_type_lbl[im], n_gm, T_CS,
                                     sampleSmall, meanReq, stdevs,
                                     output_folder)

                # Output results to a text file
                blank = '-'
                name_summary = (output_folder + '/' + name + '/' + name +
                                "_summary_selection.txt")
                with open(name_summary, "w") as f:
                    f.write(
                        "{} {}\n".format('reference hazard value = ', im_star))
                    f.write("{} {}\n".format('mean_mag_disag = ', meanMag))
                    f.write("{} {}\n".format('mean_dist_disag = ', meanDist))
                    f.write(
                        "num source event_id_ESM station_code_ESM recID_NGA "
                        "magnitude distance vs30 EC8 scale_factor\n")
                    for i in np.arange(n_gm):
                        elemento = recIdx[i]
                        if source[elemento] == 'ESM':
                            f.write(
                                "{} {} {} {} {} {} {} {} {} {:4.2f}\n".format(
                                    i + 1, source[elemento],
                                    event_id[elemento],
                                    station_code[elemento], blank,
                                    event_mw[elemento],
                                    acc_distance[elemento],
                                    station_vs30[elemento],
                                    station_ec8[elemento],
                                    finalScaleFactors[i]))
                        if source[elemento] == 'NGA-West2':
                            val = int(record_sequence_number_NGA[elemento])
                            f.write(
                                "{} {} {} {} {} {} {} {} {} {:4.2f}\n".format(
                                    i + 1, source[elemento], blank, blank,
                                    val, event_mag[elemento],
                                    acc_distance[elemento],
                                    station_vs30[elemento],
                                    station_ec8[elemento],
                                    finalScaleFactors[i]))
                f.close()

                # Output conditional spectrum to a text file
                name_cs = output_folder + '/' + name + '/' + name + "_CS.txt"
                with open(name_cs, "w") as f:
                    f.write("Period(s) lnCS(g) standard_deviation\n")
                    for i in np.arange(len(TgtPer)):
                        f.write("{:6.2f}{:6.2f}{:6.2f} \n".format(TgtPer[i],
                                                                  meanReq[i],
                                                                  stdevs[i]))
                f.close()
    return
