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

def find_ground_motion(tgt_per, tstar, avg_periods, intensity_measures, n_gm,
                       sa_known, ind_per, mean_req, n_big, simulated_spectra,
                       maxsf, event, station, allowed_index, correlated_motion,
                       selection_type, period_range, cluster, tstar1, tstar2):
    """
    Select ground motions from the database that individually match the
    statistically simulated spectra. From:
    Jayaram N, Lin T, Baker J. (2011) A Computationally Efficient Ground-Motion
    Selection Algorithm for Matching a Target Response Spectrum Mean and
    Variance. Earthq Spectra 2011;27:797-815. https://doi.org/10.1193/1.3608002.
    """
    import numpy as np
    import sys

    sample_big = np.log(sa_known[:, ind_per])

    id_sel = []
    id_spectrum_compatibility=[]
    if(selection_type=='conditional-spectrum'):
        if intensity_measures == 'AvgSA':
            id_sel_bool = np.isin(tgt_per, avg_periods)
            for i in np.arange(len(tgt_per)):
                if id_sel_bool[i]:
                    id_sel.append(i)
            id_sel = np.array(id_sel)
        else:
            id_sel = np.where(tgt_per == tstar)
        if(len(id_sel[0]) == 0):
            sys.exit('Error: tstar not included in tgt_per',tstar,tgt_per)
        ln_sa1 = np.mean(mean_req[id_sel])
        w = []
    elif(selection_type=='code-spectrum'):
        w=np.zeros(len(tgt_per))
        if tstar1>-1 and tstar2>-1:
            logical = np.logical_and(tgt_per >= tstar1, tgt_per <= tstar2)
            id_sel = np.where(logical == True)
        else:
            id_sel = np.where(tgt_per == tstar)
        if(len(id_sel[0]) == 0):
            sys.exit('Error: tstar not included in tgt_per',tstar,tgt_per)
        w[id_sel] = 1
        for itp in range(len(tgt_per)):
            if tgt_per[itp] >= period_range[0] and tgt_per[itp] <= period_range[1]:
                id_spectrum_compatibility.append(itp)
        ln_sa1 = []
    id_spectrum_compatibility = np.array(id_spectrum_compatibility)

    rec_id = np.zeros(n_gm, dtype=int)
    sample_small = []
    im_scale_fac = np.ones(n_gm)
    # Find database spectra most similar to each simulated spectrum
    for i in np.arange(n_gm):  # for each simulated spectrum
        err = np.zeros(n_big) * 1000000  # initialize error matrix
        scale_fac = np.ones(n_big)  # initialize scale factors to 1
        # compute scale factors and errors for each candidate
        # ground motion
        for j in np.arange(n_big):
            if(selection_type=='conditional-spectrum'):
                rec_value = np.exp(
                    sum(sample_big[j, id_sel]) / len(id_sel))
                # rec_value=rec_value[0]
                if rec_value == 0:
                    scale_fac[j] = 1000000
                else:
                    scale_fac[j] = np.exp(ln_sa1) / rec_value
                err[j] = sum(
                    (np.log(
                        np.exp(sample_big[j, :]) * scale_fac[j]) - np.log(
                        simulated_spectra[i, :])) ** 2)
            elif(selection_type=='code-spectrum'):
                rec_value = np.exp(
                    sum(sample_big[j, id_sel]) / len(id_sel))
                if rec_value.any() == 0:
                    scale_fac[j] = 1000000
                else:
                    scale_fac[j]=sum( w * np.log(np.exp(mean_req) / np.exp(sample_big[j, :]))) / sum(w) 
                    scale_fac[j]=np.exp(scale_fac[j])
                #Computed Mean Squared Error (MSE) of the selected record with respect to
                #the target spectrum
                err[j] = sum(w*
                    (mean_req-np.log(np.exp(sample_big[j, :]) * scale_fac[j])) ** 2)/sum(w)

        # exclude previously-selected ground motions
        if(i>0):
            err[rec_id[0:i - 1]] = 1000000
            if(correlated_motion=='no'):
            # exclude ground motions from the same stations and earthquake of 
            # previously-selected ground motions
                for l in range(i):
                    rec_idx_l = allowed_index[rec_id[l]]
                    for j in np.arange(n_big):
                        rec_idx_j = allowed_index[j]
                        if station[rec_idx_j]==station[rec_idx_l] or event[rec_idx_j]==event[rec_idx_l] or cluster[rec_idx_j]==cluster[rec_idx_l]:
                            err[j] = 1000000

        # exclude ground motions requiring too large SF
        err[scale_fac > maxsf] = 1000000
        # exclude ground motions requiring too large SF
        err[scale_fac < 1. / maxsf] = 1000000
        err[scale_fac < 1. / maxsf] = 1000000

        # find minimum-error ground motion
        rec_id[i] = np.argmin(err)
        min_err = np.min(err)
        assert (min_err < 1000), (
            'Warning: problem with target spectrum. '
            'No good matches found')
        im_scale_fac[i] = scale_fac[rec_id[i]]  # store scale factor
        sample_small.append(
            np.log(np.exp(sample_big[rec_id[i], :]) * scale_fac[
                rec_id[i]]))  # store scaled log spectrum

    return (sample_small, sample_big, id_sel, ln_sa1, rec_id,
            im_scale_fac, w, id_spectrum_compatibility)
