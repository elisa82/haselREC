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
                       maxsf):
    """
    Select ground motions from the database that individually match the
    statistically simulated spectra. From:
    Jayaram N, Lin T, Baker J. (2011) A Computationally Efficient Ground-Motion
    Selection Algorithm for Matching a Target Response Spectrum Mean and
    Variance. Earthq Spectra 2011;27:797-815. https://doi.org/10.1193/1.3608002.
    """
    import numpy as np

    sample_big = np.log(sa_known[:, ind_per])

    id_sel = []
    if intensity_measures == 'AvgSA':
        id_sel_bool = np.isin(tgt_per, avg_periods)
        for i in np.arange(len(tgt_per)):
            if id_sel_bool[i]:
                id_sel.append(i)
        id_sel = np.array(id_sel)
    else:
        id_sel = np.where(tgt_per == tstar)
    ln_sa1 = np.mean(mean_req[id_sel])

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

        # exclude previously-selected ground motions
        err[rec_id[0:i - 1]] = 1000000
        # exclude ground motions requiring too large SF
        err[scale_fac > maxsf] = 1000000
        # exclude ground motions requiring too large SF
        err[scale_fac < 1. / maxsf] = 1000000

        # find minimum-error ground motion
        rec_id[i] = np.argmin(err)
        min_err = np.min(err)
        assert (min_err < 1000), (
            'Warning: problem with simulated spectrum. '
            'No good matches found')
        im_scale_fac[i] = scale_fac[rec_id[i]]  # store scale factor
        sample_small.append(
            np.log(np.exp(sample_big[rec_id[i], :]) * scale_fac[
                rec_id[i]]))  # store scaled log spectrum

    return (sample_small, sample_big, id_sel, ln_sa1, rec_id,
            im_scale_fac)
