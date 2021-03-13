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

def optimize_ground_motion(n_loop, n_gm, sample_small, n_big, id_sel, ln_sa1,
                             maxsf, sample_big, tgt_per, mean_req, stdevs,
                             weights, penalty, rec_id, im_scale_fac):
    """
    Executes incremental changes to the initially selected ground motion set to
    further optimize its fit to the target spectrum distribution. From:
    Jayaram N, Lin T, Baker J. (2011). A Computationally Efficient Ground-Motion
    Selection Algorithm for Matching a Target Response Spectrum Mean and
    Variance. Earthq Spectra 2011;27:797-815. https://doi.org/10.1193/1.3608002.
    """

    import numpy as np

    sample_small = np.array(sample_small)
    print(
        'Please wait...This algorithm takes a few minutes '
        'depending on the number of records to be selected')

    for _ in np.arange(n_loop):
        # consider replacing each ground motion in the selected set
        for i in np.arange(n_gm):
            min_dev = 100000

            sample_small = np.delete(sample_small, i, 0)
            rec_id = np.delete(rec_id, i)
            scale_fac = np.ones(n_big)
            min_id = []

            # Try to add a new spectrum to the subset list
            for j in np.arange(n_big):
                rec_value = np.exp(
                    sum(sample_big[j, id_sel]) / len(id_sel))
                if rec_value == 0:
                    scale_fac[j] = 1000000
                else:
                    scale_fac[j] = np.exp(ln_sa1) / rec_value
                added1 = np.reshape(
                    (sample_big[j, :] + np.log(scale_fac[j])),
                    (1, len(tgt_per)))
                sample_small = np.concatenate(
                    (sample_small, added1))  # add candidate to set
                # Compute deviations from target
                dev_mean = np.mean(sample_small, axis=0) - mean_req
                dev_sig = np.std(sample_small, axis=0) - stdevs
                dev_total = weights[0] * sum(dev_mean ** 2) + weights[
                    1] * sum(dev_sig ** 2)

                # Penalize bad spectra
                # (set penalty to zero if this is not required)
                if penalty != 0:
                    for m in np.arange(len(sample_small)):
                        dev_total = dev_total + sum(
                            np.absolute(
                                np.exp(sample_small[m, :]) > np.exp(
                                    mean_req + 3 * stdevs))) \
                                   * penalty

                if scale_fac[j] > maxsf or scale_fac[j] < 1. / maxsf:
                    dev_total = dev_total + 1000000

                # Should cause improvement and record should not
                # be repeated
                if dev_total < min_dev and not any(rec_id == j):
                    min_id = np.zeros(1, dtype=int)
                    min_id[0] = j
                    min_dev = dev_total
                end = len(sample_small)
                sample_small = sample_small[0:end - 1, :]

            # Add new element in the right slot
            im_scale_fac[i] = scale_fac[min_id]
            end = len(sample_small)
            added2 = np.reshape(
                (sample_big[min_id, :] + np.log(scale_fac[min_id])),
                (1, len(tgt_per)))
            if i > 0:
                sample_small = np.concatenate((sample_small[0:i, :],
                                               added2,
                                               sample_small[i:end,
                                               :]))
                rec_id = np.concatenate(
                    (rec_id[0:i], min_id, rec_id[i:end]))
            else:
                sample_small = np.concatenate(
                    (added2, sample_small[i:end, :]))
                rec_id = np.concatenate((min_id, rec_id[i:end]))
    return rec_id, im_scale_fac, sample_small
