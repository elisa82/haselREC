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

def simulate_spectra(random_seed, n_trials, mean_req, cov_req, stdevs, n_gm,
                     weights):
    """
    Statistically simulates response spectra from the target distribution. From:
    Jayaram N, Lin T, Baker J. (2011). A Computationally Efficient Ground-Motion
    Selection Algorithm for Matching a Target Response Spectrum Mean and
    Variance. Earthq Spectra 2011;27:797-815. https://doi.org/10.1193/1.3608002.
    """
    # Import libraries
    import numpy as np
    from scipy.stats import skew

    random = np.random.RandomState(random_seed)
    # simulate response spectra from the target mean and covariance matrix
    dev_total_sim = []
    spettri = []
    for _ in np.arange(n_trials):
        spectra_sample = np.exp(random.multivariate_normal(mean_req, cov_req,
                                                           n_gm))
        spettri.append(spectra_sample)
        # evaluate simulation
        sample_mean_err = np.mean(np.log(spectra_sample), axis=0) - mean_req
        sample_std_err = np.std(np.log(spectra_sample), axis=0) - stdevs
        sample_skewness_err = skew(np.log(spectra_sample), axis=0, bias=True)
        dev_total_sim.append(
            weights[0] * sum(sample_mean_err ** 2) + weights[1] **
            sum(sample_std_err ** 2) + weights[2] *
            sum(sample_skewness_err ** 2))

    # find the simulated spectra that best match the target
    best_sample = dev_total_sim.index(min(dev_total_sim))
    # return the best set of simulations
    simulated_spectra = spettri[best_sample]
    return simulated_spectra
