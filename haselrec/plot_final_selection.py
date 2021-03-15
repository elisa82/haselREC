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

def plot_final_selection(name, lbl, n_gm, t_cs, sample_small, mean_req, stdevs,
                         output_folder):
    """
    Three plots are generated::

        1) <IM>-site_<num_site>-poe-<num_poe>_spectra.pdf
        2) <IM>-site_<num_site>-poe-<num_poe>_spectra_gms.pdf
        3) <IM>-site_<num_site>-poe-<num_poe>_dispersion.pdf

    where:
        - `<IM>` is the required intensity measure
        - `<num_site>` is the site number
        - `<num_poe>` is the probability of exceedance number

    The plots represent:

        1) the response spectra of selected ground motions (:code:`nGM`
           green lines) and their distribution (black lines), and the target
           conditional spectrum (red lines). Solid lines: average spectrum,
           dashed lines: average spectrum values plus and minus 2 standard
           deviations.

        2) the distribution of the response spectra of selected ground
           motions (black lines), and the target conditional spectrum
           (red lines). Solid lines: average spectrum, dashed lines: average
           spectrum values plus and minus 2 standard deviations.

        3) the dispersion of the response spectra of selected ground
           motions (black line) and the dispersion of the target conditional
           spectrum (red lines).
    """
    # Import libraries
    import numpy as np
    import matplotlib.pyplot as plt

    meanrecorded = np.mean(np.exp(sample_small), axis=0)
    meanrecorded_p2sigma = np.percentile(np.exp(sample_small), 50 + 34.1 + 13.6,
                                         axis=0)
    meanrecorded_n2sigma = np.percentile(np.exp(sample_small), 50 - 34.1 - 13.6,
                                         axis=0)
    meanrecorded_eps = (np.log(meanrecorded_p2sigma) - np.log(
        meanrecorded_n2sigma)) / (2 * 1.96)
    indexes = [i for i, x in enumerate(t_cs) if x > 0]
    stdevs_gt0 = []
    t_cs_gt0 = []
    mean_req_gt0 = []
    meanrecorded_eps_gt0 = []
    meanrecorded_p2sigma_gt0 = []
    meanrecorded_n2sigma_gt0 = []
    meanrecorded_gt0 = []
    for i in indexes:
        t_cs_gt0.append(t_cs[i])
        meanrecorded_eps_gt0.append(meanrecorded_eps[i])
        meanrecorded_p2sigma_gt0.append(meanrecorded_p2sigma[i])
        meanrecorded_n2sigma_gt0.append(meanrecorded_n2sigma[i])
        meanrecorded_gt0.append(meanrecorded[i])
        stdevs_gt0.append(stdevs[i])
        mean_req_gt0.append(mean_req[i])
    stdevs_gt0 = np.asarray(stdevs_gt0)
    t_cs_gt0 = np.asarray(t_cs_gt0)
    mean_req_gt0 = np.asarray(mean_req_gt0)
    meanrecorded_eps_gt0 = np.asarray(meanrecorded_eps_gt0)
    meanrecorded_p2sigma_gt0 = np.asarray(meanrecorded_p2sigma_gt0)
    meanrecorded_n2sigma_gt0 = np.asarray(meanrecorded_n2sigma_gt0)
    meanrecorded_gt0 = np.asarray(meanrecorded_gt0)

    # Spectra with ground motions
    plt.figure(figsize=(1.5 * 2.36, 2.36))
    plt.rcParams.update({'font.size': 8})
    for i in np.arange(n_gm):
        sample_small_gt0 = []
        for j in indexes:
            sample_small_gt0.append(sample_small[i, j])
        plt.loglog(t_cs_gt0, np.exp(sample_small_gt0), 'g', linewidth=.5)
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0), 'r', label='CMS', linewidth=1.0)
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0 + 2 * stdevs_gt0), '--r',
               label=r'CMS $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0 - 2 * stdevs_gt0), '--r',
               linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_gt0, 'k', label='Selected', linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_p2sigma_gt0, '--k',
               label=r'Selected $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_n2sigma_gt0, '--k', linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel('Acceleration [g]')
    plt.xlim(min(t_cs_gt0), max(t_cs_gt0))
    plt.ylim(1e-2, 1e1)
    plt.yscale('log')
    plt.xscale('log')
    #number=int(name[11])+1 #per AvgSA
    #number = int(name[9]) + 1 #per PGA
    #plt.title('site '+str(number)+' - '+lbl)
    plt.grid(True)
    plt.legend()
    #plt.savefig(output_folder + '/' + name + '/' + name + '_spectra_gms.png',
    #            bbox_inches='tight')
    plt.savefig(output_folder + '/' + name + '/' + name + '_spectra_gms.pdf',
                bbox_inches='tight')
    plt.close()

    # Spectra
    plt.figure(figsize=(1.5 * 2.36, 2.36))
    plt.rcParams.update({'font.size': 8})
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0), 'r', label='CMS', linewidth=1.0)
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0 + 2 * stdevs_gt0), '--r',
               label=r'CMS $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(t_cs_gt0, np.exp(mean_req_gt0 - 2 * stdevs_gt0), '--r',
               linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_gt0, 'k', label='Selected', linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_p2sigma_gt0, '--k',
               label=r'Selected $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(t_cs_gt0, meanrecorded_n2sigma_gt0, '--k', linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel('Acceleration [g]')
    plt.xlim(min(t_cs_gt0), max(t_cs_gt0))
    plt.ylim(1e-2, 1e1)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_folder + '/' + name + '/' + name + '_spectra.pdf',
                bbox_inches='tight')
    plt.close()

    # Dispersion
    plt.figure(figsize=(1.5 * 2.36, 2.36))
    plt.rcParams.update({'font.size': 8})
    plt.plot(t_cs_gt0, stdevs_gt0, 'r', label='CMS', linewidth=1.0)
    plt.plot(t_cs_gt0, meanrecorded_eps_gt0, 'k', label='Selected',
             linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel('Dispersion')
    plt.xlim(min(t_cs_gt0), max(t_cs_gt0))
    plt.ylim(0, 1)
    plt.xscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_folder + '/' + name + '/' + name + '_dispersion.pdf',
                bbox_inches='tight')
    plt.close()
