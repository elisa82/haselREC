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


def compute_avgsa(avg_periods, sctx, rctx, dctx, bgmpe, corr_type):
    """
    """
    # Import libraries
    from openquake.hazardlib import imt, const
    import numpy as np
    from lib.im_correlation import baker_jayaram_correlation
    from lib.im_correlation import akkar_correlation

    mean_list = []
    stddvs_list = []
    # Loop over averaging periods
    for period in avg_periods:
        # compute mean and standard deviation
        p = imt.SA(period)
        s = [const.StdDev.TOTAL]
        mean, std = bgmpe().get_mean_and_stddevs(sctx, rctx, dctx, p, s)
        mean_list.append(mean)
        stddvs_list.append(std[0])  # Support only for total!

    mean_avgsa = 0.
    stddvs_avgsa = 0.

    for i1 in np.arange(len(avg_periods)):
        mean_avgsa += mean_list[i1]
        for i2 in np.arange(len(avg_periods)):
            rho = []
            if corr_type == 'baker_jayaram':
                rho = baker_jayaram_correlation(avg_periods[i1],
                                                avg_periods[i2])
            if corr_type == 'akkar':
                rho = akkar_correlation(avg_periods[i1], avg_periods[i2])

            stddvs_avgsa += rho * stddvs_list[i1] * stddvs_list[i2]

    mean_avgsa *= (1. / len(avg_periods))
    stddvs_avgsa *= (1. / len(avg_periods)) ** 2
    stddvs_avgsa = np.sqrt(stddvs_avgsa)
    return [np.exp(mean_avgsa), stddvs_avgsa]


def compute_rho_avgsa(per, avg_periods, sctx, rctx, dctx, stddvs_avgsa, bgmpe,
                      corr_type):
    """
    """
    # Import libraries
    from openquake.hazardlib import imt, const
    from lib.im_correlation import baker_jayaram_correlation
    from lib.im_correlation import akkar_correlation

    sum_numeratore = 0
    for i1 in avg_periods:
        rho = []
        if corr_type == 'baker_jayaram':
            rho = baker_jayaram_correlation(per, i1)
        if corr_type == 'akkar':
            rho = akkar_correlation(per, i1)
        s = [const.StdDev.TOTAL]
        mean1, std1 = bgmpe().get_mean_and_stddevs(sctx, rctx, dctx,
                                                   imt.SA(i1), s)
        sum_numeratore = sum_numeratore + rho * std1[0]

    denominatore = len(avg_periods) * stddvs_avgsa
    rho_avgsa = sum_numeratore / denominatore
    return rho_avgsa
