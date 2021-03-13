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

def compute_rho_avgsa(per, avg_periods, sctx, rctx, dctx, stddvs_avgsa, bgmpe,
                      corr_type):
    """
    """
    # Import libraries
    from openquake.hazardlib import imt, const, gsim
    from .modified_akkar_correlation_model import ModifiedAkkarCorrelationModel

    sum_numeratore = 0
    for i1 in avg_periods:
        rho = []
        if corr_type == 'baker_jayaram':
            rho = gsim.mgmpe.generic_gmpe_avgsa. \
                    BakerJayaramCorrelationModel([per, i1])(0, 1)
        if corr_type == 'akkar':
            rho = ModifiedAkkarCorrelationModel([per, i1])(0, 1)
        s = [const.StdDev.TOTAL]
        mean1, std1 = bgmpe().get_mean_and_stddevs(sctx, rctx, dctx,
                                                   imt.SA(i1), s)
        sum_numeratore = sum_numeratore + rho * std1[0]

    denominatore = len(avg_periods) * stddvs_avgsa
    rho_avgsa = sum_numeratore / denominatore
    return rho_avgsa
