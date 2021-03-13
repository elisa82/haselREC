# Original Copyright (C) 2012-2021 GEM Foundation
# Modified by Andrea Francia and Elisa Zuccolo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# this code was adapted and combined from these openquake files:
#
# - https://github.com/gem/oq-engine/blob/master/openquake/hazardlib/gsim/mgmpe/akkar_coeff_table.py
# - https://github.com/gem/oq-engine/blob/master/openquake/hazardlib/gsim/mgmpe/generic_gmpe_avgsa.py
#
# we can't use the AkkarCorrelationModel from openquake directly because it does
# not consider the zero period

import pathlib

from openquake.hazardlib.gsim.mgmpe.generic_gmpe_avgsa import BaseAvgSACorrelationModel
from scipy.interpolate import interp1d

this = pathlib.Path(__file__)

coeff_table = []
with open(this.parent / 'modified_akkar_coeff_table.csv') as f:
    for row in f:
        coeff_table.append([float(col) for col in row.split(',')])

akkar_periods = [0.0,
           0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.11, 0.12, 0.13, 0.14,
           0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3,
           0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.55, 0.6,
           0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
           1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]
import numpy as np

class ModifiedAkkarCorrelationModel(BaseAvgSACorrelationModel):
    """
    Read the period-dependent correlation coefficient matrix as in:
    Akkar S., Sandikkaya MA., Ay BO., 2014, Compatible ground-motion
    prediction equations for damping scaling factors and vertical to
    horizontal spectral amplitude ratios for the broader Europe region,
    Bull Earthquake Eng, 12, pp. 517-547.
    """
    def build_correlation_matrix(self):
        """
        Constructs the correlation matrix by two-step linear interpolation
        from the correlation table
        """
        irho = np.array(coeff_table)
        iper = np.array(akkar_periods)
        if np.any(self.avg_periods < iper[0]) or\
                np.any(self.avg_periods > iper[-1]):
            raise ValueError("'avg_periods' contains values outside of the "
                             "range supported by the Akkar et al. (2014) "
                             "correlation model")
        ipl1 = interp1d(iper, irho, axis=1)
        ipl2 = interp1d(iper, ipl1(self.avg_periods), axis=0)
        self.rho = ipl2(self.avg_periods)

    def get_correlation(self, t1, t2):
        """
        Computes the correlation coefficient for the specified periods.

        :param float t1:
            First period of interest.

        :param float t2:
            Second period of interest.

        :return float:
            The predicted correlation coefficient.
        """
        periods = np.array(akkar_periods)
        rho = np.array(coeff_table)
        if t1 < periods[0] or t1 > periods[-1]:
            raise ValueError("t1 %.3f is out of valid period range (%.3f to "
                             "%.3f" % (t1, periods[0], periods[-1]))

        if t2 < periods[0] or t2 > periods[-1]:
            raise ValueError("t1 %.3f is out of valid period range (%.3f to "
                             "%.3f" % (t2, periods[0], periods[-1]))
        iloc1 = np.searchsorted(periods, t1)
        iloc2 = np.searchsorted(periods, t2)
        if iloc1:
            rho1 = rho[iloc1 - 1, :] + (t1 - periods[iloc1 - 1]) *\
                ((periods[iloc1] - periods[iloc1 - 1]) /
                 (rho[iloc1, :] - rho[iloc1 - 1, :]))
        else:
            rho1 = rho[0, :]
        if iloc2:
            rho2 = rho1[iloc2 - 1] + (t2 - periods[iloc2 - 1]) *\
                ((periods[iloc2] - periods[iloc2 - 1]) /
                 (rho1[iloc2] - rho1[iloc2 - 1]))
        else:
            rho2 = rho1[0]
        return rho2

