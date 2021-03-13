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

def compute_cs(t_cs, bgmpe, sctx, rctx, dctx, im_type, t_star, rrup, mag,
               avg_periods, corr_type, im_star, gmpe_input):
    """
    Compute the conditional spectrum according to the procedure outlined
    in Baker JW, Lee C. An Improved Algorithm for Selecting Ground Motions
    to Match a Conditional Spectrum. J Earthq Eng 2018;22:708-23.
    https://doi.org/10.1080/13632469.2016.1264334.

    When the IM and the GMM are defined for the maximum of the two horizontal
    components, the The Boore and Kishida (2017) relationship is applied to
    convert the maximum of the two horizontal components into `RotD50`. This is
    done only for `PGA` and `SA`.
    """
    import numpy as np
    import sys
    from openquake.hazardlib import imt, const, gsim
    from .compute_avgSA import compute_rho_avgsa

    # Use the same periods as the available spectra to construct the
    # conditional spectrum

    p = []
    s = [const.StdDev.TOTAL]
    if im_type == 'AvgSA':
        _ = gsim.get_available_gsims()
        p = imt.AvgSA()
        mgmpe = gsim.mgmpe.generic_gmpe_avgsa.GenericGmpeAvgSA \
            (gmpe_name=gmpe_input, avg_periods=avg_periods, corr_func=corr_type)
        mu_im_cond, sigma_im_cond = mgmpe.get_mean_and_stddevs(sctx, rctx, dctx,
                                                               p, s)
    else:
        if im_type == 'PGA':
            p = imt.PGA()
        else:
            p = imt.SA(t_star)
        s = [const.StdDev.TOTAL]
        mu_im_cond, sigma_im_cond = bgmpe().get_mean_and_stddevs(
            sctx, rctx, dctx, p, s)
    sigma_im_cond = sigma_im_cond[0]

    if (bgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT ==
            'Greater of two horizontal'):
        if im_type == 'PGA' or im_type == 'SA':
            from shakelib.conversions.imc.boore_kishida_2017 import \
                BooreKishida2017

            bk17 = BooreKishida2017(const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                    const.IMC.RotD50)
            mu_im_cond = bk17.convertAmps(p, mu_im_cond, rrup,
                                          float(mag))
            sigma_im_cond = bk17.convertSigmas(p, sigma_im_cond[0])
        else:
            sys.exit('Error: conversion between intensity measures is not '
                     'possible for AvgSA')

    # Compute how many standard deviations the PSHA differs from
    # the GMPE value
    epsilon = (np.log(im_star) - mu_im_cond) / sigma_im_cond

    mu_im = np.zeros(len(t_cs))
    sigma_im = np.zeros(len(t_cs))
    rho_t_tstar = np.zeros(len(t_cs))
    mu_im_im_cond = np.zeros(len(t_cs))

    for i in range(len(t_cs)):
        # Get the GMPE ouput for a rupture scenario
        if t_cs[i] == 0.:
            p = imt.PGA()
        else:
            p = imt.SA(t_cs[i])
        s = [const.StdDev.TOTAL]
        mu0, sigma0 = bgmpe().get_mean_and_stddevs(sctx, rctx, dctx, p, s)

        if (bgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT ==
                'Greater of two horizontal'):
            if im_type == 'PGA' or im_type == 'SA':
                from shakelib.conversions.imc.boore_kishida_2017 \
                    import BooreKishida2017

                bk17 = BooreKishida2017(const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                        const.IMC.RotD50)

                mu0 = bk17.convertAmps(p, mu0, rrup, float(mag))
                sigma0 = bk17.convertSigmas(p, sigma0[0])

        mu_im[i] = mu0[0]
        sigma_im[i] = sigma0[0][0]
        rho = None
        if im_type == 'AvgSA':
            rho = compute_rho_avgsa(t_cs[i], avg_periods, sctx,
                                    rctx, dctx, sigma_im_cond,
                                    bgmpe, corr_type)
            rho = rho[0]

        else:
            if corr_type == 'baker_jayaram':
                rho = gsim.mgmpe.generic_gmpe_avgsa. \
                    BakerJayaramCorrelationModel([t_cs[i], t_star])(0, 1)
            if corr_type == 'akkar':
                rho = gsim.mgmpe.generic_gmpe_avgsa. \
                    AkkarCorrelationModel([t_cs[i], t_star])(0, 1)
        rho_t_tstar[i] = rho
        # Get the value of the CMS
        mu_im_im_cond[i] = \
            mu_im[i] + rho_t_tstar[i] * epsilon[0] * sigma_im[i]

    # Compute covariances and correlations at all periods
    cov = np.zeros((len(t_cs), len(t_cs)))
    for i in np.arange(len(t_cs)):
        for j in np.arange(len(t_cs)):
            var1 = sigma_im[i] ** 2
            var2 = sigma_im[j] ** 2
            var_tstar = sigma_im_cond ** 2

            sigma_corr = []
            if corr_type == 'baker_jayaram':
                sigma_corr = gsim.mgmpe.generic_gmpe_avgsa. \
                    BakerJayaramCorrelationModel([t_cs[i], t_cs[j]])(0, 1) * \
                    np.sqrt(var1 * var2)
            if corr_type == 'akkar':
                sigma_corr = gsim.mgmpe.generic_gmpe_avgsa. \
                    AkkarCorrelationModel([t_cs[i], t_cs[j]])(0, 1) * \
                    np.sqrt(var1 * var2)
            sigma11 = np.matrix(
                [[var1, sigma_corr], [sigma_corr, var2]])
            sigma22 = np.array(var_tstar)
            sigma12 = np.array(
                [rho_t_tstar[i] * np.sqrt(var1 * var_tstar),
                 rho_t_tstar[j] * np.sqrt(var_tstar * var2)])
            sigma_cond = sigma11 - sigma12 * 1. / (
                sigma22) * sigma12.T
            cov[i, j] = sigma_cond[0, 1]

    # find covariance values of zero and set them to a small number
    # so that random number generation can be performed
    cov[np.absolute(cov) < 1e-10] = 1e-10
    stdevs = np.sqrt(np.diagonal(cov))

    return mu_im_im_cond, cov, stdevs
