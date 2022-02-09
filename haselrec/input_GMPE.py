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

def inizialize_gmm(index, gmpe_input, rjb, mag, z_hyp_input, dip_input, rake,
                   upper_sd_input, lower_sd_input, azimuth_input, fhw, vs30type,
                   vs30_input, z2pt5_input, z1pt0_input, site_code):
    """
    Defines all the parameters for the computation of a Ground Motion Model. If
    not defined by the user as input parameters, most parameters (dip,
    hypocentral depth, fault width, ztor, azimuth, source-to-site distances
    based on extended sources, z2pt5, z1pt0) are defined  according to the
    relationships included in:
    Kaklamanos J, Baise LG, Boore DM. (2011) Estimating unknown input parameters
    when implementing the NGA ground-motion prediction equations in engineering
    practice. Earthquake Spectra 27: 1219-1235.
    https://doi.org/10.1193/1.3650372.
    """

    import sys
    from openquake.hazardlib import gsim
    import numpy as np

    bgmpe = None
    for name_gmpe, gmpes in gsim.get_available_gsims().items():
        if name_gmpe == gmpe_input:
            bgmpe = gmpes
    if bgmpe is None:
        sys.exit('The GMM is not found')

    sctx = gsim.base.SitesContext()
    rctx = gsim.base.RuptureContext()
    dctx = gsim.base.DistancesContext()

    # -----------------------------------------------------------------------------
    # Initialise contexts

    dip, z_hyp, width, ztor, azimuth = compute_source_params(mag, z_hyp_input,
                                                             dip_input, rake,
                                                             upper_sd_input,
                                                             lower_sd_input,
                                                             azimuth_input, fhw)

    [rx, rrup, ry] = compute_dists(rjb, mag, z_hyp_input, dip_input, rake,
                                   upper_sd_input,
                                   lower_sd_input, azimuth_input, fhw)

    [vs30, vs30measured, z1pt0, z2pt5] = compute_soil_params(vs30_input,
                                                             z2pt5_input,
                                                             z1pt0_input,
                                                             gmpe_input,
                                                             vs30type,
                                                             index)

    setattr(rctx, 'width', width)
    setattr(rctx, 'ztor', ztor)
    setattr(rctx, 'dip', dip)
    setattr(dctx, 'rx', rx)
    setattr(dctx, 'rrup', rrup)
    setattr(dctx, 'ry0', ry)
    z1pt0 = z1pt0 + np.zeros(rjb.shape)
    setattr(sctx, 'z1pt0', z1pt0)
    z2pt5 = z2pt5 + np.zeros(rjb.shape)
    setattr(sctx, 'z2pt5', z2pt5)
    setattr(sctx, 'vs30measured', vs30measured)
    setattr(rctx, 'mag', mag)
    setattr(rctx, 'hypo_depth', z_hyp)
    setattr(rctx, 'rake', rake)
    setattr(rctx, 'occurrence_rate', 0.)
    setattr(dctx, 'rjb', rjb)
    vs30 = vs30 + np.zeros(rjb.shape)
    setattr(sctx, 'vs30', vs30)
    sc = site_code[index] + np.zeros(rjb.shape)
    setattr(sctx, 'sids', sc)

    return bgmpe, sctx, rctx, dctx, vs30, rrup


def compute_source_params(mag, z_hyp_input, dip_input, rake, upper_sd_input,
                          lower_sd_input, azimuth_input, fhw):
    import numpy as np

    if z_hyp_input is None:
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            z_hyp = 5.63 + 0.68 * mag
        else:
            z_hyp = 11.24 - 0.2 * mag
    else:
        z_hyp = z_hyp_input

    if dip_input is None:
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            dip = 90
        elif rake > 0:
            dip = 40
        else:
            dip = 50
    else:
        dip = dip_input

    if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
        # strike slip
        width = 10.0 ** (-0.76 + 0.27 * mag)
    elif rake > 0:
        # thrust/reverse
        width = 10.0 ** (-1.61 + 0.41 * mag)
    else:
        # normal
        width = 10.0 ** (-1.14 + 0.35 * mag)

    if upper_sd_input is None:
        upper_sd = 0
    else:
        upper_sd = upper_sd_input

    if lower_sd_input is None:
        lower_sd = 500
    else:
        lower_sd = lower_sd_input

    source_vertical_width = width * np.sin(np.radians(dip))
    ztor = max(z_hyp - 0.6 * source_vertical_width, upper_sd)
    if (ztor + source_vertical_width) > lower_sd:
        source_vertical_width = lower_sd - ztor
        width = source_vertical_width / np.sin(np.radians(dip))

    azimuth = []
    if azimuth_input is None:
        if fhw == 1:
            azimuth = 50
        elif fhw == -1:
            azimuth = -50
    else:
        azimuth = azimuth_input

    return dip, z_hyp, width, ztor, azimuth


def compute_dists(rjb, mag, z_hyp_input, dip_input, rake, upper_sd_input,
                  lower_sd_input, azimuth_input, fhw):
    """
    """
    import numpy as np

    dip, z_hyp, width, ztor, azimuth = compute_source_params(mag, z_hyp_input,
                                                             dip_input, rake,
                                                             upper_sd_input,
                                                             lower_sd_input,
                                                             azimuth_input, fhw)

    if rjb == 0:
        rx = 0.5 * width * np.cos(np.radians(dip))
    else:
        if dip == 90:
            rx = rjb * np.sin(np.radians(azimuth))
        else:
            if (0 <= azimuth < 90) or (
                    90 < azimuth <= 180):
                if (rjb * np.abs(np.tan(np.radians(azimuth))) <= width * np.cos(
                        np.radians(dip))):
                    rx = rjb * np.abs(np.tan(np.radians(azimuth)))
                else:
                    rx = rjb * np.tan(np.radians(azimuth)) * np.cos(
                        np.radians(azimuth) - np.arcsin(
                            width * np.cos(np.radians(dip)) * np.cos(
                                np.radians(azimuth)) / rjb))
            elif azimuth == 90:  # we assume that Rjb>0
                rx = rjb + width * np.cos(np.radians(dip))
            else:
                rx = rjb * np.sin(np.radians(azimuth))

    if azimuth == 90 or azimuth == -90:
        ry = 0
    elif azimuth == 0 or azimuth == 180 or azimuth == -180:
        ry = rjb
    else:
        ry = np.abs(rx * 1. / np.tan(np.radians(azimuth)))

    if dip == 90:
        rrup = np.sqrt(np.square(rjb) + np.square(ztor))
    else:
        rrup1 = []
        if rx < ztor * np.tan(np.radians(dip)):
            rrup1 = np.sqrt(np.square(rx) + np.square(ztor))
        if (ztor * np.tan(np.radians(dip)) <= rx <= ztor * np.tan(
                np.radians(dip)) + width * 1. / np.cos(np.radians(dip))):
            rrup1 = rx * np.sin(np.radians(dip)) + ztor * np.cos(
                np.radians(dip))
        if (rx > ztor * np.tan(np.radians(dip)) + width * 1. / np.cos(
                np.radians(dip))):
            rrup1 = np.sqrt(
                np.square(rx - width * np.cos(np.radians(dip))) + np.square(
                    ztor + width * np.sin(np.radians(dip))))
        rrup = np.sqrt(np.square(rrup1) + np.square(ry))
    return rx, rrup, ry


def compute_soil_params(vs30_input, z2pt5_input, z1pt0_input, gmpe_input,
                        vs30type, index):
    import numpy as np

    vs30 = float(vs30_input[index])

    vs30measured = True
    if vs30type[index] == 'inferred':
        vs30measured = False

    z1pt0 = []
    if z1pt0_input is None:
        if gmpe_input == 'AbrahamsonEtAl2014' or \
                ((gmpe_input == 'CampbellBozorgnia2008' or
                  gmpe_input == 'CampbellBozorgnia2014')
                 and z2pt5_input is None):
            if vs30 < 180:
                z1pt0 = np.exp(6.745)
            elif 180 <= vs30 <= 500:
                z1pt0 = np.exp(6.745 - 1.35 * np.log(vs30 / 180))
            else:
                z1pt0 = np.exp(5.394 - 4.48 * np.log(vs30 / 500))

        elif gmpe_input == 'ChiouYoungs2014':
            z1pt0 = np.exp(
                28.5 - 3.82 / 8 * np.log(vs30 ** 8 + 378.7 ** 8))
    else:
        z1pt0 = float(z1pt0_input[index])

    z2pt5 = []
    if gmpe_input == 'CampbellBozorgnia2008' or \
            gmpe_input == 'CampbellBozorgnia2014':
        if z2pt5_input is None:
            z2pt5 = 519 + 3.595 * z1pt0
        else:
            z2pt5 = float(z2pt5_input[index])
    return vs30, vs30measured, z1pt0, z2pt5
