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

def create_output_files(output_folder, name, im_star, mean_mag, mean_dist, n_gm,
                        rec_idx, source, event_id, station_code, event_mw,
                        acc_distance, station_vs30, station_ec8,
                        final_scale_factors, tgt_per, mean_req, stdevs,
                        record_sequence_number_nga, event_mag):
    """
    Two `.txt` files are generated::

        1) <IM>-site_<num_site>-poe-<num_poe>_CS.txt
        2) <IM>-site_<num_site>-poe-<num_poe>_summary_selection.txt

    where:
        - `<IM>` is the required intensity measure
        - `<num_site>` is the site number
        - `<num_poe>` is the probability of exceedance number

    The files contain:

        1) the CS. It has 3 columns: period (s), ln(CS) (g), standard deviation
           Example::

            Period(s) lnCS(g) standard_deviation
                 0.00  -0.53                0.00
                 0.01  -0.52                0.02
                 0.10  -0.08                0.24
                 0.20   0.18                0.31
                 0.30   0.11                0.39
                 0.40   0.03                0.46
                 0.50  -0.12                0.52

        2) A summary about the selection. It contains 3 rows with information
           about the conditioning value and the mean magnitude and distance from
           the disaggregation analysis + :code:`nGM` rows (one for each record)
           with the following information: sequential number used to identify
           recordings from the selection, source database, event ID and station
           code (for ESM recordings) or recording ID (for NGA-West2 recordings),
           magnitude of the earthquake, source-to-station distance,
           `vs30` of the station, EC8 soil category, applied scale factor.

           Example::

            reference hazard value =  0.5879783000000001
            mean_mag_disag =  6.347683361014376
            mean_dist_disag =  15.909873748773544
            num source event_id_ESM station_code_ESM recID_NGA magnitude distance vs30 EC8 scale_factor
            1 NGA-West2 - - 170 6.53 29.07 192.05 nan 2.62
            2 NGA-West2 - - 171 6.53 19.44 264.57 nan 1.88
            3 NGA-West2 - - 179 6.53 27.13 208.91 nan 1.54
            4 NGA-West2 - - 183 6.53 28.09 206.08 nan 1.12
            5 NGA-West2 - - 159 6.53 2.62 242.05 nan 2.46
            6 NGA-West2 - - 180 6.53 27.8 205.63 nan 1.43
            7 NGA-West2 - - 184 6.53 27.23 202.26 nan 1.35
            8 NGA-West2 - - 181 6.53 27.47 203.22 nan 1.31
            9 NGA-West2 - - 266 6.33 36.67 242.05 nan 4.86
            10 NGA-West2 - - 165 6.53 18.88 242.05 nan 2.22

    """
    import numpy as np

    # Output results to a text file
    blank = '-'
    name_summary = (output_folder + '/' + name + '/' + name +
                    "_summary_selection.txt")
    with open(name_summary, "w") as f:
        f.write(
            "{} {}\n".format('reference hazard value = ', im_star))
        f.write("{} {}\n".format('mean_mag_disag = ', mean_mag))
        f.write("{} {}\n".format('mean_dist_disag = ', mean_dist))
        f.write(
            "num source event_id_ESM station_code_ESM recID_NGA "
            "magnitude distance vs30 EC8 scale_factor\n")
        for i in np.arange(n_gm):
            elemento = rec_idx[i]
            if source[elemento] == 'ESM':
                f.write(
                    "{} {} {} {} {} {} {} {} {} {:4.2f}\n".format(
                        i + 1, source[elemento],
                        event_id[elemento],
                        station_code[elemento], blank,
                        event_mw[elemento],
                        acc_distance[elemento],
                        station_vs30[elemento],
                        station_ec8[elemento],
                        final_scale_factors[i]))
            if source[elemento] == 'NGA-West2':
                val = int(record_sequence_number_nga[elemento])
                f.write(
                    "{} {} {} {} {} {} {} {} {} {:4.2f}\n".format(
                        i + 1, source[elemento], blank, blank,
                        val, event_mag[elemento],
                        acc_distance[elemento],
                        station_vs30[elemento],
                        station_ec8[elemento],
                        final_scale_factors[i]))

    # Output conditional spectrum to a text file
    name_cs = output_folder + '/' + name + '/' + name + "_CS.txt"
    with open(name_cs, "w") as f:
        f.write("Period(s) lnCS(g) standard_deviation\n")
        for i in np.arange(len(tgt_per)):
            f.write("{:6.2f}{:6.2f}{:6.2f} \n".format(tgt_per[i],
                                                      mean_req[i],
                                                      stdevs[i]))
    return
