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
