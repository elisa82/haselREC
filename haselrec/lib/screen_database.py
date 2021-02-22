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

def screen_database(database_path, allowed_database, allowed_recs_vs30,
                    allowed_recs_mag, allowed_recs_d, allowed_ec8_code,
                    target_periods, n_gm, allowed_depth, vs30):
    """
    Screen the database of candidate ground motion to select only appropriate
    ground motions
    """
    # Import libraries
    import numpy as np
    import pandas as pd

    known_per = np.array(
        [0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25,
         0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9,
         1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8,
         9, 10])

    dbacc = pd.read_csv(database_path, sep=';', engine='python')

    event_id = dbacc['event_id']
    event_mw = dbacc['Mw']
    event_mag = dbacc['M']
    record_sequence_number_nga = dbacc['record_sequence_number_NGA']
    station_ec8 = dbacc['ec8_code']
    station_vs30 = dbacc['vs30_m_sec']
    acc_distance = dbacc['epi_dist']
    station_code = dbacc['station_code']
    event_depth = dbacc['ev_depth_km']
    # sensor_depth = dbacc['sensor_depth_m']
    is_free_field_esm = dbacc['proximity_code']
    is_free_field_nga = dbacc['GMX_first']
    source = dbacc['source']
    epi_lon = dbacc['epi_lon']
    # epi_lat = dbacc['epi_lat']

    if allowed_recs_vs30 is None:
        if vs30 >= 800.0:
            allowed_recs_vs30 = [800.0, 3000.0]
        elif 360. <= vs30 < 800.:
            allowed_recs_vs30 = [360.0, 800.0]
        elif 180. <= vs30 < 360.:
            allowed_recs_vs30 = [180.0, 360.0]
        else:
            allowed_recs_vs30 = [0.0, 180.0]

    if allowed_ec8_code is None:
        if vs30 >= 800.0:
            allowed_ec8_code = 'A'
        elif 360. <= vs30 < 800.:
            allowed_ec8_code = 'B'
        elif 180. <= vs30 < 360.:
            allowed_ec8_code = 'C'
        else:
            allowed_ec8_code = 'D'

    # Match periods (known periods and target periods for error computations) 
    # save the indices of the matched periods in known_per
    ind_per = np.zeros((len(target_periods), 1), dtype=int)
    for i in np.arange(len(target_periods)):
        ind_per[i] = np.argmin(np.abs(known_per - target_periods[i]))

    # Remove any repeated values from TgtPer and redefine TgtPer as periods 
    # provided in databases
    ind_per = np.unique(ind_per)
    rec_per = known_per[ind_per]

    sa_list = []
    allowed_index = []
    for i in np.arange(len(event_id)):
        rotd50 = np.array([dbacc['rotD50_pga'][i], dbacc['rotD50_T0_010'][i],
                           dbacc['rotD50_T0_025'][i], dbacc['rotD50_T0_040'][i],
                           dbacc['rotD50_T0_050'][i], dbacc['rotD50_T0_070'][i],
                           dbacc['rotD50_T0_100'][i], dbacc['rotD50_T0_150'][i],
                           dbacc['rotD50_T0_200'][i], dbacc['rotD50_T0_250'][i],
                           dbacc['rotD50_T0_300'][i], dbacc['rotD50_T0_350'][i],
                           dbacc['rotD50_T0_400'][i], dbacc['rotD50_T0_450'][i],
                           dbacc['rotD50_T0_500'][i], dbacc['rotD50_T0_600'][i],
                           dbacc['rotD50_T0_700'][i], dbacc['rotD50_T0_750'][i],
                           dbacc['rotD50_T0_800'][i], dbacc['rotD50_T0_900'][i],
                           dbacc['rotD50_T1_000'][i], dbacc['rotD50_T1_200'][i],
                           dbacc['rotD50_T1_400'][i], dbacc['rotD50_T1_600'][i],
                           dbacc['rotD50_T1_800'][i], dbacc['rotD50_T2_000'][i],
                           dbacc['rotD50_T2_500'][i], dbacc['rotD50_T3_000'][i],
                           dbacc['rotD50_T3_500'][i], dbacc['rotD50_T4_000'][i],
                           dbacc['rotD50_T5_000'][i], dbacc['rotD50_T6_000'][i],
                           dbacc['rotD50_T7_000'][i], dbacc['rotD50_T8_000'][i],
                           dbacc['rotD50_T9_000'][i],
                           dbacc['rotD50_T10_000'][i]])
        sa_geo = None
        if source[i] == 'ESM':
            sa_geo = rotd50 / 981  # in g
        if source[i] == 'NGA-West2':
            sa_geo = rotd50  # already in g
        sa_list.append(sa_geo)
        if all(v > 0 for v in sa_geo):
            # print('Need to test if the screening of database is ok')
            if source[i] in allowed_database:
                if (source[i] == 'ESM' and is_free_field_esm[i] == 0) or \
                        (source[i] == 'NGA-West2' and
                         is_free_field_nga[i] == "I"):
                    if ((source[i] == 'ESM' and
                         allowed_recs_mag[0] <= event_mw[i] <=
                         allowed_recs_mag[1]) or
                            (source[i] == 'NGA-West2' and
                             allowed_recs_mag[0] <= event_mag[i] <=
                             allowed_recs_mag[1])):
                        if (allowed_depth[0] <= event_depth[i] <=
                                allowed_depth[1]):
                            if (allowed_recs_d[0] <= acc_distance[i] <=
                                    allowed_recs_d[1]):
                                if np.isnan(station_vs30[i]):
                                    if not pd.isnull(station_ec8[i]):
                                        if (station_ec8[i][0] in
                                                allowed_ec8_code or
                                                allowed_ec8_code == 'All'):
                                            if source[i] == 'ESM':
                                                allowed_index.append(i)
                                            if source[i] == 'NGA-West2':
                                                if 'ESM' in allowed_database:
                                                    if (epi_lon[i] < -31 or
                                                            epi_lon[i] > 70):
                                                        allowed_index.append(i)
                                                else:
                                                    allowed_index.append(i)
                                else:
                                    if (allowed_recs_vs30[0] <= station_vs30[i]
                                            < allowed_recs_vs30[1]):
                                        if source[i] == 'ESM':
                                            allowed_index.append(i)
                                        if source[i] == 'NGA-West2':
                                            if 'ESM' in allowed_database:
                                                if (epi_lon[i] < -31 or
                                                        epi_lon[i] > 70):
                                                    allowed_index.append(i)
                                            else:
                                                allowed_index.append(i)

    sa = np.vstack(sa_list)
    sa_known = sa[allowed_index]

    # count number of allowed spectra
    n_big = len(allowed_index)
    print(['Number of allowed ground motions = ', n_big])
    assert (n_big >= n_gm), \
        'Warning: there are not enough allowable ground motions'

    return [sa_known, ind_per, rec_per, n_big, allowed_index, event_id,
            station_code, source, record_sequence_number_nga, event_mw,
            event_mag, acc_distance, station_vs30, station_ec8]
