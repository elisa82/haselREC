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

def screen_database(database_path, allowed_database, allowed_recs_vs30,
                    radius_dist, dist_range_input, radius_mag, mean_dist, mean_mag,
                    allowed_ec8_code, target_periods, n_gm, allowed_depth,
                    vs30, bgmpe):
    """
    Screen the database of candidate ground motion to select only appropriate
    ground motions. The screening criteria are:

        - database (NGAWest2, ESM or both); When both databases are considered,
          ground motions from the NGA-West2 database are retained only if
          recorded at stations located outside the geographical area covered by
          the ESM database;
        - magnitude range, defined as a symmetric interval around the mean
          magnitude from the disaggregation analysis;
        - distance range, defined as a symmetric interval around the mean
          distance from the disaggregation analysis;
        - range of allowed `vs30`. If not defined, it is set by the code
          according to the `vs30` of the site, following the `vs30` limit values
          associated to EC8 soil categories;
        - range of allowed EC8 codes. If not defined, they are set according to
          the vs30 of the site. If both `vs30` and EC8 soil classes criteria are
          specified, preference is given to the `vs30`; therefore, EC8 soil
          category criterium is considered only if the `vs30` of the station is
          not specified;
        - range of allowed focal depths;
        - only free-field ground motions are retained.
    """
    # Import libraries
    import numpy as np
    import pandas as pd

    known_per = np.array(
        [0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25,
         0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9,
         1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0,
         9.0, 10])

    # Match periods (known periods and target periods for error computations) 
    # save the indices of the matched periods in known_per
    ind_per = np.zeros((len(target_periods), 1), dtype=int)
    for i in np.arange(len(target_periods)):
        ind_per[i] = np.argmin(np.abs(known_per - target_periods[i]))

    # Remove any repeated values from TgtPer and redefine TgtPer as periods 
    # provided in databases
    ind_per = np.unique(ind_per)
    rec_per = known_per[ind_per]

    if dist_range_input is None:
        dist_range=[-99999.,99999.]
    else:
        dist_range=dist_range_input
    if ( (mean_dist - radius_dist) < dist_range[0] ):
        dist_min_km = dist_range[0]
    else:
        dist_min_km = mean_dist - radius_dist
    if ( (mean_dist + radius_dist) > dist_range[1] ):
        dist_max_km = dist_range[1]
    else:
        dist_max_km = mean_dist + radius_dist

    allowed_recs_d = [dist_min_km, dist_max_km]
    if(bgmpe()=='[Ambraseys1996]'):
        mean_mw=np.exp(1.421+0.108*mean_mag)-1.863
    else:
        mean_mw=mean_mag
    allowed_recs_mag = [mean_mw - radius_mag, mean_mw + radius_mag]

    if allowed_recs_vs30 is None:
        if vs30 >= 800.0:
            allowed_recs_vs30 = [800.0, 3000.0]
        elif 360. < vs30 < 800.:
            allowed_recs_vs30 = [360.0, 800.0]
        elif vs30 ==360.:
            allowed_recs_vs30 = [180.0, 800.0]
        elif 180. < vs30 < 360.:
            allowed_recs_vs30 = [180.0, 360.0]
        elif vs30 == 180.:
            allowed_recs_vs30 = [0.0, 360.0]
        else:
            allowed_recs_vs30 = [0.0, 180.0]

    if allowed_ec8_code is None:
        if vs30 >= 800.0:
            allowed_ec8_code = 'A'
        elif 360. < vs30 < 800.:
            allowed_ec8_code = 'B'
        elif vs30 ==360.:
            allowed_ec8_code = ['B', 'C']
        elif 180. < vs30 < 360.:
            allowed_ec8_code = 'C'
        elif vs30 == 180.:
            allowed_ec8_code = ['C', 'D']
        else:
            allowed_ec8_code = 'D'

    sa_geo_list = []
    sa_hor1_list = []
    sa_hor2_list = []
    allowed_index = []

    event_id_tot = []
    station_code_tot = []
    event_mw_tot = []
    event_mag_tot = []
    acc_distance_tot = [] 
    station_vs30_tot = []
    station_ec8_tot = []
    source_tot = []
    record_sequence_number_nga_tot = []

    if 'ESM' in allowed_database:

        ESM_path=database_path+'/ESM_flatfile_SA.csv'
        dbacc = pd.read_csv(ESM_path, sep=';', engine='python')
        event_id = dbacc['event_id']
        event_mw = dbacc['Mw']
        station_ec8 = dbacc['ec8_code']
        station_vs30 = dbacc['vs30_m_sec']
        acc_distance = dbacc['epi_dist']
        station_code = dbacc['station_code']
        event_depth = dbacc['ev_depth_km']
        # sensor_depth = dbacc['sensor_depth_m']
        is_free_field = dbacc['proximity_code']
        epi_lon = dbacc['ev_longitude']
        epi_lat = dbacc['ev_latitude']

        for i in np.arange(len(event_id)):
            event_id_tot.append(event_id[i])
            station_code_tot.append(station_code[i])
            event_mw_tot.append(event_mw[i])
            event_mag_tot.append(np.nan)
            acc_distance_tot.append(acc_distance[i])
            station_vs30_tot.append(station_vs30[i])
            station_ec8_tot.append(station_ec8[i])
            source_tot.append('ESM')
            record_sequence_number_nga_tot.append(np.nan)

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

            hor1 = np.array([dbacc['U_pga'][i], dbacc['U_T0_010'][i],
                               dbacc['U_T0_025'][i], dbacc['U_T0_040'][i],
                               dbacc['U_T0_050'][i], dbacc['U_T0_070'][i],
                               dbacc['U_T0_100'][i], dbacc['U_T0_150'][i],
                               dbacc['U_T0_200'][i], dbacc['U_T0_250'][i],
                               dbacc['U_T0_300'][i], dbacc['U_T0_350'][i],
                               dbacc['U_T0_400'][i], dbacc['U_T0_450'][i],
                               dbacc['U_T0_500'][i], dbacc['U_T0_600'][i],
                               dbacc['U_T0_700'][i], dbacc['U_T0_750'][i],
                               dbacc['U_T0_800'][i], dbacc['U_T0_900'][i],
                               dbacc['U_T1_000'][i], dbacc['U_T1_200'][i],
                               dbacc['U_T1_400'][i], dbacc['U_T1_600'][i],
                               dbacc['U_T1_800'][i], dbacc['U_T2_000'][i],
                               dbacc['U_T2_500'][i], dbacc['U_T3_000'][i],
                               dbacc['U_T3_500'][i], dbacc['U_T4_000'][i],
                               dbacc['U_T5_000'][i], dbacc['U_T6_000'][i],
                               dbacc['U_T7_000'][i], dbacc['U_T8_000'][i],
                               dbacc['U_T9_000'][i],
                               dbacc['U_T10_000'][i]])

            hor2 = np.array([dbacc['V_pga'][i], dbacc['V_T0_010'][i],
                               dbacc['V_T0_025'][i], dbacc['V_T0_040'][i],
                               dbacc['V_T0_050'][i], dbacc['V_T0_070'][i],
                               dbacc['V_T0_100'][i], dbacc['V_T0_150'][i],
                               dbacc['V_T0_200'][i], dbacc['V_T0_250'][i],
                               dbacc['V_T0_300'][i], dbacc['V_T0_350'][i],
                               dbacc['V_T0_400'][i], dbacc['V_T0_450'][i],
                               dbacc['V_T0_500'][i], dbacc['V_T0_600'][i],
                               dbacc['V_T0_700'][i], dbacc['V_T0_750'][i],
                               dbacc['V_T0_800'][i], dbacc['V_T0_900'][i],
                               dbacc['V_T1_000'][i], dbacc['V_T1_200'][i],
                               dbacc['V_T1_400'][i], dbacc['V_T1_600'][i],
                               dbacc['V_T1_800'][i], dbacc['V_T2_000'][i],
                               dbacc['V_T2_500'][i], dbacc['V_T3_000'][i],
                               dbacc['V_T3_500'][i], dbacc['V_T4_000'][i],
                               dbacc['V_T5_000'][i], dbacc['V_T6_000'][i],
                               dbacc['V_T7_000'][i], dbacc['V_T8_000'][i],
                               dbacc['V_T9_000'][i],
                               dbacc['V_T10_000'][i]])

            sa_geo_ESM = None
            sa_hor1_ESM = None
            sa_hor2_ESM = None
            sa_geo_ESM = rotd50 / 981  # in g
            sa_hor1_ESM = hor1 / 981  # in g
            sa_hor2_ESM = hor2 / 981  # in g
            sa_geo_list.append(sa_geo_ESM)
            sa_hor1_list.append(sa_hor1_ESM)
            sa_hor2_list.append(sa_hor2_ESM)
            if all(v > 0 for v in sa_geo_ESM):
                if (is_free_field[i] == 0):
                    if (allowed_recs_mag[0] <= event_mw[i] <= allowed_recs_mag[1]):
                        if (allowed_depth[0] <= event_depth[i] <= allowed_depth[1]):
                            if (allowed_recs_d[0] <= acc_distance[i] <= allowed_recs_d[1]):
                                if np.isnan(station_vs30[i]):
                                    if not pd.isnull(station_ec8[i]):
                                        if (station_ec8[i][0] in
                                                allowed_ec8_code or
                                                allowed_ec8_code == 'All'):
                                                allowed_index.append(i)
                                else:
                                    if (allowed_recs_vs30[0] <= station_vs30[i]
                                            < allowed_recs_vs30[1]):
                                            allowed_index.append(i)

    if 'NGA-West2' in allowed_database:
        if('ESM' in allowed_database):
            index0_database=len(event_id)
        else:
            index0_database=0

        NGAWest2_path=database_path+'/Updated_NGA_West2_Flatfile_RotD50_d050_public_version.csv'
        dbacc = pd.read_csv(NGAWest2_path, sep=';', engine='python')
        event_id = dbacc['EQID']
        event_mag = dbacc['Earthquake Magnitude']
        record_sequence_number_nga = dbacc['Record Sequence Number']
        station_vs30 = dbacc['Vs30 (m/s) selected for analysis']
        acc_distance = dbacc['EpiD (km)']
        station_code = dbacc['Station ID  No.']
        event_depth = dbacc['Hypocenter Depth (km)']
        is_free_field = dbacc["GMX's C1"]
        epi_lon = dbacc['Hypocenter Longitude (deg)']
        epi_lat = dbacc['Hypocenter Latitude (deg)']

        for i in np.arange(index0_database,index0_database+len(event_id)):
            event_id_tot.append(event_id[i])
            station_code_tot.append(station_code[i])
            event_mw_tot.append(np.nan)
            event_mag_tot.append(event_mag[i])
            acc_distance_tot.append(acc_distance[i])
            station_vs30_tot.append(station_vs30[i])
            station_ec8_tot.append(np.nan)
            source_tot.append('NGA-West2')
            record_sequence_number_nga_tot.append(record_sequence_number_nga[i])
            rotd50 = np.array([dbacc['PGA (g)'][i], dbacc['T0.010S'][i],
                               dbacc['T0.025S'][i], dbacc['T0.040S'][i],
                               dbacc['T0.050S'][i], dbacc['T0.070S'][i],
                               dbacc['T0.100S'][i], dbacc['T0.150S'][i],
                               dbacc['T0.200S'][i], dbacc['T0.250S'][i],
                               dbacc['T0.300S'][i], dbacc['T0.350S'][i],
                               dbacc['T0.400S'][i], dbacc['T0.450S'][i],
                               dbacc['T0.500S'][i], dbacc['T0.600S'][i],
                               dbacc['T0.700S'][i], dbacc['T0.750S'][i],
                               dbacc['T0.800S'][i], dbacc['T0.900S'][i],
                               dbacc['T1.000S'][i], dbacc['T1.200S'][i],
                               dbacc['T1.400S'][i], dbacc['T1.600S'][i],
                               dbacc['T1.800S'][i], dbacc['T2.000S'][i],
                               dbacc['T2.500S'][i], dbacc['T3.000S'][i],
                               dbacc['T3.500S'][i], dbacc['T4.000S'][i],
                               dbacc['T5.000S'][i], dbacc['T6.000S'][i],
                               dbacc['T7.000S'][i], dbacc['T8.000S'][i],
                               dbacc['T9.000S'][i],
                               dbacc['T10.000S'][i]])

            sa_geo_NGAW2 = None
            sa_geo_NGAW2 = rotd50 # already in g
            sa_geo_list.append(sa_geo_NGAW2)
            if all(v > 0 for v in sa_geo_NGAW2):
                if(is_free_field[i] == "I"):
                    if(allowed_recs_mag[0] <= event_mag[i] <= allowed_recs_mag[1]):
                        if (allowed_depth[0] <= event_depth[i] <= allowed_depth[1]):
                            if (allowed_recs_d[0] <= acc_distance[i] <= allowed_recs_d[1]):
                                if (allowed_recs_vs30[0] <= station_vs30[i] < allowed_recs_vs30[1]):
                                    if 'ESM' in allowed_database:
                                        if (epi_lon[i] < -31 or epi_lon[i] > 70):
                                            allowed_index.append(i)
                                    else:
                                        allowed_index.append(i)

    sa = np.vstack(sa_geo_list)
    sa_known = sa[allowed_index]

    # count number of allowed spectra
    n_big = len(allowed_index)
    print(['Number of allowed ground motions = ', n_big])
    assert (n_big >= n_gm), \
        'Warning: there are not enough allowable ground motions'

    event_id_tot=np.asarray(event_id_tot)
    station_code_tot=np.asarray(station_code_tot)
    event_mw_tot=np.asarray(event_mw_tot)
    event_mag_tot=np.asarray(event_mag_tot)
    acc_distance_tot=np.asarray(acc_distance_tot)
    station_vs30_tot=np.asarray(station_vs30_tot)
    station_ec8_tot=np.asarray(station_ec8_tot)
    source_tot=np.asarray(source_tot)
    record_sequence_number_nga_tot=np.asarray(record_sequence_number_nga_tot)

    return [sa_known, ind_per, rec_per, n_big, allowed_index, event_id_tot,
            station_code_tot, source_tot, record_sequence_number_nga_tot, event_mw_tot,
            event_mag_tot, acc_distance_tot, station_vs30_tot, station_ec8_tot]
