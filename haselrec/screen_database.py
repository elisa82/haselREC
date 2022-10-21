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
                    vs30, comp, radius_dist_type, radius_mag_type, selection_type):
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
    import h5py
    import sys
    from obspy.core import UTCDateTime

    known_per = np.array(
        [0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25,
         0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9,
         1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0,
         9.0, 10])
    sa = []
    if comp=='single-component' or selection_type=='code-spectrum':
        if 'NGA-West2' in allowed_database or 'KiK-net' in allowed_database:
            known_per=known_per[1:]
    allowed_index = []
    comp_allowed = []

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

    if radius_dist_type == 'both':
        allowed_recs_d = [dist_min_km, dist_max_km]
    elif radius_dist_type == 'right':
        allowed_recs_d = [mean_dist, dist_max_km]
    elif radius_dist_type == 'left':
        allowed_recs_d = [dist_min_km, mean_dist]
    else:
        sys.exit('The radius_dist_type is not allowed. It can be "both","right" or "left"')
    if radius_mag_type == 'both':
        allowed_recs_mag = [mean_mag - radius_mag, mean_mag + radius_mag]
    elif radius_mag_type == 'right':
        allowed_recs_mag = [mean_mag, mean_mag + radius_mag]
    elif radius_mag_type == 'left':
        allowed_recs_mag = [mean_mag - radius_mag, mean_mag]
    else:
        sys.exit('The radius_mag_type is not allowed. It can be "both","right" or "left"')    

    if allowed_recs_vs30 is None:
        if len(vs30) > 0:
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
        if len(vs30) > 0:
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

    event_id_tot = []
    station_code_tot = []
    event_mw_tot = []
    event_mag_tot = []
    acc_distance_tot = [] 
    station_vs30_tot = []
    station_ec8_tot = []
    source_tot = []
    record_sequence_number_tot = []
    fminNS2_tot = []
    fmaxNS2_tot = []
    fminEW2_tot = []
    fmaxEW2_tot = []
    epi_lon_tot = []
    epi_lat_tot = []
    origin_time_tot = []

    num_events_esm=0
    num_events_nga=0
    num_events_kiknet=0

    if 'ESM' in allowed_database:
        ESM_path=database_path+'/ESM_flatfile_SA.csv'
        dbacc = pd.read_csv(ESM_path, sep=';', engine='python')
        event_id_esm = dbacc['event_id']
        num_events_esm = len(event_id_esm)
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

        for i in np.arange(len(event_id_esm)):
            event_id_tot.append(event_id_esm[i])
            station_code_tot.append(station_code[i])
            event_mw_tot.append(event_mw[i])
            event_mag_tot.append(np.nan)
            acc_distance_tot.append(acc_distance[i])
            station_vs30_tot.append(station_vs30[i])
            station_ec8_tot.append(station_ec8[i])
            source_tot.append('ESM')
            record_sequence_number_tot.append(np.nan)
            fminNS2_tot.append(np.nan)
            fmaxNS2_tot.append(np.nan)
            fminEW2_tot.append(np.nan)
            fmaxEW2_tot.append(np.nan)
            epi_lon_tot.append(epi_lon[i])
            epi_lat_tot.append(epi_lat[i])
            origin_time_tot.append(UTCDateTime(int(dbacc['event_time'][i][6:10]), int(dbacc['event_time'][i][3:5]), int(dbacc['event_time'][i][0:2]), int(dbacc['event_time'][i][11:13]), int(dbacc['event_time'][i][14:16])))

            rec_ok=0
            #bad quality recording
            if event_id_esm[i] == 'IT-2013-0003' and station_code[i] == 'PSC': 
                rec_ok=0
            else:

                if (is_free_field[i] == 0):
                    if (allowed_recs_mag[0] <= event_mw[i] <= allowed_recs_mag[1]):
                        if (allowed_depth[0] <= event_depth[i] <= allowed_depth[1]):
                            if (allowed_recs_d[0] <= acc_distance[i] <= allowed_recs_d[1]):
                                if not np.isnan(station_vs30[i]):
                                    if (allowed_recs_vs30[0] <= station_vs30[i]
                                            < allowed_recs_vs30[1]):
                                        rec_ok=1
                                else:
                                    if allowed_ec8_code and not pd.isnull(station_ec8[i]):
                                        if (station_ec8[i][0] in
                                            allowed_ec8_code or
                                            allowed_ec8_code == 'All'):
                                                rec_ok=1
                                            
            if rec_ok==1:
                if(comp=='two-component' and selection_type=='conditional-spectrum'):
                    sa_entry = np.array([dbacc['rotD50_pga'][i], dbacc['rotD50_T0_010'][i],
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
                                       dbacc['rotD50_T9_000'][i], dbacc['rotD50_T10_000'][i]])

                    if all(v>0 for v in sa_entry):
                        sa_geo = sa_entry / 981  # in g
                        comp_allowed.append(0)
                        allowed_index.append(i)
                        sa.append(sa_geo)
                elif(comp=='single-component' or selection_type=='code-spectrum'):
                    sa_entry1 = np.array([np.abs(dbacc['U_pga'][i]), dbacc['U_T0_010'][i],
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
                                       dbacc['U_T9_000'][i], dbacc['U_T10_000'][i]])

                    sa_entry2 = np.array([np.abs(dbacc['V_pga'][i]), dbacc['V_T0_010'][i],
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
                                       dbacc['V_T9_000'][i], dbacc['V_T10_000'][i]])
                    if comp == 'single-component':
                        if all(v>0 for v in sa_entry1):
                            sa_entry1 = sa_entry1 / 981  # in g
                            allowed_index.append(i)
                            comp_allowed.append(1)
                            if 'NGA-West2' in allowed_database or 'KiK-net' in allowed_database:
                                sa_entry1=sa_entry1[1:]
                            sa.append(sa_entry1)
                        if all(v>0 for v in sa_entry2):
                            sa_entry2 = sa_entry2 / 981  # in g
                            allowed_index.append(i)
                            comp_allowed.append(2)
                            if 'NGA-West2' in allowed_database or 'KiK-net' in allowed_database:
                                sa_entry2=sa_entry2[1:]
                            sa.append(sa_entry2)
                    elif comp == 'two-component':
                        if all(v>0 for v in sa_entry1):
                            sa_entry1 = sa_entry1 / 981  # in g
                            sa_entry2 = sa_entry2 / 981  # in g
                            allowed_index.append(i)
                            comp_allowed.append(0)
                            if 'NGA-West2' in allowed_database or 'KiK-net' in allowed_database:
                                sa_entry1=sa_entry1[1:]
                                sa_entry2=sa_entry2[1:]
                            sa.append(np.sqrt(sa_entry1**2+sa_entry2**2))

    if 'NGA-West2' in allowed_database:
        index0_database=num_events_esm

        NGAWest2_path=database_path+'/Updated_NGA_West2_Flatfile_RotD50_d050_public_version.csv'
        dbacc = pd.read_csv(NGAWest2_path, sep=';', engine='python')
        event_id_nga = dbacc['EQID']
        num_events_nga = len(event_id_nga)
        event_mag = dbacc['Earthquake Magnitude']
        record_sequence_number_nga = dbacc['Record Sequence Number']
        station_vs30 = dbacc['Vs30 (m/s) selected for analysis']
        acc_distance = dbacc['EpiD (km)']
        station_code = dbacc['Station ID  No.']
        event_depth = dbacc['Hypocenter Depth (km)']
        is_free_field = dbacc["GMX's C1"]
        epi_lon = dbacc['Hypocenter Longitude (deg)']
        epi_lat = dbacc['Hypocenter Latitude (deg)']

        for i in np.arange(len(event_id_nga)):
            event_id_tot.append(event_id_nga[i])
            station_code_tot.append(station_code[i])
            event_mw_tot.append(np.nan)
            event_mag_tot.append(event_mag[i])
            acc_distance_tot.append(acc_distance[i])
            station_vs30_tot.append(station_vs30[i])
            station_ec8_tot.append(np.nan)
            source_tot.append('NGA-West2')
            record_sequence_number_tot.append(int(record_sequence_number_nga[i]))
            fminNS2_tot.append(np.nan)
            fmaxNS2_tot.append(np.nan)
            fminEW2_tot.append(np.nan)
            fmaxEW2_tot.append(np.nan)
            epi_lon_tot.append(epi_lon[i])
            epi_lat_tot.append(epi_lat[i])
            year=str(dbacc['YEAR'][i])
            mody=str(dbacc['MODY'][i])
            hrmn=str(dbacc['HRMN'][i])
            if(year=='-999'):
                year='1000'
            if(mody=='-999'):
                mody='0615' #arbitrary date in the middle of the year
            if(hrmn=='-999'):
                hrmn='0000'
            if(len(mody)==3):
                mody='0'+str(dbacc['MODY'][i])
            if(len(hrmn)==3):
                hrmn='0'+str(dbacc['HRMN'][i])
            elif(len(hrmn)==2):
                hrmn='00'+str(dbacc['HRMN'][i])
            else:
                hrmn='000'+hrmn
            origin_time_tot.append(UTCDateTime(int(year), int(mody[0:2]), int(mody[2:4]), int(hrmn[0:2]), int(hrmn[2:4])))
            rec_ok=0
            if(is_free_field[i] == "I"):
                if(allowed_recs_mag[0] <= event_mag[i] <= allowed_recs_mag[1]):
                    if (allowed_depth[0] <= event_depth[i] <= allowed_depth[1]):
                        if (allowed_recs_d[0] <= acc_distance[i] <= allowed_recs_d[1]):
                            if (allowed_recs_vs30[0] <= station_vs30[i] < allowed_recs_vs30[1]):
                                if (record_sequence_number_nga[i] != 789):
                                    inside=0
                                    if 'ESM' in allowed_database:
                                        if (epi_lon[i] > -9 and epi_lon[i] < 64) and (epi_lat[i] > 28 and epi_lat[i] < 52):
                                            inside=1
                                    if 'KiK-net' in allowed_database: 
                                        if (epi_lon[i] > 120 and epi_lon[i] < 160) and (epi_lat[i] > 18 and epi_lat[i] < 55):
                                            inside=1
                                    if inside == 0:
                                        rec_ok=1

            if comp=='two-component' and selection_type=='conditional-spectrum':
                if rec_ok == 1:
                    sa_entry = np.array([dbacc['PGA (g)'][i], dbacc['T0.010S'][i],
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
                                       dbacc['T9.000S'][i], dbacc['T10.000S'][i]])
                    if all(v>0 for v in sa_entry):
                        sa_geo = sa_entry # already in g
                        comp_allowed.append(0)
                        allowed_index.append(i+index0_database)
                        sa.append(sa_geo)
            elif comp=='single-component' or selection_type=='code-spectrum':
                if rec_ok == 1:
                    sa_entry1 = np.array([dbacc['Sa1_T0.010S'][i],
                                       dbacc['Sa1_T0.025S'][i], dbacc['Sa1_T0.040S'][i],
                                       dbacc['Sa1_T0.050S'][i], dbacc['Sa1_T0.070S'][i],
                                       dbacc['Sa1_T0.100S'][i], dbacc['Sa1_T0.150S'][i],
                                       dbacc['Sa1_T0.200S'][i], dbacc['Sa1_T0.250S'][i],
                                       dbacc['Sa1_T0.300S'][i], dbacc['Sa1_T0.350S'][i],
                                       dbacc['Sa1_T0.400S'][i], dbacc['Sa1_T0.450S'][i],
                                       dbacc['Sa1_T0.500S'][i], dbacc['Sa1_T0.600S'][i],
                                       dbacc['Sa1_T0.700S'][i], dbacc['Sa1_T0.750S'][i],
                                       dbacc['Sa1_T0.800S'][i], dbacc['Sa1_T0.900S'][i],
                                       dbacc['Sa1_T1.000S'][i], dbacc['Sa1_T1.200S'][i],
                                       dbacc['Sa1_T1.400S'][i], dbacc['Sa1_T1.600S'][i],
                                       dbacc['Sa1_T1.800S'][i], dbacc['Sa1_T2.000S'][i],
                                       dbacc['Sa1_T2.500S'][i], dbacc['Sa1_T3.000S'][i],
                                       dbacc['Sa1_T3.500S'][i], dbacc['Sa1_T4.000S'][i],
                                       dbacc['Sa1_T5.000S'][i], dbacc['Sa1_T6.000S'][i],
                                       dbacc['Sa1_T7.000S'][i], dbacc['Sa1_T8.000S'][i],
                                       dbacc['Sa1_T9.000S'][i], dbacc['Sa1_T10.000S'][i]])
                    sa_entry2 = np.array([dbacc['Sa2_T0.010S'][i],
                                       dbacc['Sa2_T0.025S'][i], dbacc['Sa2_T0.040S'][i],
                                       dbacc['Sa2_T0.050S'][i], dbacc['Sa2_T0.070S'][i],
                                       dbacc['Sa2_T0.100S'][i], dbacc['Sa2_T0.150S'][i],
                                       dbacc['Sa2_T0.200S'][i], dbacc['Sa2_T0.250S'][i],
                                       dbacc['Sa2_T0.300S'][i], dbacc['Sa2_T0.350S'][i],
                                       dbacc['Sa2_T0.400S'][i], dbacc['Sa2_T0.450S'][i],
                                       dbacc['Sa2_T0.500S'][i], dbacc['Sa2_T0.600S'][i],
                                       dbacc['Sa2_T0.700S'][i], dbacc['Sa2_T0.750S'][i],
                                       dbacc['Sa2_T0.800S'][i], dbacc['Sa2_T0.900S'][i],
                                       dbacc['Sa2_T1.000S'][i], dbacc['Sa2_T1.200S'][i],
                                       dbacc['Sa2_T1.400S'][i], dbacc['Sa2_T1.600S'][i],
                                       dbacc['Sa2_T1.800S'][i], dbacc['Sa2_T2.000S'][i],
                                       dbacc['Sa2_T2.500S'][i], dbacc['Sa2_T3.000S'][i],
                                       dbacc['Sa2_T3.500S'][i], dbacc['Sa2_T4.000S'][i],
                                       dbacc['Sa2_T5.000S'][i], dbacc['Sa2_T6.000S'][i],
                                       dbacc['Sa2_T7.000S'][i], dbacc['Sa2_T8.000S'][i],
                                       dbacc['Sa2_T9.000S'][i], dbacc['Sa2_T10.000S'][i]])

                    if comp == 'single-component':
                        if all(v>0 for v in sa_entry1):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(1)
                            sa.append(sa_entry1)
                        if all(v>0 for v in sa_entry2):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(2)
                            sa.append(sa_entry2)
                    elif comp == 'two-component':
                        if all(v>0 for v in sa_entry1):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(0)
                            sa.append(np.sqrt(sa_entry1**2+sa_entry2**2))

    if 'KiK-net' in allowed_database:
        index0_database=num_events_esm+num_events_nga

        kiknet_path=database_path+'/attributes.csv'
        dbacc = pd.read_csv(kiknet_path, sep=',', engine='python')
        address = dbacc['Adress']
        num_events_kiknet = len(address)
        event_mag = dbacc['MT_Magnitude_']
        record_number = dbacc['Code']
        station_vs30 = dbacc['Vs30']
        acc_distance = dbacc['repi_0']
        station_code = dbacc['station']
        event_depth = dbacc['depth']
        tectonic_environment = dbacc['Tectonic_Garcia_']
        fminNS2 = dbacc['fLow_NS2']
        fmaxNS2 = dbacc['fHigh_NS2']
        fminEW2 = dbacc['fLow_EW2']
        fmaxEW2 = dbacc['fHigh_EW2']
        epi_lon = dbacc['lon']
        epi_lat = dbacc['lat']

        fileName = database_path+'/Database_small.hdf5'
        fn = h5py.File(fileName,'r')
        a = ['Address']
        a.extend(list(fn['T'][:]))

        for i in np.arange(len(address)):
            event_id_tot.append(address[i][3:17])
            station_code_tot.append(station_code[i])
            event_mw_tot.append(np.nan)
            event_mag_tot.append(event_mag[i])
            acc_distance_tot.append(acc_distance[i])
            station_vs30_tot.append(station_vs30[i])
            station_ec8_tot.append(np.nan)
            source_tot.append('KiK-net')
            record_sequence_number_tot.append(record_number[i])
            fminNS2_tot.append(fminNS2[i])
            fmaxNS2_tot.append(fmaxNS2[i])
            fminEW2_tot.append(fminEW2[i])
            fmaxEW2_tot.append(fmaxEW2[i])
            epi_lon_tot.append(epi_lon[i])
            epi_lat_tot.append(epi_lat[i])
            origin_time_tot.append(UTCDateTime(int(dbacc['originDate'][i].split('/')[2]),int(dbacc['originDate'][i].split('/')[0]),int(dbacc['originDate'][i].split('/')[1]),int(dbacc['originTime'][i].split(':')[0]),int(dbacc['originTime'][i].split(':')[1]),int(str(dbacc['originTime'][i]).split(':')[2][0:2])))

            index_kiknet=[1,4,10,16,21,27,33,38,41,45,48,51,54,57,59,62,63,64,66,68,70,72,74,76,78,81,84,87,90,95,97,99,101,103,105]
            rec_ok=0
            if(tectonic_environment[i] == 2):
                if(allowed_recs_mag[0] <= event_mag[i] <= allowed_recs_mag[1]):
                    if (allowed_depth[0] <= event_depth[i] <= allowed_depth[1]):
                        if (allowed_recs_d[0] <= acc_distance[i] <= allowed_recs_d[1]):
                            if (allowed_recs_vs30[0] <= station_vs30[i] < allowed_recs_vs30[1]):
                                if (record_number[i] != 'TCGH141103230712'):
                                    rec_ok = 1

            if comp=='two-component' and selection_type=='conditional-spectrum':
                sys.exit('Error: KiK-net database does not have RotD50 component. The geometric mean will be added in the future')
            elif comp=='single-component' or selection_type=='code-spectrum':
                if rec_ok == 1:
                    sa_entry1=[]
                    sa_entry2=[]
                    direction='EW2'
                    b = [address.iloc[i][2:-2]+ '/' + direction]
                    b.extend(list(fn[b[0]]['PSa'][:]))
                    b=np.asarray(b)
                    for kk in range(len(index_kiknet)):
                        sa_entry1.append(b[index_kiknet[kk]].astype(np.float64))
                    sa_entry1  = np.asarray(sa_entry1)
                    direction='NS2'
                    b = [address.iloc[i][2:-2]+ '/' + direction]
                    b.extend(list(fn[b[0]]['PSa'][:]))
                    b=np.asarray(b)
                    for kk in range(len(index_kiknet)):
                        sa_entry2.append(b[index_kiknet[kk]].astype(np.float64))
                    sa_entry2  = np.asarray(sa_entry2)
                    if comp == 'single-component':
                        if all(v>0 for v in sa_entry1):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(1)
                            sa.append(sa_entry1)
                        if all(v>0 for v in sa_entry2):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(2)
                            sa.append(sa_entry2)
                    elif comp == 'two-component':
                        if all(v>0 for v in sa_entry1):
                            allowed_index.append(i+index0_database)
                            comp_allowed.append(0)
                            sa.append(np.sqrt(sa_entry1**2+sa_entry2**2))


    sa_known=np.vstack(sa)

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
    record_sequence_number_tot=np.asarray(record_sequence_number_tot)
    fminNS2_tot=np.asarray(fminNS2_tot)
    fmaxNS2_tot=np.asarray(fmaxNS2_tot)
    fminEW2_tot=np.asarray(fminEW2_tot)
    fmaxEW2_tot=np.asarray(fmaxEW2_tot)
    epi_lon_tot=np.asarray(epi_lon_tot)
    epi_lat_tot=np.asarray(epi_lat_tot)
    origin_time_tot=np.asarray(origin_time_tot)

    cluster,flagvector = declustering(event_id_tot,event_mw_tot,
                                event_mag_tot,origin_time_tot,
                                epi_lon_tot,epi_lat_tot)
    return (sa_known, ind_per, rec_per, n_big, allowed_index, event_id_tot,
           station_code_tot, source_tot, record_sequence_number_tot, event_mw_tot,
           event_mag_tot, acc_distance_tot, station_vs30_tot, station_ec8_tot,
           comp_allowed, fminNS2_tot, fmaxNS2_tot, fminEW2_tot, fmaxEW2_tot, cluster)

def declustering(event_id_tot,event_mw_tot,event_mag_tot,origin_time_tot,epi_lon_tot,epi_lat_tot):

    import numpy as np

    neq = len(event_id_tot)
    mag=np.zeros((neq))
    sw_space=np.zeros((neq))
    sw_time=np.zeros((neq))
    julian_date=np.zeros((neq))
    for i in np.arange(neq):
        mag[i] = define_magnitude(i,event_mw_tot,event_mag_tot)
        # Get space and time windows corresponding to each event
        sw_space[i], sw_time[i] = GK(mag[i])
        julian_date[i] = to_julian_date(origin_time_tot[i].year, origin_time_tot[i].month, origin_time_tot[i].day)
    # Initial Position Identifier
    eqid = np.arange(0, neq, 1)
    # Pre-allocate cluster index vectors
    vcl = np.zeros(neq, dtype=int)
    # Sort magnitudes into descending order
    id0 = np.flipud(np.argsort(mag, kind='heapsort'))
    longitude = epi_lon_tot[id0]
    latitude = epi_lat_tot[id0]
    sw_space = sw_space[id0]
    sw_time = sw_time[id0]
    julian_date = julian_date[id0]
    eqid = eqid[id0]
    flagvector = np.zeros(neq, dtype=int)

    # Begin cluster identification
    clust_index = 0
    for i in range(0, neq - 1):
        if vcl[i] == 0:
            # Find Events inside both fore- and aftershock time windows
            dt = julian_date - julian_date[i]
            vsel = np.logical_and(
                vcl == 0,
                np.logical_and(
                    dt >= (-sw_time[i]),
                    dt <= sw_time[i]))
            # Of those events inside time window,
            # find those inside distance window
            vsel1 = haversine(longitude[vsel],
                              latitude[vsel],
                              longitude[i],
                              latitude[i]) <= sw_space[i]
            vsel[vsel] = vsel1[:, 0]
            temp_vsel = np.copy(vsel)
            temp_vsel[i] = False
            if any(temp_vsel):
                # Allocate a cluster number
                vcl[vsel] = clust_index + 1
                flagvector[vsel] = 1
                # For those events in the cluster before the main event,
                # flagvector is equal to -1
                temp_vsel[dt >= 0.0] = False
                flagvector[temp_vsel] = -1
                flagvector[i] = 0
                clust_index += 1

    # Re-sort the catalog_matrix into original order
    id1 = np.argsort(eqid, kind='heapsort')
    eqid = eqid[id1]
    vcl = vcl[id1]
    flagvector = flagvector[id1]

    return vcl, flagvector


def GK(mag):
    #Gardner-Knopoff (1974)
    dR = 10**(0.1238*mag + 0.983) #in km
    if mag >=6.5: #in days
        dT = 10**(0.032*mag+2.7389)
    else:
        dT = 10**(0.5409*mag-0.547)
    return dR,dT

def define_magnitude(index,event_mw,event_mag):
    
    import numpy as np

    if not np.isnan(event_mw[index]):
        mag=event_mw[index]
    elif not np.isnan(event_mag[index]):
        mag=event_mag[index]
    else:
        mag = 0
    return mag

def to_julian_date(I, J, K):

    """
    """
    JD=K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4

    return JD


def haversine(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
    """
    Allows to calculate geographical distance
    using the haversine formula.

    :param lon1: longitude of the first set of locations
    :type lon1: numpy.ndarray
    :param lat1: latitude of the frist set of locations
    :type lat1: numpy.ndarray
    :param lon2: longitude of the second set of locations
    :type lon2: numpy.float64
    :param lat2: latitude of the second set of locations
    :type lat2: numpy.float64
    :keyword radians: states if locations are given in terms of radians
    :type radians: bool
    :keyword earth_rad: radius of the earth in km
    :type earth_rad: float
    :returns: geographical distance in km
    :rtype: numpy.ndarray
    """

    import numpy as np

    if not radians:
        cfact = np.pi / 180.
        lon1 = cfact * lon1
        lat1 = cfact * lat1
        lon2 = cfact * lon2
        lat2 = cfact * lat2

    # Number of locations in each set of points
    if not np.shape(lon1):
        nlocs1 = 1
        lon1 = np.array([lon1])
        lat1 = np.array([lat1])
    else:
        nlocs1 = np.max(np.shape(lon1))
    if not np.shape(lon2):
        nlocs2 = 1
        lon2 = np.array([lon2])
        lat2 = np.array([lat2])
    else:
        nlocs2 = np.max(np.shape(lon2))
    # Pre-allocate array
    distance = np.zeros((nlocs1, nlocs2))
    i = 0
    while i < nlocs2:
        # Perform distance calculation
        dlat = lat1 - lat2[i]
        dlon = lon1 - lon2[i]
        aval = (np.sin(dlat / 2.) ** 2.) + (np.cos(lat1) * np.cos(lat2[i]) *
                                            (np.sin(dlon / 2.) ** 2.))
        distance[:, i] = (2. * earth_rad * np.arctan2(np.sqrt(aval),
                                                      np.sqrt(1 - aval))).T
        i += 1
    return distance
