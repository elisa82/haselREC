# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# OpenSel is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenSel is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenSel. If not, see <http://www.gnu.org/licenses/>.

def create_esm_acc(folder):
    """
    """
    # Import libraries
    import glob
    from obspy.core import Stats
    import numpy as np

    filename_in = ''
    for l in range(1, 3):
        if folder.find('ESM/GR') > -1:
            file_ew = folder + '/*2.D.*'
            file_ns = folder + '/*3.D.*'
        else:
            file_ew = folder + '/*E.D.*'
            file_ns = folder + '/*N.D.*'
        if l == 1:
            filename_in = glob.glob(file_ew)[0]
        if l == 2:
            filename_in = glob.glob(file_ns)[0]

        headers = {}

        # read file
        fh = open(filename_in, 'rt')
        for i in range(64):
            key, value = fh.readline().strip().split(':', 1)
            headers[key.strip()] = value.strip()

        header = Stats()

        header['dyna'] = {}

        header['network'] = headers['NETWORK']
        header['station'] = headers['STATION_CODE']
        header['location'] = headers['LOCATION']
        header['channel'] = headers['STREAM']
        try:
            # use toUTCDateTime to convert from DYNA format
            header['starttime'] \
                = to_utc_date_time(headers
                                   ['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS'])
        except:
            header['starttime'] = to_utc_date_time('19700101_000000')
        header['sampling_rate'] = 1 / float(headers['SAMPLING_INTERVAL_S'])
        header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
        header['npts'] = int(headers['NDATA'])
        header['calib'] = 1  # not in file header

        # DYNA dict float data
        header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(
            headers['EVENT_LATITUDE_DEGREE'])
        header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(
            headers['EVENT_LONGITUDE_DEGREE'])
        header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
        header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
        header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
        header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
        header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(
            headers['STATION_LATITUDE_DEGREE'])
        header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(
            headers['STATION_LONGITUDE_DEGREE'])
        header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
        header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(
            headers['EPICENTRAL_DISTANCE_KM'])
        header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(
            headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
        header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
        header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(
            headers['INSTRUMENTAL_FREQUENCY_HZ'])
        header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(
            headers['INSTRUMENTAL_DAMPING'])
        header['dyna']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

        # data type is acceleration
        if headers['DATA_TYPE'] == "ACCELERATION" \
                or headers['DATA_TYPE'] == "ACCELERATION RESPONSE SPECTRUM":
            header['dyna']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
            header['dyna']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
        # data type is velocity
        if headers['DATA_TYPE'] == "VELOCITY" \
                or headers['DATA_TYPE'] == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
            header['dyna']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
            header['dyna']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
        # data type is displacement
        if headers['DATA_TYPE'] == "DISPLACEMENT" \
                or headers['DATA_TYPE'] == "DISPLACEMENT RESPONSE SPECTRUM":
            header['dyna']['PGD_CM'] = strtofloat(headers['PGD_CM'])
            header['dyna']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])

        header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(
            headers['LOW_CUT_FREQUENCY_HZ'])
        header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(
            headers['HIGH_CUT_FREQUENCY_HZ'])

        # DYNA dict int data
        header['dyna']['STATION_ELEVATION_M'] = strtoint(
            headers['STATION_ELEVATION_M'])
        header['dyna']['SENSOR_DEPTH_M'] = strtoint(headers['SENSOR_DEPTH_M'])
        header['dyna']['N_BIT_DIGITAL_CONVERTER'] = strtoint(
            headers['N_BIT_DIGITAL_CONVERTER'])
        header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

        # DYNA dict string data
        header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
        header['dyna']['EVENT_ID'] = headers['EVENT_ID']
        header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
        header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
        header['dyna']['MAGNITUDE_W_REFERENCE'] = headers[
            'MAGNITUDE_W_REFERENCE']
        header['dyna']['MAGNITUDE_L_REFERENCE'] = headers[
            'MAGNITUDE_L_REFERENCE']
        header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
        header['dyna']['STATION_NAME'] = headers['STATION_NAME']
        header['dyna']['SITE_CLASSIFICATION_EC8'] = headers[
            'SITE_CLASSIFICATION_EC8']
        header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers[
            'MORPHOLOGIC_CLASSIFICATION']
        header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers[
            'DATE_TIME_FIRST_SAMPLE_PRECISION']
        header['dyna']['UNITS'] = headers['UNITS']
        header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
        header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers[
            'INSTRUMENT_ANALOG/DIGITAL']
        header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
        header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
        header['dyna']['LATE_NORMAL_TRIGGERED'] = headers[
            'LATE/NORMAL_TRIGGERED']
        header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
        header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
        header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
        header['dyna']['PROCESSING'] = headers['PROCESSING']
        header['dyna']['DATA_LICENSE'] = headers['DATA_LICENSE']
        header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers[
            'DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
        header['dyna']['DATA_CITATION'] = headers['DATA_CITATION']
        header['dyna']['DATA_CREATOR'] = headers['DATA_CREATOR']
        header['dyna']['ORIGINAL_DATA_MEDIATOR_CITATION'] = headers[
            'ORIGINAL_DATA_MEDIATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_MEDIATOR'] = headers[
            'ORIGINAL_DATA_MEDIATOR']
        header['dyna']['ORIGINAL_DATA_CREATOR_CITATION'] = headers[
            'ORIGINAL_DATA_CREATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_CREATOR'] = headers[
            'ORIGINAL_DATA_CREATOR']
        header['dyna']['USER1'] = headers['USER1']
        header['dyna']['USER2'] = headers['USER2']
        header['dyna']['USER3'] = headers['USER3']
        header['dyna']['USER4'] = headers['USER4']
        header['dyna']['USER5'] = headers['USER5']

        # read data
        acc_data = np.loadtxt(fh, dtype='float32')
        fh.close()

        time1 = []
        time2 = []
        inp_acc1 = []
        inp_acc2 = []
        npts1 = []
        npts2 = []
        time = []
        for j in range(0, header['npts']):
            t = j * header['delta']
            time.append(t)

        if l == 1:
            inp_acc1 = np.asarray(acc_data) / 981  # in g
            # comp1=header['channel']
            npts1 = header['npts']
            time1 = time
        if l == 2:
            inp_acc2 = np.asarray(acc_data) / 981  # in g
            # comp2=header['channel']
            npts2 = header['npts']
            time2 = time

    return time1, time2, inp_acc1, inp_acc2, npts1, npts2


def to_utc_date_time(value):
    """
    """
    # Import libraries
    from obspy.core import UTCDateTime

    try:
        date, time = value.split('_')
    except ValueError:
        date = value

    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])

    hour = int(time[0:2])
    mins = int(time[2:4])
    secs = float(time[4:])

    return UTCDateTime(year, month, day, hour, mins) + secs


def strtofloat(sf):
    """
    """
    try:
        x = float(sf)
    except:
        return None
    return x


def strtoint(sf):
    """
    """
    try:
        x = int(sf)
    except:
        return None
    return x


def create_nga_acc(num_rec, path_nga_folder):
    """
    """
    # Import libraries
    import numpy as np

    # desc1 = ""
    # desc2 = ""
    for i in range(1, 3):
        file_acc = path_nga_folder + '/RSN' + str(num_rec) + '_' + str(
            i) + '.AT2'
        acc_data = []
        time = []
        with open(file_acc, 'r') as f:
            content = f.readlines()
        counter = 0
        for x in content:
            # if counter == 1:
            # if i == 1:
            # desc1 = x
            # desc1 = desc1[0:(len(x) - 1)]
            # if i == 2:
            # desc2 = x
            # desc2 = desc2[0:(len(x) - 1)]
            if counter == 3:
                row4val = x
                val = row4val.split()
                npts_comma = val[1]
                npts = int(npts_comma[0:len(npts_comma) - 1])
                if i == 1:
                    npts1 = npts
                if i == 2:
                    npts2 = npts
                dt = float(val[3])
                for j in range(0, npts):
                    t = j * dt
                    time.append(t)
                if i == 1:
                    time1 = time
                if i == 2:
                    time2 = time
            elif counter > 3:
                data = str(x).split()
                for value in data:
                    a = float(value)
                    acc_data.append(a)
                    if i == 1:
                        inp_acc1 = np.asarray(acc_data)
                    if i == 2:
                        inp_acc2 = np.asarray(acc_data)
            counter = counter + 1
    return time1, time2, inp_acc1, inp_acc2, npts1, npts2
