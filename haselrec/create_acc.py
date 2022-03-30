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

def create_esm_acc(folder, event, station, num):

    """

    ESM recordings can be stored in advance or automatically downloaded from
    internet using a token file (:code:`token.txt`). To obtain the token file you need
    at first to register at: `https://esm-db.eu/` and then you can run the command::

        curl -X POST -F 'message={"user_email": "email","user_password": "password"}
        ' "https://esm-db.eu/esmws/generate-signed-message/1/query" > token.txt

    """

    # Import libraries
    import glob
    from obspy.core import Stats
    import numpy as np
    import os
    from zipfile import ZipFile
    import requests
    import sys

    if not os.path.isdir(folder):
        zip_output = 'output_' + str(num) + '.zip'

        params = (
            ('eventid', event),
            ('data-type', 'ACC'),
            ('station', station),
            ('format', 'ascii'),
        )

        files = {
            'message': ('path/to/token.txt', open('token.txt', 'rb')),
        }

        headers = {'Authorization': 'token {}'.format('token.txt')}

        url = 'https://esm-db.eu/esmws/eventdata/1/query'

        req = requests.post(url=url, params=params, files=files)

        if req.status_code == 200:
            with open(zip_output, "wb") as zf:
                zf.write(req.content)
        else:
            if req.status_code == 403:
                sys.exit('Problem with ESM download. Maybe the token is no longer valid')
            else:
                sys.exit('Problem with ESM download. Status code: '+str(req.status_code))

        with ZipFile(zip_output, 'r') as zipObj:
            zipObj.extractall(folder)
        os.remove(zip_output)

    time1 = []
    time2 = []
    inp_acc1 = []
    inp_acc2 = []
    npts1 = []
    npts2 = []

    filename_in = ''
    for i in range(1, 3):
        if folder.find('ESM/GR') > -1:
            file_ew = folder + '/*2.D.*'
            file_ns = folder + '/*3.D.*'
        else:
            file_ew = folder + '/*E.D.*'
            file_ns = folder + '/*N.D.*'
        if i == 1:
            filename_in = glob.glob(file_ew)[0]
        if i == 2:
            filename_in = glob.glob(file_ns)[0]

        headers = {}

        # read file
        fh = open(filename_in, 'rt')
        for j in range(64):
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
        except ValueError:
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

        time = []
        for j in range(0, header['npts']):
            t = j * header['delta']
            time.append(t)

        if i == 1:
            inp_acc1 = np.asarray(acc_data) / 981  # in g
            # comp1=header['channel']
            npts1 = header['npts']
            time1 = time
        if i == 2:
            inp_acc2 = np.asarray(acc_data) / 981  # in g
            # comp2=header['channel']
            npts2 = header['npts']
            time2 = time

    return time1, time2, inp_acc1, inp_acc2, npts1, npts2


def to_utc_date_time(value):
    """
    """
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
    except ValueError:
        return None
    return x


def strtoint(sf):
    """

    """
    try:
        x = int(sf)
    except ValueError:
        return None
    return x


def create_nga_acc(num_rec, path_nga_folder):

    """
    All NGA-West2 recordings have to be stored in advance in a folder and
    renamed as::

        - RSN<NUM>_1.AT2 (1st horiz comp)
        - RSN<NUM>_2.AT2 (2nd horiz comp)
        - RSN<NUM>_3.AT2 (vertical component, not used)

    where `<NUM>` is the record sequence number of NGA recordings
    """

    # Import libraries
    import numpy as np

    # desc1 = ""
    # desc2 = ""
    time1 = []
    time2 = []
    inp_acc1 = []
    inp_acc2 = []
    npts1 = []
    npts2 = []
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
            counter = counter + 1
        
        if i == 1:
            inp_acc1 = np.asarray(acc_data)
        if i == 2:
            inp_acc2 = np.asarray(acc_data)

    return time1, time2, inp_acc1, inp_acc2, npts1, npts2

def create_kiknet_acc(recid, path_kiknet_folder, 
        fminNS2, fmaxNS2, fminEW2, fmaxEW2):

    """
    KiK-net acc are stored within Database_small.hdf5 file
    """

    # Import libraries
    import numpy as np
    from obspy.core import Trace,UTCDateTime
    import re
    from obspy.signal import filter

    # desc1 = ""
    # desc2 = ""
    time1 = []
    time2 = []
    inp_acc1 = []
    inp_acc2 = []
    npts1 = []
    npts2 = []
    for i in range(1, 3):
        if i==1:
            comp = 'EW2'
            fmin = fminEW2
            fmax = fmaxEW2
        elif i==2:
            comp = 'NS2'
            fmin = fminNS2
            fmax = fmaxNS2

        file_acc = path_kiknet_folder + '/' + str(recid) + '/' + str(recid) + '.' + comp
        hdrnames = ['Origin Time', 'Lat.', 'Long.', 'Depth. (km)', 'Mag.',
                    'Station Code', 'Station Lat.', 'Station Long.',
                    'Station Height(m)', 'Record Time', 'Sampling Freq(Hz)',
                    'Duration Time(s)', 'Dir.', 'Scale Factor', 'Max. Acc. (gal)',
                    'Last Correction', 'Memo.']
        acc_data = []
        time = []
        with open(file_acc, 'r') as f:
            content = f.readlines()
        counter = 0
        for line in content:
            if counter < 17:
                if not line.startswith(hdrnames[counter]):
                    sys.exit("Expected line to start with %s but got %s "
                                % (hdrnames[counter], line))
                else:
                    flds=line.split()

            if(counter==0):
                origin_time = flds[2] + ' ' + flds[3]
                origin_time = UTCDateTime.strptime(origin_time, '%Y/%m/%d %H:%M:%S')
                # All times are in Japanese standard time which is 9 hours ahead of UTC
                origin_time -= 9 * 3600.

            elif(counter==1):
                lat = float(flds[1])
            
            elif(counter==2):
                lon = float(flds[1])

            elif(counter==3):
                dp = float(flds[2])

            elif(counter==4):
                mag = float(flds[1])

            elif(counter==5):
                stnm = flds[2]

            elif(counter==6):
                stla = float(flds[2])

            elif(counter==7):
                stlo = float(flds[2])

            elif(counter==8):
                stel = float(flds[2])

            elif(counter==9):
                record_time = flds[2] + ' ' + flds[3]
                # A 15 s delay is added to the record time by the
                # the K-NET and KiK-Net data logger
                record_time = UTCDateTime.strptime(record_time, '%Y/%m/%d %H:%M:%S') - 15.0
                # All times are in Japanese standard time which is 9 hours ahead of UTC
                record_time -= 9 * 3600.

            elif(counter==10):
                freqstr = flds[2]
                m = re.search('[0-9]*', freqstr)
                freq = int(m.group())

            elif(counter==11):
                duration = float(flds[2])

            elif(counter==12):
                channel = flds[1].replace('-', '')
                kiknetcomps = {'1': 'NS1', '2': 'EW1', '3': 'UD1',
                       '4': 'NS2', '5': 'EW2', '6': 'UD2'}
                if channel.strip() in kiknetcomps.keys():  # kiknet directions are 1-6
                    channel = kiknetcomps[channel.strip()]

            elif(counter==13):
                eqn = flds[2]
                num, denom = eqn.split('/')
                num = float(re.search('[0-9]*', num).group())
                denom = float(denom)
                # convert the calibration from gal to m/s^2
                calib = 0.01 * num / denom

            elif(counter==14):
                accmax = float(flds[3])

            elif(counter==15):
                last_correction=flds[2] + ' ' + flds[3]
                last_correction = UTCDateTime.strptime(last_correction, '%Y/%m/%d %H:%M:%S')
                # All times are in Japanese standard time which is 9 hours ahead of UTC
                last_correction -= 9 * 3600.

            elif counter > 16:
                data = str(line).split()
                for value in data:
                    a = float(value)
                    acc_data.append(a)
            counter = counter + 1

        data = np.array(acc_data)
        tr = Trace(data)
        tr.detrend("linear")
        tr.taper(max_percentage=0.05, type='cosine', side='both')
        filter_order=4
        pad = np.zeros(int(round(1.5*filter_order/fmin*freq)))
        tr.data = np.concatenate([pad, tr.data, pad])
        fN = freq/2
        if fmax < fN:
            tr.data=filter.bandpass(tr.data,freqmin=fmin, freqmax=fmax, df=freq, corners=4, zerophase=True )
        else:
            tr.data=filter.highpass(tr.data,freq=fmin, df=freq, corners=4, zerophase=True )
        tr.data = tr.data[len(pad):len(tr.data)-len(pad)]
        tr.data = tr.data * calib / 9.81 #in g

        npts=len(tr.data)

        time = []
        for j in range(0, npts):
            t = j * 1/freq
            time.append(t)
        time=np.asarray(time)
        if i == 1:
            inp_acc1 = tr.data
            npts1 = npts
            time1 = time
        if i == 2:
            inp_acc2 = tr.data
            npts2 = npts
            time2 = time
    return time1, time2, inp_acc1, inp_acc2, npts1, npts2
