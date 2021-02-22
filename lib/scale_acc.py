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

def scale_acc(n_gm, nga, path_nga, path_esm, source, event, station, name,
              output_folder, sf):
    """
    Scales selected accelerograms and writes selected accelerograms
    """
    # Import libraries
    import numpy as np
    import os
    from lib.create_acc import create_nga_acc
    from lib.create_acc import create_esm_acc
    from zipfile import ZipFile
    import requests

    # Read accelerograms, save them and apply scaling factor
    for i in np.arange(n_gm):
        time1 = []
        time2 = []
        inp_acc1 = []
        inp_acc2 = []
        npts1 = 0
        npts2 = 0
        if source[i] == 'NGA-West2':
            val = int(nga[i])
            [time1, time2, inp_acc1, inp_acc2, npts1, npts2] = \
                create_nga_acc(val, path_nga)
        elif source[i] == 'ESM':
            folder_esm = path_esm + '/' + event[i] + '-' + station[i]
            if not os.path.isdir(folder_esm):
                zip_output = 'output_' + str(i) + '.zip'
                
                params = (
                    ('eventid', event[i]),
                    ('data-type', 'ACC'),
                    ('station', station[i]),
                    ('format', 'ascii'),
                )

                files = {
                    'message': ('path/to/token.txt', open('token.txt', 'rb')),
                }

                headers={'Authorization': 'token {}'.format('token.txt')}

                url = 'https://esm-db.eu/esmws/eventdata/1/query'

                req = requests.post(url=url, params=params, files=files)

                if req.status_code == 200:
                    with open(zip_output, "wb") as zf:
                        zf.write(req.content)
                else:
                    sys.exit()

                with ZipFile(zip_output, 'r') as zipObj:
                    zipObj.extractall(folder_esm)
                os.remove(zip_output)
            [time1, time2, inp_acc1, inp_acc2, npts1, npts2] = \
                create_esm_acc(folder_esm)

        # Create the filenames
        file_time_scaled_acc_out_1 = (output_folder + '/' + name +
                                      '/GMR_time_scaled_acc_' + str(i + 1) +
                                      '_1.txt')
        file_time_scaled_acc_out_2 = (output_folder + '/' + name +
                                      '/GMR_time_scaled_acc_' + str(i + 1) +
                                      '_2.txt')

        with open(file_time_scaled_acc_out_1, "w", newline='') as f1:
            for j in np.arange(npts1):
                f1.write("{:10.3f} {:15.10f}\n".format(time1[j],
                                                       inp_acc1[j] * sf[i]))
        with open(file_time_scaled_acc_out_2, "w", newline='') as f2:
            for j in np.arange(npts2):
                f2.write("{:10.3f} {:15.10f}\n".format(time2[j],
                                                       inp_acc2[j] * sf[i]))
