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

def scale_acc(n_gm, recid, path_nga, path_esm, source, event, station,
              output_folder, sf, comp, path_kiknet, 
              fminNS2, fmaxNS2, fminEW2, fmaxEW2, vertical_component):

    """

    Scales selected accelerograms and writes them in the output folder.

    :code:`nGM` x 2 files are created. Their name is::

        GMR_time_scaled_acc_<GMnum>_<comp>.txt

    where:
        - <`GMnum`> is a sequential number ranging from 1 to :code:`nGM`
        - <`comp`> can be `1` or `2` and indicates the horizontal component of motion

    Each file contains the selected scaled accelerograms, expressed with
    2 columns: time series `(s)` and accelerations `(g)`

    """
    # Import libraries
    import numpy as np
    from .create_acc import create_nga_acc
    from .create_acc import create_esm_acc
    from .create_acc import create_kiknet_acc


    # Read accelerograms, save them and apply scaling factor
    for i in np.arange(n_gm):
        time1 = []
        time2 = []
        time3 = []
        inp_acc1 = []
        inp_acc2 = []
        inp_acc3 = []
        npts1 = 0
        npts2 = 0
        npts3 = 0

        if source[i] == 'NGA-West2':
            val = int(recid[i])
            [time1, time2, time3, inp_acc1, inp_acc2, inp_acc3, npts1, 
                    npts2, npts3] = \
            create_nga_acc(val, path_nga)

        elif source[i] == 'ESM':
            folder_esm = path_esm + '/' + event[i] + '-' + station[i]
            [time1, time2, time3, inp_acc1, inp_acc2, inp_acc3, npts1, 
                    npts2, npts3] = \
            create_esm_acc(folder_esm,event[i],station[i],i)

        if source[i] == 'KiK-net':
            val = str(recid[i])
            [time1, time2, time3, inp_acc1, inp_acc2, inp_acc3, npts1, 
                    npts2, npts3] = \
            create_kiknet_acc(val, path_kiknet, 
                    fminNS2[i], fmaxNS2[i], fminEW2[i], fmaxEW2[i])

        if(comp[i]=='-'):
        # Create the filenames
            file_time_scaled_acc_out_1 = (output_folder +
                                      '/GMR_time_scaled_acc_' + str(i + 1) +
                                      '_1.txt')
            file_time_scaled_acc_out_2 = (output_folder +
                                      '/GMR_time_scaled_acc_' + str(i + 1) +
                                      '_2.txt')

            with open(file_time_scaled_acc_out_1, "wb") as f1:
                for j in np.arange(npts1):
                    f1.write("{:10.3f} {:15.10f}\n".format(time1[j],
                                                       inp_acc1[j] * sf[i]))
            with open(file_time_scaled_acc_out_2, "wb") as f2:
                for j in np.arange(npts2):
                    f2.write("{:10.3f} {:15.10f}\n".format(time2[j],
                                                       inp_acc2[j] * sf[i]))
        else:
            file_time_scaled_acc_out = (output_folder + 
                                      '/GMR_time_scaled_acc_' + str(i + 1) +
                                      '.txt')
            if(comp[i]==1):
                time=time1
                inp_acc=inp_acc1
                npts=npts1
            elif(comp[i]==2):
                time=time2
                inp_acc=inp_acc2
                npts=npts2
            with open(file_time_scaled_acc_out, "wb") as f1:
                for j in np.arange(npts):
                    f1.write("{:10.3f} {:15.10f}\n".format(time[j],
                                                       inp_acc[j] * sf[i]))

        
        if vertical_component == 'yes':
            time=time3
            inp_acc=inp_acc3
            npts=npts3
            file_time_scaled_vert_acc_out = (output_folder + 
                                      '/GMR_time_scaled_vertical_acc_' + str(i + 1) +
                                      '.txt')
            with open(file_time_scaled_vert_acc_out, "wb") as f1:
                for j in np.arange(npts):
                    f1.write("{:10.3f} {:15.10f}\n".format(time[j],
                                                       inp_acc[j] * sf[i]))

