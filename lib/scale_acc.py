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

def scale_acc(nGM,NGA,path_NGA,path_ESM,source,event,station,name,output_folder,SF):
     # Import libraries
    import numpy as np
    import os, sys
    from lib.create_acc import create_NGA_acc
    from lib.create_acc import create_ESM_acc

    #Read accelegrams, save them and apply scaling factor
    dts = []
    durs =[]
    names1 =[]
    names2 =[]
    for i in np.arange(nGM):
        time1=[]
        time2=[]
        inp_acc1=[]
        inp_acc2=[]
        npts1=0
        npts2=0
        comp1=''
        comp2=''
        desc1=''
        desc2=''
        if(source[i]=='NGA-West2'):
            val=int(NGA[i])
            [desc1,desc2,time1,time2,inp_acc1,inp_acc2,npts1,npts2]=create_NGA_acc(val,path_NGA)
            desc1='%'+desc1
            desc2='%'+desc2
        elif(source[i]=='ESM'):
            folder_ESM=path_ESM+'/'+event[i]+'-'+station[i]
            if not os.path.isdir(folder_ESM):
                zip_output='output_'+str(i)+'.zip'
                command='curl -X POST -F "message=@token.txt" "https://esm-db.eu/esmws/eventdata/1/query?eventid='+event[i]+'&data-type=ACC&station='+station[i]+'&format=ascii" -o '+zip_output
                os.system(command)
                command='unzip -o '+zip_output+' -d '+folder_ESM
                os.system(command)
                command='rm '+zip_output
                os.system(command)
            [time1,time2,inp_acc1,inp_acc2,npts1,npts2,comp1,comp2]=create_ESM_acc(folder_ESM)
            desc1='%'+event[i]+' '+station[i]+' '+comp1
            desc2='%'+event[i]+' '+station[i]+' '+comp2

        # Get the time steps and durations
        #dts.append(time1[1]-time1[0])
        #durs.append(time1[-1])

        # Create the filenames
        file_time_scaled_acc_out_1=output_folder+'/'+name+'/GMR_time_scaled_acc_'+str(i+1)+'_1.txt'
        file_time_scaled_acc_out_2=output_folder+'/'+name+'/GMR_time_scaled_acc_'+str(i+1)+'_2.txt'

        with open(file_time_scaled_acc_out_1, "w",newline='') as f1:
            for j in np.arange(npts1):
                f1.write("{:10.3f} {:15.10f}\n".format(time1[j],inp_acc1[j]*SF[i]))
        with open(file_time_scaled_acc_out_2, "w",newline='') as f2:
            for j in np.arange(npts2):
                f2.write("{:10.3f} {:15.10f}\n".format(time2[j],inp_acc2[j]*SF[i]))

        #file_scaled_acc_out_1=output_folder+'/'+name+'/GMR_scaled_acc_'+str(i+1)+'_1.txt'
        #file_scaled_acc_out_2=output_folder+'/'+name+'/GMR_scaled_acc_'+str(i+1)+'_2.txt'
        #names1.append('GMR_scaled_acc_'+str(i+1)+'_1.txt')
        #names2.append('GMR_scaled_acc_'+str(i+1)+'_2.txt')

        #with open(file_scaled_acc_out_1, "w",newline='') as f1:
        #    for j in np.arange(npts1):
        #        f1.write("{:15.10f}\n".format(inp_acc1[j]*SF[i]))
        #f1.close()
        #with open(file_scaled_acc_out_2, "w",newline='') as f2:
        #    for j in np.arange(npts2):
        #        f2.write("{:15.10f}\n".format(inp_acc2[j]*SF[i]))
        #f2.close()


#        # Print the time steps and the durations also
#        file_dts=output_folder+'/'+name+'/GMR_dts.txt'
#        file_durs=output_folder+'/'+name+'/GMR_durs.txt'
#        file_names1=output_folder+'/'+name+'/GMR_names1.txt'
#        file_names2=output_folder+'/'+name+'/GMR_names2.txt'
#
#        with open(file_dts, "w",newline='') as f1:
#            for j in np.arange(len(dts)):
#                f1.write("{:15.10f}\n".format(dts[j]))
#        f1.close()
#        with open(file_durs, "w",newline='') as f2:
#            for j in np.arange(len(durs)):
#                f2.write("{:15.10f}\n".format(durs[j]))
#        f2.close()
#        with open(file_names1, "w",newline='') as f1:
#            for j in np.arange(len(names1)):
#                f1.write("{:s}\n".format(names1[j]))
#        f1.close()
#        with open(file_names2, "w",newline='') as f2:
#            for j in np.arange(len(names2)):
#                f2.write("{:s}\n".format(names2[j]))
#        f2.close()
