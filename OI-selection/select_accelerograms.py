from openquake.hazardlib.gsim.mgmpe import akkar_coeff_table as act
from openquake.hazardlib import gsim, imt, const
import numpy as np
import pandas as pd
import os, sys
import csv
import matplotlib.pyplot as plt
import matplotlib
import time
from scipy.stats import skew
import requests
import glob
from obspy.core import Stream, Trace, UTCDateTime, Stats
import re
import string


scale_accelerograms=0 #(0=no, 1=yes)

selection_type=0 #=1 uses AvgSA, =0 uses T1 (period at which spectra should be scaled and matched)
T1=0 #conditioning period. Used only if selection_type=0, for PGA uses T1=0
site_code=[1]
rlz_code=[1]

#PSHA and disaggregation parameters
path_results='../OQ-hazard/SA/outs'
num_disagg=36
num_classical=36
probability_of_exceedance=0.02
investigation_time=50

#parameters for conditional spectrum computation
TgtPer=[0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1,2,3,4]
corr_type='akkar'   #baker_jayaram or akkar
GMPE='Akkar_Bommer_2010' #for now only available Akkar_Bommer_2010 Chiou_Youngs_2014 Boore_et_al_2014
avg_periods=[0.2,0.3,0.4,0.5,0.75,1.0]
rake=0.
hypo_depth=10.
AR=1
Vs30=300

#database
database_path='database_flatfile_NGA2_only_filtered3.csv'
#screen for suitable ground motions
allowedRecs_Vs30=[180,360]     # upper and lower bound of allowable Vs30 values
allowedEC8code='C'
maxsf=5.0 #The maximum allowable scale factor
radius_dist=50 #km
radius_mag=0.25
maximum_depth=30

#parameters for statistically simulation of response spectra
nTrials = 5 #number of iterations of the initial spectral simulation step to perform
np.random.seed(333)
weights=[1.0,2.0,0.3] #[Weights for error in mean, standard deviation and skewness] Used to find the simulated spectra that best match the target from the statistically simulated response spectra
nGM=30 #number of records to select ==> number of records to select since the code search the database spectra most similar to each simulated spectrum

#otimization parameters. Execution of incremental changes to the initially selected ground motion set to further optimise its fit to the target spectrum distribution. 
nLoop=2 #Number of loops of optimization to perform
penalty=10 #>0 to penalize selected spectra moire than 3 sigma from the target at any period, =0 otherwise.

path_NGA_folder='NGA2records' #NGA recordings have to be stored
output_ESM='ESM'
#If not, found in the folder, ESM recording are authomatically downloaded from internet, need to generate the file token.txt
#At first you need to register at: http://tex.mi.ingv.it/
#curl -X POST -F 'message={"user_email": "elisa.zuccolo@eucentre.it","user_password": "password"}' "http://tex.mi.ingv.it/esmws/generate-signed-message/1/query" > token.txt

# -----------------------------------------------------------------------------
# Compute conditional hazard spectrum

def akkar_correlation(t1, t2):
    """
    Read the period-dependent correlation coefficient matrix rho_H(T0,T) given in 
    Akkar S, Sandikkaya MA, Ay BO (2014) Compatible ground-motion prediction equations
    for damping scaling factors and vertical-to-horizontal spectral amplitude ratios
    for the broader Europe region, Bull Earthquake Eng 12:517-547
    """

    return act.coeff_table[act.periods.index(t1)][act.periods.index(t2)]


def baker_jayaram_correlation(t1, t2):
    """
    NOTE: subroutine taken from: https://usgs.github.io/shakemap/shakelib

    Produce inter-period correlation for any two spectral periods.

    Based upon:
    Baker, J.W. and Jayaram, N., "Correlation of spectral acceleration
    values from NGA ground motion models," Earthquake Spectra, (2007).

    Args:
        t1, t2 (float):
            The two periods of interest.

    Returns:
        rho (float): The predicted correlation coefficient

    """
    t_min = min(t1, t2)
    t_max = max(t1, t2)

    c1 = 1.0 - np.cos(np.pi / 2.0 - np.log(t_max / max(t_min, 0.109)) * 0.366)

    if t_max < 0.2:
        c2 = 1.0 - 0.105 * (1.0 - 1.0 / (1.0 + np.exp(100.0 * t_max - 5.0)))*(t_max - t_min) / (t_max - 0.0099)
    else:
        c2 = 0

    if t_max < 0.109:
        c3 = c2
    else:
        c3 = c1

    c4 = c1 + 0.5 * (np.sqrt(c3) - c3) * (1.0 + np.cos(np.pi * t_min / 0.109))

    if t_max <= 0.109:
        rho = c2
    elif t_min > 0.109:
        rho = c1
    elif t_max < 0.2:
        rho = min(c2, c4)
    else:
        rho = c4

    return rho

def toUTCDateTime(value):
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
    try:
        x = float(sf)
    except:
        return None
    return x

def strtoint(sf):
    try:
        x = int(sf)
    except:
        return None
    return x

def create_ESM_acc(folder):
    for l in range(1,3):
        if(folder.find('ESM/GR') >-1):
            file_EW=folder+'/*2.D.*'
            file_NS=folder+'/*3.D.*'
        else:
            file_EW=folder+'/*E.D.*'
            file_NS=folder+'/*N.D.*'
        if(l==1):
            filename_in=glob.glob(file_EW)[0]
        if(l==2):
            filename_in=glob.glob(file_NS)[0]

        acc_data=[]
        time=[]

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
            header['starttime'] = toUTCDateTime(headers['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS']) # use toUTCDateTime to convert from DYNA format
        except:
            header['starttime'] = toUTCDateTime('19700101_000000')
        header['sampling_rate'] = 1/float(headers['SAMPLING_INTERVAL_S'])
        header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
        header['npts'] = int(headers['NDATA'])
        header['calib'] = 1 # not in file header

        ##DYNA dict float data
        header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(headers['EVENT_LATITUDE_DEGREE'])
        header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(headers['EVENT_LONGITUDE_DEGREE'])
        header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
        header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
        header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
        header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
        header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(headers['STATION_LATITUDE_DEGREE'])
        header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(headers['STATION_LONGITUDE_DEGREE'])
        header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
        header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(headers['EPICENTRAL_DISTANCE_KM'])
        header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
        header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
        header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(headers['INSTRUMENTAL_FREQUENCY_HZ'])
        header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(headers['INSTRUMENTAL_DAMPING'])
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

        header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(headers['LOW_CUT_FREQUENCY_HZ'])
        header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(headers['HIGH_CUT_FREQUENCY_HZ'])

        ##DYNA dict int data
        header['dyna']['STATION_ELEVATION_M'] = strtoint(headers['STATION_ELEVATION_M'])
        header['dyna']['SENSOR_DEPTH_M'] = strtoint(headers['SENSOR_DEPTH_M'])
        header['dyna']['N_BIT_DIGITAL_CONVERTER'] =  strtoint(headers['N_BIT_DIGITAL_CONVERTER'])
        header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

        ##DYNA dict string data
        header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
        header['dyna']['EVENT_ID'] = headers['EVENT_ID']
        header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
        header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
        header['dyna']['MAGNITUDE_W_REFERENCE'] = headers['MAGNITUDE_W_REFERENCE']
        header['dyna']['MAGNITUDE_L_REFERENCE'] = headers['MAGNITUDE_L_REFERENCE']
        header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
        header['dyna']['STATION_NAME'] = headers['STATION_NAME']
        header['dyna']['SITE_CLASSIFICATION_EC8'] = headers['SITE_CLASSIFICATION_EC8']
        header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers['MORPHOLOGIC_CLASSIFICATION']
        header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers['DATE_TIME_FIRST_SAMPLE_PRECISION']
        header['dyna']['UNITS'] = headers['UNITS']
        header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
        header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers['INSTRUMENT_ANALOG/DIGITAL']
        header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
        header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
        header['dyna']['LATE_NORMAL_TRIGGERED'] = headers['LATE/NORMAL_TRIGGERED']
        header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
        header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
        header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
        header['dyna']['PROCESSING'] = headers['PROCESSING']
        header['dyna']['DATA_LICENSE'] = headers['DATA_LICENSE']
        header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers['DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
        header['dyna']['DATA_CITATION'] = headers['DATA_CITATION']
        header['dyna']['DATA_CREATOR'] = headers['DATA_CREATOR']
        header['dyna']['ORIGINAL_DATA_MEDIATOR_CITATION'] = headers['ORIGINAL_DATA_MEDIATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_MEDIATOR'] = headers['ORIGINAL_DATA_MEDIATOR']
        header['dyna']['ORIGINAL_DATA_CREATOR_CITATION'] = headers['ORIGINAL_DATA_CREATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_CREATOR'] = headers['ORIGINAL_DATA_CREATOR']
        header['dyna']['USER1'] = headers['USER1']
        header['dyna']['USER2'] = headers['USER2']
        header['dyna']['USER3'] = headers['USER3']
        header['dyna']['USER4'] = headers['USER4']
        header['dyna']['USER5'] = headers['USER5']

        # read data
        acc_data = np.loadtxt(fh, dtype='float32')
        fh.close()

        for j in range (0,header['npts']):
            t = j * header['delta']
            time.append(t)

        if(l==1):
            inp_acc1=np.asarray(acc_data)/981 #in g
            comp1=header['channel']
            npts1=header['npts']
            time1=time
        if(l==2):
            inp_acc2=np.asarray(acc_data)/981 #in g
            comp2=header['channel']
            npts2=header['npts']
            time2=time

    return time1,time2,inp_acc1,inp_acc2,npts1,npts2,comp1,comp2


def create_NGA_acc(num_rec,path_NGA_folder):
    desc1=""
    desc2=""
    for i in range(1,3):
        file_acc=path_NGA_folder+'/RSN'+str(num_rec)+'_'+str(i)+'.AT2'
        acc_data=[]
        time=[]
        with open(file_acc,'r') as f:
            content = f.readlines()
        counter = 0
        row4Val=""
        for x in content:
            if counter == 1:
                if(i==1):
                    desc1 = x
                    desc1 = desc1[0:(len(x)-1)]
                if(i==2):
                    desc2 = x
                    desc2 = desc2[0:(len(x)-1)]
            elif counter == 3:
                row4Val = x
                val = row4Val.split()
                npts_comma =val[1]
                npts=int(npts_comma[0:len(npts_comma)-1])
                if(i==1):
                    npts1=npts
                if(i==2):
                    npts2=npts
                dt = float(val[3])
                for j in range (0,npts):
                    t = j * dt
                    time.append(t)
                if(i==1):
                    time1=time
                if(i==2):
                    time2=time
            elif counter > 3:
                data = str(x).split()
                for value in data:
                    a = float(value)
                    acc_data.append(a)
                    if(i==1):
                        inp_acc1=np.asarray(acc_data)
                    if(i==2):
                        inp_acc2=np.asarray(acc_data)
            counter = counter + 1
    return desc1,desc2,time1,time2,inp_acc1,inp_acc2,npts1,npts2

def screen_database(database_path,allowedRecs_Vs30,allowedRecs_Mag,allowedRecs_D,allowedEC8code,TgtPer,nGM,maximum_depth):
    dbacc=pd.read_csv(database_path,sep=';',engine='python')
    knownPer=np.array([0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,5,6,7,8,9,10])

    event_id=dbacc['event_id']
    event_mw=dbacc['Mw']
    event_mag=dbacc['M']
    record_sequence_number_NGA=dbacc['record_sequence_number_NGA']
    station_ec8=dbacc['ec8_code']
    station_vs30=dbacc['vs30_m_sec']
    acc_distance=dbacc['epi_dist']
    station_code=dbacc['station_code']
    event_depth=dbacc['ev_depth_km']
    sensor_depth=dbacc['sensor_depth_m']
    is_free_field_ESM=dbacc['proximity_code']
    is_free_field_NGA=dbacc['GMX_first']
    source=dbacc['source']
    epi_lon=dbacc['epi_lon']
    epi_lat=dbacc['epi_lat']

    # Match periods (known periods and target periods for error computations) 
    # save the indicies of the matched periods in knownPer
    indPer = np.zeros((len(TgtPer),1), dtype=int);
    for i in np.arange(len(TgtPer)):
        indPer[i]=np.argmin(np.absolute(knownPer-TgtPer[i]))

    # Remove any repeated values from TgtPer and redefine TgtPer as periods 
    # provided in databases
    indPer = np.unique(indPer);
    recPer = knownPer[indPer];
    
    SA_list=[]
    allowedIndex=[]
    for i in np.arange(len(event_id)):
        rotD50=np.array([dbacc['rotD50_pga'][i],dbacc['rotD50_T0_010'][i],dbacc['rotD50_T0_025'][i],dbacc['rotD50_T0_040'][i],dbacc['rotD50_T0_050'][i],dbacc['rotD50_T0_070'][i],dbacc['rotD50_T0_100'][i],dbacc['rotD50_T0_150'][i],dbacc['rotD50_T0_200'][i],dbacc['rotD50_T0_250'][i],dbacc['rotD50_T0_300'][i],dbacc['rotD50_T0_350'][i],dbacc['rotD50_T0_400'][i],dbacc['rotD50_T0_450'][i],dbacc['rotD50_T0_500'][i],dbacc['rotD50_T0_600'][i],dbacc['rotD50_T0_700'][i],dbacc['rotD50_T0_750'][i],dbacc['rotD50_T0_800'][i],dbacc['rotD50_T0_900'][i],dbacc['rotD50_T1_000'][i],dbacc['rotD50_T1_200'][i],dbacc['rotD50_T1_400'][i],dbacc['rotD50_T1_600'][i],dbacc['rotD50_T1_800'][i],dbacc['rotD50_T2_000'][i],dbacc['rotD50_T2_500'][i],dbacc['rotD50_T3_000'][i],dbacc['rotD50_T3_500'][i],dbacc['rotD50_T4_000'][i],dbacc['rotD50_T5_000'][i],dbacc['rotD50_T6_000'][i],dbacc['rotD50_T7_000'][i],dbacc['rotD50_T8_000'][i],dbacc['rotD50_T9_000'][i],dbacc['rotD50_T10_000'][i]])
        if(source[i]=='ESM'):
            SA_geo=rotD50/981 #in g
        if(source[i]=='NGA-West2'):
            SA_geo=rotD50 #already in g
        SA_list.append(SA_geo)
        if(all(v > 0 for v in SA_geo)):
            if(source[i]=='ESM'):
                if(event_depth[i]<=maximum_depth):
                    if(is_free_field_ESM[i]==0):
                        if(event_mw[i]>=allowedRecs_Mag[0] and event_mw[i]<=allowedRecs_Mag[1]):
                            if(acc_distance[i]>=allowedRecs_D[0] and acc_distance[i]<=allowedRecs_D[1]):
                                if(allowedEC8code=='None'):
                                    allowedIndex.append(i)
                                else:
                                    if(station_ec8[i]==allowedEC8code or station_ec8[i]==allowedEC8code+"*"):
                                        allowedIndex.append(i)
                                    elif(station_vs30[i]>=allowedRecs_Vs30[0] and station_vs30[i]<allowedRecs_Vs30[1]):
                                            allowedIndex.append(i)

            if(source[i]=='NGA-West2'):
                if(epi_lon[i]<-31 or epi_lon[i]>70):
                    if(event_depth[i]<=maximum_depth):
                        if(is_free_field_NGA[i]=="I"):
                            if(event_mag[i]>=allowedRecs_Mag[0] and event_mag[i]<=allowedRecs_Mag[1]):
                                if(acc_distance[i]>=allowedRecs_D[0] and acc_distance[i]<=allowedRecs_D[1]):
                                    if(station_vs30[i]>=allowedRecs_Vs30[0] and station_vs30[i]<allowedRecs_Vs30[1]):
                                        allowedIndex.append(i)
    SA=np.vstack(SA_list)
    SaKnown=SA[allowedIndex]

    # count number of allowed spectra
    nBig = len(allowedIndex)
    print(['Number of allowed ground motions = ', nBig])
    assert (nBig >= nGM), 'Warning: there are not enough allowable ground motions'

    return SaKnown,indPer,recPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance

def compute_avgSA(avg_periods,sctx, rctx, dctx):
    mean_list = []
    stddvs_list = []
    # Loop over averaging periods
    for period in avg_periods:
        # compute mean and standard deviation
        P=imt.SA(period=period)
        S=[const.StdDev.TOTAL]
        mean,std = bgmpe.get_mean_and_stddevs(sctx, rctx, dctx,P,S)
        mean_list.append(mean)
        stddvs_list.append(std[0]) # Support only for total!

    mean_avgsa = 0.
    stddvs_avgsa = 0.

    for i1 in np.arange(len(avg_periods)):
        mean_avgsa += mean_list[i1]
        for i2 in np.arange(len(avg_periods)):
            if(corr_type=='baker_jayaram'):
                rho = baker_jayaram_correlation(avg_periods[i1],avg_periods[i2])
            if(corr_type=='akkar'):
                rho = akkar_correlation(avg_periods[i1],avg_periods[i2])

            stddvs_avgsa += rho*stddvs_list[i1] * stddvs_list[i2]

    mean_avgsa *= (1./len(avg_periods))
    stddvs_avgsa *= (1./len(avg_periods))**2
    stddvs_avgsa=np.sqrt(stddvs_avgsa)
    return [np.exp(mean_avgsa),stddvs_avgsa]

def compute_rho_avgSA(per,avg_periods,sctx,rctx,dctx,stddvs_avgsa):
    sum_numeratore=0
    for i1 in avg_periods:
        if(corr_type=='baker_jayaram'):
            rho=baker_jayaram_correlation(per,i1)
        if(corr_type=='akkar'):
            rho=akkar_correlation(per,i1)
        S=[const.StdDev.TOTAL]
        mean1,std1 = bgmpe.get_mean_and_stddevs(sctx, rctx, dctx, imt.SA(period=i1),S)
        sum_numeratore=sum_numeratore+rho*std1[0]

    denominatore=len(avg_periods)*stddvs_avgsa
    rho_avgSA=sum_numeratore/denominatore
    return rho_avgSA

def simulate_spectra(nTrials,meanReq,covReq,stdevs,nGM):
# simulate response spectra from the target mean and covariance matrix
    devTotalSim =[]
    spettri=[]
    for j in np.arange(nTrials):
        SpectraSample=np.exp(np.random.multivariate_normal(meanReq,covReq,nGM))
        spettri.append(SpectraSample)
        # evaluate simulation
        sampleMeanErr=np.mean(np.log(SpectraSample),axis=0)-meanReq
        sampleStdErr=np.std(np.log(SpectraSample),axis=0)-stdevs
        sampleSkewnessErr=skew(np.log(SpectraSample),axis=0,bias=True)
        devTotalSim.append(weights[0]*sum(sampleMeanErr**2)+weights[1]**sum(sampleStdErr**2)+weights[2]*sum(sampleSkewnessErr**2))

    bestSample=devTotalSim.index(min(devTotalSim)) # find the simulated spectra that best match the target 
    simulated_spectra = spettri[bestSample] #return the best set of simulations

    return simulated_spectra


#****************************************************************************************************************

for ii in np.arange(len(site_code)):

    site = site_code[ii]
    rlz = rlz_code[ii]

    for poe in [probability_of_exceedance]:
        if(selection_type==1):
            disagg_results='rlz-'+str(rlz)+'-AvgSA-sid-'+str(site)+'-poe-0_Mag_Dist_'+str(num_disagg)+'.csv'
            name='AvgSA-site_'+str(site)+'-poe-'+str(poe)
        elif(selection_type==0):
            if T1 == 0: #need to handle T1=PGA
                disagg_results='rlz-'+str(rlz)+'-PGA-sid-'+str(site)+'-poe-0_Mag_Dist_'+str(num_disagg)+'.csv'
                name='PGA-site_'+str(site)+'-poe-'+str(poe)
            else:
                disagg_results='rlz-'+str(rlz)+'-SA('+Tstar+')-sid-'+str(site)+'-poe-0_Mag_Dist_'+str(num_disagg)+'.csv'
                name='SA('+Tstar+')-site_'+str(site)+'-poe-'+str(poe)
            print('This needs to be fixed, dont know why. maybe to automate the filename inputs?')
        else: #TODO
            print('TODO')
        print(name)

        file_with_OQ_acc_value='hazard_map-mean_'+str(num_classical)+'.csv'

        if(selection_type==1):
            selected_column='AvgSA-'+str(probability_of_exceedance)
        else: #TODO
            print('TODO')

        print(selected_column)

        #Retrieve disaggregation results

        meanLst = [],[] 
        df=pd.read_csv(''.join([path_results,'/',disagg_results]),skiprows=1) 
        df['rate'] = -np.log(1-df['poe'])/investigation_time 
        df['rate_norm'] = df['rate']/ df['rate'].sum() 
        mode=df.sort_values(by='rate_norm',ascending=False)[0:1] 
        meanMag=np.sum(df['mag']*df['rate_norm']) 
        meanDist=np.sum(df['dist']*df['rate_norm']) 
        print(meanMag,meanDist)
        allowedRecs_D=[meanDist-radius_dist,meanDist+radius_dist]
        allowedRecs_Mag=[meanMag-radius_mag,meanMag+radius_mag]

        #Retrieve conditioning value
        df=pd.read_csv(''.join([path_results,'/',file_with_OQ_acc_value]),skiprows=1)
        output_oq=df[selected_column]
        output_oq=output_oq[0]

# -----------------------------------------------------------------------------
        # Initialise GSIMs

        _ = gsim.get_available_gsims()

        if(GMPE=='Akkar_Bommer_2010'):
            bgmpe = gsim.akkar_bommer_2010.AkkarBommer2010()
        if(GMPE=='Chiou_Youngs_2014'):
            bgmpe = gsim.chiou_youngs_2014.ChiouYoungs2014()
        if(GMPE=='Boore_et_al_2014'):
            bgmpe = gsim.boore_2014.BooreEtAl2014()
        sctx = gsim.base.SitesContext()
        rctx = gsim.base.RuptureContext()
        dctx = gsim.base.DistancesContext()

# -----------------------------------------------------------------------------
        # Initialise contexts
        if(GMPE=='Chiou_Youngs_2014'):
            area=WC1994().get_median_area(meanMag,rake)
            length=np.sqrt(area*AR)
            width=area/length
            source_vertical_width=width*np.sin(np.radians(dip))
            ztor=max(hypo_depth-source_vertical_width/2.,upper_sd)
            if((ztor+source_vertical_width)>lower_sd):
                source_vertical_width=lower_sd-ztor
                width=source_vertical_width/np.sin(np.radians(dip))

            rjb=meanDist

            if(rjb==0):
                rx=0.5*width*np.cos(np.radians(dip))
            else:
                if(dip==90):
                    rx=rjb*np.sin(np.radians(azimuth))
                else:
                    if(azimuth>=0 and azimuth<90):
                        if(rjb*np.abs(np.tan(np.radians(azimuth)))<=width*np.cos(np.radians(dip))):
                            rx=rjb*np.abs(np.tan(np.radians(azimuth)))
                        else:
                            rx=rjb*np.tan(np.radians(azimuth))*np.cos(np.radians(azimuth)-np.arcsin(width*np.cos(np.radians(dip))*np.cos(np.radians(azimuth))/rjb))

            if(dip==90):
                rrup=np.sqrt(np.square(rjb)+np.square(ztor))
            else:
                if(rx<ztor*np.tan(np.radians(dip))):
                    rrup1=np.sqrt(np.square(rx)+np.square(ztor))
                if(rx>=ztor*np.tan(np.radians(dip)) and rx<=ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
                    rrup1=rx*np.sin(np.radians(dip))+ztor*np.cos(np.radians(dip))
                if(rx>ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
                    rrup1=np.sqrt(np.square(rx-width*np.cos(np.radians(dip)))+np.square(ztor+width*np.sin(np.radians(dip))))
                if(azimuth==90 or azimuth==-90):
                    ry=0
                elif(azimuth==0 or azimuth==180 or azimuth==-180):
                    ry=rjb
                else:
                    ry=np.abs(rx*1./np.tan(np.radians(azimuth)))
                rrup=np.sqrt(np.square(rrup1)+np.square(ry))

            print(rjb,rx,rrup)

        Dist = np.arange(meanDist,meanDist+1,1.)
        if(GMPE=='Chiou_Youngs_2014'):
            rx=np.arange(rx,rx+1,1.)
            rrup=np.arange(rrup,rrup+1,1.)
            setattr(rctx, 'ztor', ztor)
            setattr(rctx, 'dip', dip)
            setattr(dctx, 'rx', rx)
            setattr(dctx, 'rrup', rrup)
            setattr(sctx, 'vs30measured',0)
            setattr(sctx, 'z1pt0', z1pt0)
        Vs30 = Vs30+np.zeros(Dist.shape)
        setattr(rctx, 'mag', meanMag)
        setattr(rctx, 'rake', rake)
        setattr(rctx, 'hypo_depth', hypo_depth)
        setattr(dctx, 'rjb', Dist)
        setattr(sctx, 'vs30', Vs30)

        [SaKnown,indPer,TgtPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance]=screen_database(database_path,allowedRecs_Vs30,allowedRecs_Mag,allowedRecs_D,allowedEC8code,TgtPer,nGM,maximum_depth)

        TgtMean=[]
        rho=[]
        mean=[]
        sigma=[]
        if(selection_type==1):
            [mean_SaTcond,stddvs_SaTcond]=compute_avgSA(avg_periods,sctx, rctx, dctx)
            epsilon=(np.log(output_oq)-np.log(mean_SaTcond))/stddvs_SaTcond
        if(selection_type==0):
            if(T1==0):
                P = imt.PGA()
            else:
                P = imt.SA(period=T1)
            S=[const.StdDev.TOTAL]
            mean_SaTcond,stddvs_SaTcond=bgmpe.get_mean_and_stddevs(sctx,rctx,dctx,P,S)
            stddvs_SaTcond=stddvs_SaTcond[0]
            epsilon=(np.log(output_oq)-mean_SaTcond)/stddvs_SaTcond

        for per in TgtPer:
            if(per==0):
                P = imt.PGA()
            else:
                P=imt.SA(period=per)
            S=[const.StdDev.TOTAL]
            bMean_SA, bStDev_SA = bgmpe.get_mean_and_stddevs(sctx, rctx, dctx, P, S)
            mean.append(bMean_SA)
            sigma.append(bStDev_SA[0])
            if(selection_type==1):
                rho_per=compute_rho_avgSA(per,avg_periods,sctx,rctx,dctx,stddvs_SaTcond)
                rho.append(rho_per[0])
            if(selection_type==0):
                if(corr_type=='baker_jayaram'):
                    rho_per = baker_jayaram_correlation(per,T1)
                if(corr_type=='akkar'):
                    rho_per = akkar_correlation(per,T1)
                rho.append(rho_per)
            spectrum=bMean_SA+rho_per*bStDev_SA[0]*epsilon
            # (Log) Response Spectrum Mean: TgtMean
            TgtMean.append(spectrum[0])

        TgtMean=np.array(TgtMean)

#### Compute covariances and correlations at all periods
        TgtCovs = np.zeros((len(TgtPer),len(TgtPer)))
        for i in np.arange(len(TgtPer)):
            for j in np.arange(len(TgtPer)):
                Ti = TgtPer[i]
                Tj = TgtPer[j]
                varT = stddvs_SaTcond**2
                varT=varT[0]
                var1 = sigma[i]**2
                var1=var1[0]
                var2 = sigma[j]**2
                var2=var2[0]
        
                if(corr_type=='baker_jayaram'):
                    sigmaCorr=baker_jayaram_correlation(Ti, Tj)*np.sqrt(var1*var2)
                if(corr_type=='akkar'):
                    sigmaCorr=akkar_correlation(Ti, Tj)*np.sqrt(var1*var2)
                sigma11 = np.matrix([[var1, sigmaCorr], [sigmaCorr, var2]])
                sigma22 = np.array(varT)
                sigma12 =np.array([[rho[i]*np.sqrt(var1*varT)],[rho[j]*np.sqrt(var2*varT)]])
                sigmaCond = sigma11 - sigma12*1./(sigma22)*sigma12.T
                TgtCovs[i,j] = sigmaCond.item((0,1))

# find covariance values of zero and set them to a small number so that random number generation can be performed
        TgtCovs[np.absolute(TgtCovs) < 1e-10] = 1e-10

        meanReq=TgtMean
        covReq=TgtCovs
        stdevs=np.sqrt(np.diagonal(covReq))

        simulated_spectra=simulate_spectra(nTrials,meanReq,covReq,stdevs,nGM)

        sampleBig =np.log(SaKnown[:,indPer])

        if(selection_type==1):
            id_avgSA=[]
            id_avgSA_bool=np.isin(TgtPer,avg_periods) 
            for i in np.arange(len(TgtPer)):
                if(id_avgSA_bool[i] == True):
                    id_avgSA.append(i)
            id_avgSA=np.array(id_avgSA)
            lnSa1=np.mean(meanReq[id_avgSA])
        if(selection_type==0):
            id_T1=[]
            id_T1=np.where(TgtPer==T1)
            lnSa1=np.mean(meanReq[id_T1])

        recID = np.zeros(nGM,dtype=int)
        sampleSmall = []
        IMScaleFac = np.ones(nGM)
# Find database spectra most similar to each simulated spectrum
        for i in np.arange(nGM): # for each simulated spectrum
            err=np.zeros(nBig)*1000000 #initialize error matrix
            scaleFac = np.ones(nBig) # initialize scale factors to 1
            #compute scale factors and errors for each candidate ground motion
            for j in np.arange(nBig):
                if(selection_type==1):
                    rec_value=np.exp(sum(sampleBig[j,id_avgSA])/len(id_avgSA))
                if(selection_type==0):
                    rec_value=np.exp(sampleBig[j,id_T1])
                if (rec_value == 0):
                    scaleFac[j] = 1000000
                else:
                    scaleFac[j]=np.exp(lnSa1)/rec_value
                err[j] = sum((np.log(np.exp(sampleBig[j,:])*scaleFac[j]) - np.log(simulated_spectra[i,:]))**2)

            err[recID[0:i-1]] = 1000000 # exclude previously-selected ground motions
            err[scaleFac > maxsf] = 1000000 # exclude ground motions requiring too large of a scale factor
            err[scaleFac < 1./maxsf] = 1000000 # exclude ground motions requiring too large of a scale factor

            # find minimum-error ground motion
            recID[i]=np.argmin(err)
            min_err=np.min(err)
            assert (min_err < 1000), 'Warning: problem with simulated spectrum. No good matches found'
            IMScaleFac[i] = scaleFac[recID[i]] #store scale factor
            sampleSmall.append(np.log(np.exp(sampleBig[recID[i],:])*scaleFac[recID[i]])) #store scaled log spectrum

        sampleSmall=np.array(sampleSmall)

# Further optimize the ground motion selection

        print('Please wait...This algorithm takes a few minutes depending on the number of records to be selected')

        for k in np.arange(nLoop):
            for i in np.arange(nGM): # consider replacing each ground motion in the selected set
                minDev=100000

                sampleSmall=np.delete(sampleSmall,i, 0)
                recID=np.delete(recID,i)

                # Try to add a new spectrum to the subset list
                for j in np.arange(nBig):
                    if(selection_type==1):
                        rec_value=np.exp(sum(sampleBig[j,id_avgSA])/len(id_avgSA))
                    if(selection_type==0):
                        rec_value=np.exp(sampleBig[j,id_T1])
                    if (rec_value == 0):
                        scaleFac[j] = 1000000
                    else:
                        scaleFac[j] = np.exp(lnSa1)/rec_value
                    added1=np.reshape((sampleBig[j,:]+np.log(scaleFac[j])),(1,len(TgtPer)))
                    sampleSmall=np.concatenate((sampleSmall,added1)) #add candidate to set
                    # Compute deviations from target
                    devMean=np.mean(sampleSmall,axis=0) - meanReq
                    devSkew=skew(sampleSmall,axis=0,bias=True)
                    devSig=np.std(sampleSmall,axis=0) - stdevs
                    devTotal = weights[0]*sum(devMean**2)+weights[1]*sum(devSig**2)

                    # Penalize bad spectra (set penalty to zero if this is not required)
                    if(penalty != 0):
                        for m in np.arange(len(sampleSmall)):
                            devTotal = devTotal + sum(np.absolute(np.exp(sampleSmall[m,:])>np.exp(meanReq+3*stdevs)))*penalty

                    if(scaleFac[j] > maxsf or scaleFac[j] < 1./maxsf):
                        devTotal = devTotal + 1000000

                    # Should cause improvement and record should not be repeated
                    if (devTotal < minDev and not any(recID==j)):
                        minID = np.zeros(1,dtype=int)
                        minID[0]=j
                        minDev = devTotal
                    end=len(sampleSmall)
                    sampleSmall = sampleSmall[0:end-1,:]
        
                # Add new element in the right slot
                IMScaleFac[i] = scaleFac[minID]
                end=len(sampleSmall)
                added2=np.reshape((sampleBig[minID,:]+np.log(scaleFac[minID])),(1,len(TgtPer)))
                if(i>0):
                    sampleSmall = np.concatenate((sampleSmall[0:i,:],added2,sampleSmall[i:end,:]))
                    recID = np.concatenate((recID[0:i],minID,recID[i:end]))
                else:
                    sampleSmall = np.concatenate((added2,sampleSmall[i:end,:]))
                    recID = np.concatenate((minID,recID[i:end]))

        # Output information
        finalRecords = recID
        recIdx = [allowedIndex[i] for i in finalRecords]
        finalScaleFactors = IMScaleFac
        meanrecorded=np.mean(np.exp(sampleSmall),axis=0)

        for i in np.arange(nGM):
            plt.loglog(TgtPer,np.exp(sampleSmall[i,:]),'g')
        plt.loglog(TgtPer,np.exp(meanReq),'r',label='target')
        plt.loglog(TgtPer,np.exp(meanReq+2*stdevs),'--r',label='target+2*sigma')
        plt.loglog(TgtPer,np.exp(meanReq-2*stdevs),'--r',label='target-2*sigma')
        plt.loglog(TgtPer,meanrecorded,'b',label='mean recorded')
        plt.xlabel('Period (s)')
        plt.ylabel('Spectral acceleration (g)')
        plt.xlim((0.01, np.max(TgtPer)))
        plt.legend()

        plt.grid(True)

        if not os.path.exists(name):
            os.makedirs(name)

        name_fig=name+'/'+name+'_selection.png'
        plt.savefig(name_fig, bbox_inches='tight')

        plt.close()

        # Output results to a text file
        blank='-'
        name_summary=name+'/'+name+"_summary_selection.txt"
        with open(name_summary, "w") as f:
            f.write("{} {}\n".format('reference hazard value = ',output_oq))
            f.write("{} {}\n".format('mean_mag_disag = ',meanMag))
            f.write("{} {}\n".format('mean_dist_disag = ',meanDist))
            f.write("num source event_id_ESM station_code_ESM recID_NGA magnitude distance scale_factor\n")
            for i in np.arange(nGM):
                elemento=recIdx[i]
                if(source[elemento]=='ESM'):
                    f.write("{} {} {} {} {} {} {} {:6.2f} \n".format(i+1,source[elemento],event_id[elemento],station_code[elemento],blank,event_mw[elemento],acc_distance[elemento],finalScaleFactors[i]))
                if(source[elemento]=='NGA-West2'):
                    val=int(record_sequence_number_NGA[elemento])
                    f.write("{} {} {} {} {} {} {} {:6.2f} \n".format(i+1,source[elemento],blank,blank,val,event_mag[elemento],acc_distance[elemento],finalScaleFactors[i]))
        f.close()

        if(scale_accelerograms==1):

            #Read accelegrams, save them and apply scaling factor
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
                elemento=recIdx[i]
                if(source[elemento]=='NGA-West2'):
                    val=int(record_sequence_number_NGA[elemento])
                    [desc1,desc2,time1,time2,inp_acc1, inp_acc2,npts1,npts2]=create_NGA_acc(val,path_NGA_folder)
                    desc1='%'+desc1
                    desc2='%'+desc2
                if(source[elemento]=='ESM'):
                    folder_output=output_ESM+event_id[elemento]+station_code[elemento]
                    if (os.path.isdir(folder_output)==False):
                        zip_output='output_'+str(i)+'.zip'
                        command='curl -X POST -F "message=@token.txt" "http://tex.mi.ingv.it/esmws/eventdata/1/query?eventid='+event_id[elemento]+'&data-type=ACC&station='+station_code[elemento]+'&format=ascii" -o '+zip_output
                        os.system(command)
                        command='unzip -o '+zip_output+' -d '+folder_output
                        os.system(command)
                        command='rm '+zip_output
                        os.system(command)
                    [time1,time2,inp_acc1,inp_acc2,npts1,npts2,comp1,comp2]=create_ESM_acc(folder_output)
                    desc1='%'+event_id[elemento]+' '+station_code[elemento]+' '+comp1
                    desc2='%'+event_id[elemento]+' '+station_code[elemento]+' '+comp2

                # Create the filenames
                file_time_scaled_acc_out_1=name+'/GMR_time_scaled_acc_'+str(i+1)+'_1.txt'
                file_time_scaled_acc_out_2=name+'/GMR_time_scaled_acc_'+str(i+1)+'_2.txt'

                with open(file_time_scaled_acc_out_1, "w",newline='') as f1:
                    f1.write("{} \n".format(desc1))
                    for j in np.arange(npts1):
                        f1.write("{:10.3f} {:15.10f}\n".format(time1[j],inp_acc1[j]*finalScaleFactors[i]))
                f1.close()
                with open(file_time_scaled_acc_out_2, "w",newline='') as f2:
                    f2.write("{} \n".format(desc2))
                    for j in np.arange(npts2):
                        f2.write("{:10.3f} {:15.10f}\n".format(time2[j],inp_acc2[j]*finalScaleFactors[i]))
                f2.close()
