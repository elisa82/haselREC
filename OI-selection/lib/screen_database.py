def screen_database(database_path,allowed_database,allowedRecs_Vs30,allowedRecs_Mag,allowedRecs_D,allowedEC8code,minT,maxT,nGM,allowed_depth,allowedRecs_Vs30_defined,allowedEC8code_defined,Vs30):
    # Import libraries
    import numpy as np
    import pandas as pd

    knownPer=np.array([0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,5,6,7,8,9,10])
    
    dbacc=pd.read_csv(database_path,sep=';',engine='python')

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

    if(allowedRecs_Vs30_defined==0):
        if(Vs30>=800.0):
            allowedRecs_Vs30=[800.0,3000.0]
        elif(Vs30>=360. and Vs30<800.):
            allowedRecs_Vs30=[360.0,800.0]
        elif(Vs30>=180. and Vs30<360.):
            allowedRecs_Vs30=[180.0,360.0]
        else:
            allowedRecs_Vs30=[0.0,180.0]

    if(allowedEC8code_defined==0):
        if(Vs30>=800.0):
            allowedEC8code='A'
        elif(Vs30>=360. and Vs30<800.):
            allowedEC8code='B'
        elif(Vs30>=180. and Vs30<360.):
            allowedEC8code='C'
        else:
            allowedEC8code='D'

    # select periods in the range
    indPer = []
    for i in np.arange(len(knownPer)):
        if(knownPer[i]>=minT and knownPer[i]<=maxT):
            indPer.append(i)

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
            #print('Need to test if the screening of database is ok')
            if(source[i] in allowed_database):
                if((source[i]=='ESM' and is_free_field_ESM[i]==0) or (source[i]=='NGA-West2' and is_free_field_NGA[i]=="I")):
                    if((source[i]=='ESM' and event_mw[i]>=allowedRecs_Mag[0] and event_mw[i]<=allowedRecs_Mag[1]) or (source[i]=='NGA-West2' and event_mag[i]>=allowedRecs_Mag[0] and event_mag[i]<=allowedRecs_Mag[1])):
                        if(event_depth[i]>=allowed_depth[0] and event_depth[i]<=allowed_depth[1]):
                            if(acc_distance[i]>=allowedRecs_D[0] and acc_distance[i]<=allowedRecs_D[1]):
                                if(np.isnan(station_vs30[i])):
                                    if not pd.isnull(station_ec8[i]):
                                        if(station_ec8[i][0] in allowedEC8code or allowedEC8code=='All'):
                                            if(source[i]=='ESM'):
                                                allowedIndex.append(i)
                                            if(source[i]=='NGA-West2'):
                                                if(epi_lon[i]<-31 or epi_lon[i]>70):
                                                    allowedIndex.append(i)
                                else:
                                    if(station_vs30[i]>=allowedRecs_Vs30[0] and station_vs30[i]<allowedRecs_Vs30[1]):
                                        if(source[i]=='ESM'):
                                            allowedIndex.append(i)
                                        if(source[i]=='NGA-West2'):
                                            if(epi_lon[i]<-31 or epi_lon[i]>70):
                                                allowedIndex.append(i)

    SA=np.vstack(SA_list)
    SaKnown=SA[allowedIndex]

    # count number of allowed spectra
    nBig = len(allowedIndex)
    print(['Number of allowed ground motions = ', nBig])
    assert (nBig >= nGM), 'Warning: there are not enough allowable ground motions'

    return SaKnown,indPer,recPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance,station_vs30,station_ec8
