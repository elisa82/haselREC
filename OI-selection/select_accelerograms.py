#%% Import libraries
# Standard built-in libraries
import os,sys
import numpy as np
import pandas as pd
from scipy.stats import skew

# Libraries from other files in bin
sys.path.append("bin")
from openquake.hazardlib import gsim, imt, const
from im_correlation import akkar_correlation
from im_correlation import baker_jayaram_correlation
from screen_database import screen_database
from create_acc import create_ESM_acc
from create_acc import create_NGA_acc
from compute_avgSA import compute_avgSA
from compute_avgSA import compute_rho_avgSA
from simulate_spectra import simulate_spectra
from scale_acc import scale_acc
from plot_final_selection import plot_final_selection

#%% General notes
print('Usage: python select_accelerograms.py $filename')
print("### The correlation model needs to be handled a bit better")
print("### The akkar one could be interpolated in its python script definition (I think)")


print("### The outputs and printing of results need to be more detailed here")

#%% Initial setup
#fileini = sys.argv[1]
fileini = 'job_selection.ini'
print("### Define directly for now")

input={}
with open(fileini) as fp:
   line = fp.readline()
   while line:
       if line.strip().find('=')>=0:
            key, value =line.strip().split('=', 1)
            input[key.strip()] = value.strip()
       line = fp.readline()

#%% Extract input parameters
# Calculation mode
download_and_scale_acc=int(input['download_and_scale_acc']) #1=yes, 0=no


# Hazard parameters
intensity_measures = [ x.strip() for x in input['intensity_measures'].strip('{}').split(',') ]
site_code = [ x.strip() for x in input['site_code'].strip('{}').split(',') ]
rlz_code = [ x.strip() for x in input['rlz_code'].strip('{}').split(',') ]
path_hazard_results=input['path_hazard_results']
num_disagg=int(input['num_disagg'])
num_classical=int(input['num_classical'])
probability_of_exceedance_num = [ x.strip() for x in input['probability_of_exceedance_num'].strip('{}').split(',') ]
probability_of_exceedance = [ x.strip() for x in input['probability_of_exceedance'].strip('{}').split(',') ]
investigation_time=float(input['investigation_time'])


# Conditional spectrum parameters
period_range = [ x.strip() for x in input['period_range'].strip('{}').split(',') ]
minT=float(period_range[0])
maxT=float(period_range[1])


Tstar=np.zeros(len(intensity_measures))
for i in np.arange(len(intensity_measures)):
    if(intensity_measures[i][0:2]=='SA'):
        Tstar[i]=intensity_measures[i].strip('(,),SA')
        Tstar[i]=float(Tstar[i])
#    if(intensity_measures[i][0:3]=='PGA'):
#	    Tstar[i]=0.0
print('### Should we add Tstar = 0 for PGA? If we do, the akkar_correlation function gives an error because it does not support T=0')

corr_type=input['corr_type']  #baker_jayaram or akkar
if(maxT>4.0 and corr_type=='akkar'):
    sys.exit('Error: akkar correlation model is defined only for T<4s')
GMPE=input['GMPE'] #for now only available Akkar_Bommer_2010 Chiou_Youngs_2014 Boore_et_al_2014, BooreAtkinson2008 maybe and array of GMPE according with sites?
avg_periods = [ x.strip() for x in input['avg_periods'].strip('{}').split(',') ]
avg_periods= np.array(avg_periods,dtype=float)
rake=float(input['rake'])
hypo_depth=float(input['hypo_depth'])
AR=float(input['AR'])
Vs30=float(input['Vs30']) #maybe an array of Vs30 according with sites?


# Database parameters for screening recordings
database_path=input['database_path']
allowed_database = [ x.strip() for x in input['allowed_database'].strip('{}').split(',') ]
allowedRecs_Vs30 = [ x.strip() for x in input['allowedRecs_Vs30'].strip('[]').split(',') ]# upper and lower bound of allowable Vs30 values
allowedRecs_Vs30= np.array(allowedRecs_Vs30,dtype=float)
allowedEC8code = [ x.strip() for x in input['allowedEC8code'].strip('{}').split(',') ]
maxsf=float(input['maxsf']) #The maximum allowable scale factor
#Maybe we can give the possibility to give different SF for different IM, sites, ecc.. and also for radius_dist and mag
radius_dist= float(input['radius_dist']) # km
radius_mag = float(input['radius_mag'])
allowed_depth=[ x.strip() for x in input['allowed_depth'].strip('[]').split(',') ] # upper and lower bound of allowable Vs30 values
allowed_depth= np.array(allowed_depth,dtype=float)


# Selection parameters
nGM=int(input['nGM']) #number of records to select ==> number of records to select since the code search the database spectra most similar to each simulated spectrum
random_seed=int(input['random_seed']) 
nTrials=int(input['nTrials']) #number of iterations of the initial spectral simulation step to perform
weights= [ x.strip() for x in input['weights'].strip('{}').split(',') ]# [Weights for error in mean, standard deviation and skewness] Used to find the simulated spectra that best match the target from the statistically simulated response spectra
weights= np.array(weights,dtype=float)
#otimization parameters. Execution of incremental changes to the initially selected ground motion set to further optimise its fit to the target spectrum distribution.
nLoop=int(input['nLoop']) #Number of loops of optimization to perform
penalty=float(input['penalty']) #>0 to penalize selected spectra moire than 3 sigma from the target at any period, =0 otherwise.


# Accelerogram folders
path_NGA_folder=input['path_NGA_folder'] #NGA recordings have to be stored
path_ESM_folder=input['path_NGA_folder']
# If not, found in the folder, ESM recording are authomatically downloaded from internet, need to generate the file token.txt
#At first you need to register at: http://tex.mi.ingv.it/
#curl -X POST -F 'message={"user_email": "elisa.zuccolo@eucentre.it","user_password": "password"}' "http://tex.mi.ingv.it/esmws/generate-signed-message/1/query" > token.txt

# Output folder
output_folder=input['output_folder']

#%% Start the routine
print('Inputs loaded, starting selection....')
ind = 1

# For each site investigated
for ii in np.arange(len(site_code)):
    
    # Get the current site and realisation indices
    site = site_code[ii]
    rlz = rlz_code[ii]
    
    # For each hazard of poe level investigated
    for poe in np.arange(len(probability_of_exceedance_num)):
        
        # For each intensity measure investigated
        for im in np.arange(len(intensity_measures)):
                
                # Get the name of the disaggregation file to look in
                disagg_results='rlz-'+str(rlz)+'-'+intensity_measures[im]+'-sid-'+str(site)+'-poe-0_Mag_Dist_'+str(num_disagg)+'.csv'
                name=intensity_measures[im]+'-site_'+str(site)+'-poe-'+str(poe)
                selected_column=intensity_measures[im]+'-'+str(probability_of_exceedance[poe])
                file_with_OQ_acc_value='hazard_map-mean_'+str(num_classical)+'.csv'
                
                # Print some on screen feedback
                print('Processing '+name+' Case: '+str(ind)+'/'+str(len(site_code)*len(probability_of_exceedance_num)*len(intensity_measures)))
                ind += 1

                #Retrieve disaggregation results
                meanLst = [],[]
                df=pd.read_csv(''.join([path_hazard_results,'/',disagg_results]),skiprows=1)
                df['rate'] = -np.log(1-df['poe'])/investigation_time
                df['rate_norm'] = df['rate']/ df['rate'].sum()
                mode=df.sort_values(by='rate_norm',ascending=False)[0:1]
                meanMag=np.sum(df['mag']*df['rate_norm'])
                meanDist=np.sum(df['dist']*df['rate_norm'])
#               print(meanMag,meanDist)
                allowedRecs_D=[meanDist-radius_dist,meanDist+radius_dist]
                allowedRecs_Mag=[meanMag-radius_mag,meanMag+radius_mag]

                #Retrieve conditioning value
                df=pd.read_csv(''.join([path_hazard_results,'/',file_with_OQ_acc_value]),skiprows=1)
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
                if(GMPE=='BooreAtkinson2008'):
                    bgmpe = gsim.boore_atkinson_2008.BooreAtkinson2008()
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

#            print(rjb,rx,rrup)

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

                [SaKnown,indPer,TgtPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance]=screen_database(database_path,allowed_database,allowedRecs_Vs30,allowedRecs_Mag,allowedRecs_D,allowedEC8code,minT,maxT,nGM,allowed_depth)

                TgtMean=[]
                rho=[]
                mean=[]
                sigma=[]
                if(intensity_measures[im]=='AvgSA'):
                    [mean_SaTcond,stddvs_SaTcond]=compute_avgSA(avg_periods,sctx, rctx, dctx, bgmpe, corr_type)
                    mean_SaTcond=np.log(mean_SaTcond)
                else:
                    if(intensity_measures[im]=='PGA'):
                        P = imt.PGA()
                    else:
                        P = imt.SA(period=Tstar[im])
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
                    if(intensity_measures[im]=='AvgSA'):
                        rho_per=compute_rho_avgSA(per,avg_periods,sctx,rctx,dctx,stddvs_SaTcond, bgmpe, corr_type)
                        rho.append(rho_per[0])
                    else:
                        if(corr_type=='baker_jayaram'):
                            rho_per = baker_jayaram_correlation(per,Tstar[im])
                        if(corr_type=='akkar'):
                            rho_per = akkar_correlation(per,Tstar[im])
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

                simulated_spectra=simulate_spectra(nTrials,meanReq,covReq,stdevs,nGM,weights)

                sampleBig =np.log(SaKnown[:,indPer])

                id_sel=[]
                if(intensity_measures[im]=='AvgSA'):
                    id_sel_bool=np.isin(TgtPer,avg_periods)
                    for i in np.arange(len(TgtPer)):
                        if(id_sel_bool[i] == True):
                            id_sel.append(i)
                    id_sel=np.array(id_sel)
                else:
                    id_sel=np.where(TgtPer==Tstar[im])
                lnSa1=np.mean(meanReq[id_sel])

                recID = np.zeros(nGM,dtype=int)
                sampleSmall = []
                IMScaleFac = np.ones(nGM)
# Find database spectra most similar to each simulated spectrum
                for i in np.arange(nGM): # for each simulated spectrum
                    err=np.zeros(nBig)*1000000 #initialize error matrix
                    scaleFac = np.ones(nBig) # initialize scale factors to 1
                    #compute scale factors and errors for each candidate ground motion
                    for j in np.arange(nBig):
                        rec_value=np.exp(sum(sampleBig[j,id_sel])/len(id_sel))
                        print(rec_value,id_sel)
                        rec_value=rec_value[0]
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
                            rec_value=np.exp(sum(sampleBig[j,id_sel])/len(id_sel))
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

                folder=output_folder+'/'+name 
                if not os.path.exists(folder):
                    os.makedirs(folder)

                # Plot the figure
                plot_final_selection(name,nGM,TgtPer,sampleSmall,meanReq,stdevs,meanrecorded,output_folder)

                # Output results to a text file
                blank='-'
                name_summary=output_folder+'/'+name+'/'+name+"_summary_selection.txt"
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

                if(download_and_scale_acc==1):

                    scale_acc(nGM,recIdx,record_sequence_number_NGA,path_NGA_folder,path_ESM_folder,event_id,station_code,name,output_folder)
