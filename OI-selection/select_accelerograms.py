#%% Import libraries
# Standard built-in libraries
import os,sys
import numpy as np
import pandas as pd
from scipy.stats import skew

# Libraries from other files in lib
sys.path.append("lib")
from openquake.hazardlib import gsim, imt, const
from im_correlation import akkar_correlation
from im_correlation import baker_jayaram_correlation
from screen_database import screen_database
from compute_avgSA import compute_avgSA
from compute_avgSA import compute_rho_avgSA
from simulate_spectra import simulate_spectra
from scale_acc import scale_acc
from plot_final_selection import plot_final_selection

#%% General notes
print('Usage: python select_accelerograms.py job_selection.ini')
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
site_code= np.array(site_code,dtype=int)
rlz_code = [ x.strip() for x in input['rlz_code'].strip('{}').split(',') ]
rlz_code= np.array(rlz_code,dtype=int)
if(len(rlz_code)!=len(site_code)):
    sys.exit('Error: rlz_code must be an array of the same length of site_code')
path_results_classical=input['path_results_classical']
path_results_disagg=input['path_results_disagg']
num_disagg=int(input['num_disagg'])
num_classical=int(input['num_classical'])
probability_of_exceedance_num = [ x.strip() for x in input['probability_of_exceedance_num'].strip('{}').split(',') ]
probability_of_exceedance_num= np.array(probability_of_exceedance_num,dtype=int)
probability_of_exceedance = [ x.strip() for x in input['probability_of_exceedance'].strip('{}').split(',') ]
if(len(probability_of_exceedance_num)!=len(probability_of_exceedance_num)):
    sys.exit('Error: probability_of_exceedance_num must be of the same size of probability_of_exceedance')
investigation_time=float(input['investigation_time'])


# Conditional spectrum parameters
period_range = [ x.strip() for x in input['period_range'].strip('{}').split(',') ]
minT=float(period_range[0])
maxT=float(period_range[1])

Tstar=np.zeros(len(intensity_measures))
im_type = []
im_type_lbl = []
for i in np.arange(len(intensity_measures)):
	if(intensity_measures[i]=='AvgSA'):
		im_type.append('AvgSA')
		im_type_lbl.append(r'AvgSa')
	elif(intensity_measures[i][0:2]=='SA'):
		im_type.append('SA')
		im_type_lbl.append(r'Sa(T)')
		Tstar[i]=intensity_measures[i].strip('(,),SA')
		Tstar[i]=float(Tstar[i])
	elif(intensity_measures[i][0:3]=='PGA'):
		im_type.append('PGA')
		im_type_lbl.append(r'PGA')
		print('### Should we add Tstar = 0 for PGA? If we do, the akkar_correlation function gives an error because it does not support T=0')
#		Tstar[i]=0.0
	else:
		sys.exit('Error: this intensity measure type '+intensity_measures[i]+' is not supported yet')

corr_type=input['corr_type']  #baker_jayaram or akkar
if(maxT>4.0 and corr_type=='akkar'):
	sys.exit('Error: akkar correlation model is defined only for T<4s')
GMPE_input=input['GMPE'] #for now only available Akkar_Bommer_2010 Chiou_Youngs_2014 Boore_et_al_2014, BooreAtkinson2008 maybe and array of GMPE according with sites?
avg_periods = [ x.strip() for x in input['avg_periods'].strip('{}').split(',') ]
avg_periods= np.array(avg_periods,dtype=float)
#avg_periods= input['avg_periods'] # Gives an error with the above but works with this
rake=float(input['rake'])
Vs30_input=[ x.strip() for x in input['Vs30'].strip('{}').split(',') ] 
if(len(Vs30_input)!=len(site_code)):
	sys.exit('Error: Vs30 must be an array of the same length of site_code')

vs30Type=[ x.strip() for x in input['vs30Type'].strip('{}').split(',') ] 
if(len(vs30Type)!=len(site_code)):
	sys.exit('Error: Vs30 must be an array of the same length of site_code')

try:
	hypo_depth=float(input['hypo_depth'])
	hypo_defined=1
except KeyError:
	print('Warning: if used, the hypocentral depth will be defined inside the code')
	hypo_defined=0

try:
	dip_input=float(input['dip'])
	dip_defined=1
except KeyError:
	print('Warning: if used, the dip angle will be defined inside the code')
	dip_defined=0

try:
	azimuth=float(input['azimuth'])
except KeyError:
	try:
		FHW=int(input['hanging_wall_flag'])
		if(FHW==1):
			azimuth=50
		elif(FHW==-1):
			azimuth=-50
		else:
			sys.exit('Error: The hanging_wall_flag must be =1 or =-1')
	except KeyError:
		sys.exit('Error: The azimuth or the hanging_wall_flag must be defined')

z1pt0_defined=-1
z2pt5_defined=-1
z2pt5=[]
z1pt0=[]

if(GMPE_input=='CampbellBozorgnia2008' or GMPE_input=='CampbellBozorgnia2014'):
	try:
		z2pt5_input=[ x.strip() for x in input['z2pt5'].strip('{}').split(',') ]
		z2pt5_defined=1
		if(len(z2pt5_input)!=len(site_code)):
			sys.exit('Error: z2pt5 must be an array of the same length of site_code')
	except KeyError:
		print('Warning: z2pt5 will be defined inside the code')
		z2pt5_defined=0
if(GMPE_input=='AbrahamsonEtAl2014' or GMPE_input=='ChiouYoungs2014'):
	try:
		z1pt0_input=[ x.strip() for x in input['z1pt0'].strip('{}').split(',') ]
		z1pt0_defined=1
		if(len(z1pt0_input)!=len(site_code)):
			sys.exit('Error: z2pt5 must be an array of the same length of site_code')
	except KeyError:
		print('Warning: z1pt0 will be defined inside the code')
		z1pt0_defined=0

try:
	upper_sd=float(input['upper_sd'])
except KeyError:
	upper_sd=0

try:
	lower_sd=float(input['lower_sd'])
except KeyError:
	lower_sd=500

# Database parameters for screening recordings
database_path=input['database_path']
allowed_database = [ x.strip() for x in input['allowed_database'].strip('{}').split(',') ]

allowedRecs_Vs30=''
allowedEC8code=''
try:
	allowedRecs_Vs30 = [ x.strip() for x in input['allowedRecs_Vs30'].strip('[]').split(',') ]# upper and lower bound of allowable Vs30 values
	allowedRecs_Vs30= np.array(allowedRecs_Vs30,dtype=float)
	allowedRecs_Vs30_defined=1
except KeyError:
	allowedRecs_Vs30_defined=0

try:
	allowedEC8code = [ x.strip() for x in input['allowedEC8code'].strip('{}').split(',') ]
	allowedEC8code_defined=1
except KeyError:
	allowedEC8code_defined=0

try:
	maxsf_input=float(input['maxsf'])
except ValueError:
	maxsf_input=[ x.strip() for x in input['maxsf'].strip('{}').split(',') ] #The maximum allowable scale factor. They must be specified as a function of probability_of_exceedance
	maxsf_input= np.array(maxsf_input,dtype=float)
	if(len(probability_of_exceedance_num)!=len(maxsf_input)):
		sys.exit('Error: maxsf must be of the same size of probability_of_exceedance')
#Maybe we can give the possibility to give different SF for different IM, sites, ecc.. and also for radius_dist and mag
try:
	radius_dist_input=float(input['radius_dist'])
except ValueError:
	radius_dist_input=[ x.strip() for x in input['radius_dist'].strip('{}').split(',') ] #km
	radius_dist_input= np.array(radius_dist_input,dtype=float)
	if(len(probability_of_exceedance_num)!=len(radius_dist_input)):
		sys.exit('Error: radius_dist must be of the same size of probability_of_exceedance')
try:
	radius_mag_input=float(input['radius_mag'])
except ValueError:
	radius_mag_input=[ x.strip() for x in input['radius_mag'].strip('{}').split(',') ] 
	radius_mag_input= np.array(radius_mag_input,dtype=float)
	if(len(probability_of_exceedance_num)!=len(radius_mag_input)):
		sys.exit('Error: radius_mag must be of the same size of probability_of_exceedance')

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
path_ESM_folder=input['path_ESM_folder']
# If not, found in the folder, ESM recording are authomatically downloaded from internet, need to generate the file token.txt
#At first you need to register at: https://tex.mi.ingv.it/
#curl -X POST -F 'message={"user_email": "email","user_password": "password"}' "https://tex.mi.ingv.it/esmws/generate-signed-message/1/query" > token.txt

# Output folder
output_folder=input['output_folder']

#%% Start the routine
print('Inputs loaded, starting selection....')
ind = 1

# For each site investigated
for ii in site_code:

	# Get the current site and realisation indices
	site = site_code[ii]
	rlz = rlz_code[ii]

	# For each hazard of poe level investigated
	count_poe=-1
	for poe in probability_of_exceedance_num:
		count_poe=count_poe+1
		if hasattr(maxsf_input, '__len__'):
			maxsf=maxsf_input[count_poe]
		else:
			maxsf=maxsf_input
		if hasattr(radius_dist_input, '__len__'):
			radius_dist=radius_dist_input[count_poe]
		else:
			radius_dist=radius_dist_input
		if hasattr(radius_mag_input, '__len__'):
			radius_mag=radius_mag_input[count_poe]
		else:
			radius_mag=radius_mag_input

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
			df=pd.read_csv(''.join([path_results_disagg,'/',disagg_results]),skiprows=1)
			df['rate'] = -np.log(1-df['poe'])/investigation_time
			df['rate_norm'] = df['rate']/ df['rate'].sum()
			mode=df.sort_values(by='rate_norm',ascending=False)[0:1]
			meanMag=np.sum(df['mag']*df['rate_norm'])
			meanDist=np.sum(df['dist']*df['rate_norm'])
			allowedRecs_D=[meanDist-radius_dist,meanDist+radius_dist]
			allowedRecs_Mag=[meanMag-radius_mag,meanMag+radius_mag]

			#Retrieve conditioning value (this is the Sa(Tstar), for example)
			df=pd.read_csv(''.join([path_results_classical,'/',file_with_OQ_acc_value]),skiprows=1)
			output_oq=df[selected_column]
			im_star=output_oq[0]

#%% -----------------------------------------------------------------------------
			# Initialise GSIMs that are needed (i.e. the GMPEs)

			#GSIM = gsim.get_available_gsims()

			for name_gmpe, gmpes in gsim.get_available_gsims().items():
				if name_gmpe==GMPE_input:
					bgmpe=gmpes

			# Set up the site, rupture and distance contexts
			sctx = gsim.base.SitesContext()
			rctx = gsim.base.RuptureContext()
			dctx = gsim.base.DistancesContext()

# -----------------------------------------------------------------------------
			# Initialise contexts
			rjb=np.arange(meanDist,meanDist+1,1.)
			mag=meanMag
			if(hypo_defined==1):
				Z_hyp=hypo_depth
			else:
				if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
					Z_hyp=5.63+0.68*mag
				else:
					Z_hyp=11.24-0.2*mag

			if(dip_defined==1):
				dip=dip_input
			else:
				if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
					dip=90
				elif RakeAverage > 0:
					dip=40
				else:
					dip=50

			if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
				 # strike slip
				width= 10.0 ** (-0.76 + 0.27 *mag)
			elif RakeAverage > 0:
				# thrust/reverse
				width= 10.0 ** (-1.61 + 0.41 *mag)
			else:
				# normal
				width= 10.0 ** (-1.14 + 0.35 *mag)

			source_vertical_width=width*np.sin(np.radians(dip))
			ztor=max(Z_hyp-0.6*source_vertical_width,upper_sd)
			if((ztor+source_vertical_width)>lower_sd):
				source_vertical_width=lower_sd-ztor
				width=source_vertical_width/np.sin(np.radians(dip))

			if(rjb==0):
				rx=0.5*width*np.cos(np.radians(dip))
			else:
				if(dip==90):
					rx=rjb*np.sin(np.radians(azimuth))
				else:
					if (azimuth>=0 and azimuth<90) or (azimuth>90 and azimuth<=180):
						if(rjb*np.abs(np.tan(np.radians(azimuth)))<=width*np.cos(np.radians(dip))):
							rx=rjb*np.abs(np.tan(np.radians(azimuth)))
						else:
							rx=rjb*np.tan(np.radians(azimuth))*np.cos(np.radians(azimuth)-np.arcsin(width*np.cos(np.radians(dip))*np.cos(np.radians(azimuth))/rjb))
					elif (azimuth==90): #we assume that Rjb>0
						rx=rjb+width*np.cos(np.radians(dip))
					else:
						rx=rjb*np.sin(np.radians(azimuth))

			if(azimuth==90 or azimuth==-90):
				ry=0
			elif(azimuth==0 or azimuth==180 or azimuth==-180):
				ry=rjb
			else:
				ry=np.abs(rx*1./np.tan(np.radians(azimuth)))

			if(dip==90):
				rrup=np.sqrt(np.square(rjb)+np.square(ztor))
			else:
				if(rx<ztor*np.tan(np.radians(dip))):
					rrup1=np.sqrt(np.square(rx)+np.square(ztor))
				if(rx>=ztor*np.tan(np.radians(dip)) and rx<=ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
					rrup1=rx*np.sin(np.radians(dip))+ztor*np.cos(np.radians(dip))
				if(rx>ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
					rrup1=np.sqrt(np.square(rx-width*np.cos(np.radians(dip)))+np.square(ztor+width*np.sin(np.radians(dip))))
				rrup=np.sqrt(np.square(rrup1)+np.square(ry))

			Vs30=float(Vs30_input[ii])
			if(vs30Type[ii]=='inferred'):
				vs30measured=1
			if(vs30Type[ii]=='measured'):
				vs30measured=0

			if(z1pt0_defined==0 and GMPE_input=='AbrahamsonEtAl2014') or (z2pt5_defined==0):
				if(Vs30<180):
					z1pt0=np.exp(6.745)
				elif(Vs30>=180 and Vs30<=500):
					z1pt0=np.exp(6.745-1.35*np.log(Vs30/180))
				else:
					z1pt0=np.exp(5.394-4.48*np.log(Vs30/500))

			if(z1pt0_defined==0 and GMPE_input=='ChiouYoungs2014'):
				z1pt0=np.exp(28.5-3.82/8*np.log(Vs30**8+378.7**8))

			if(z1pt0_defined==1):
				z1pt0=float(z1pt0_input[ii])

			if(z2pt5_defined==0):
				z2pt5=519+3.595*z1pt0

			if(z2pt5_defined==1):
				z2pt5=float(z2pt5_input[ii])

			setattr(rctx, 'width', width)
			setattr(rctx, 'ztor', ztor)
			setattr(rctx, 'dip', dip)
			setattr(dctx, 'rx', rx)
			setattr(dctx, 'rrup', rrup)
			setattr(dctx, 'ry0', ry)
			z1pt0=z1pt0+np.zeros(rjb.shape)
			setattr(sctx, 'z1pt0', z1pt0)
			z2pt5=z2pt5+np.zeros(rjb.shape)
			setattr(sctx, 'z2pt5', z2pt5)
			setattr(sctx, 'vs30measured',vs30measured)
			setattr(rctx, 'mag', mag)
			setattr(rctx, 'hypo_depth', Z_hyp)
			setattr(rctx, 'rake', rake)
			setattr(dctx, 'rjb', rjb)
			Vs30=Vs30+np.zeros(rjb.shape)
			setattr(sctx, 'vs30', Vs30)

#%% -----------------------------------------------------------------------------
			# Screen the database of available ground motions
			[SaKnown,indPer,TgtPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance]=screen_database(database_path,allowed_database,allowedRecs_Vs30,allowedRecs_Mag,allowedRecs_D,allowedEC8code,minT,maxT,nGM,allowed_depth)

#%% -----------------------------------------------------------------------------
			# Get the GMPE ouput
			if(im_type[im]=='AvgSA'):
				[median_im_cond, sigma_im_cond] = compute_avgSA(avg_periods,sctx, rctx, dctx, bgmpe, corr_type)
				mu_im_cond = np.log(median_im_cond)
			else:
				if(im_type[im]=='PGA'):
					P = imt.PGA()
				if(im_type[im]=='SA'):
					P = imt.SA(period=Tstar[im])
				S = [const.StdDev.TOTAL]
				mu_im_cond, sigma_im_cond = bgmpe.get_mean_and_stddevs(sctx,rctx,dctx,P,S)
				sigma_im_cond = sigma_im_cond[0]

			# Compute how many standard deviations the PSHA value is above the GMPE value
			epsilon = (np.log(im_star) - mu_im_cond)/sigma_im_cond

#                if(intensity_measures[im]=='AvgSA'):
#                    [mean_SaTcond,stddvs_SaTcond]=compute_avgSA(avg_periods,sctx, rctx, dctx, bgmpe, corr_type)
#                    mean_SaTcond=np.log(mean_SaTcond)
#                else:
#                    if(intensity_measures[im]=='PGA'):
#                        P = imt.PGA()
#                    else:
#                        P = imt.SA(period=Tstar[im])
#                    S=[const.StdDev.TOTAL]
#                    mean_SaTcond,stddvs_SaTcond=bgmpe.get_mean_and_stddevs(sctx,rctx,dctx,P,S)
#                    stddvs_SaTcond=stddvs_SaTcond[0]
#                epsilon=(np.log(im_star)-mean_SaTcond)/stddvs_SaTcond

#%% -----------------------------------------------------------------------------
			# Create the conditional spectrum
#			T_CS = np.arange(0.01,10.0,0.01)
			T_CS = TgtPer # Use the same periods as the available spectra to construct the conditional spectrum
			mu_im = np.zeros(len(T_CS))
			sigma_im = np.zeros(len(T_CS))
			rho_T_Tstar = np.zeros(len(T_CS))
			mu_im_im_cond = np.zeros(len(T_CS))

			for i in range(len(T_CS)):
				# Get the GMPE ouput for a rupture scenario
				if im_type[im] == 'PGA':
					P = imt.PGA()
				elif im_type[im] == 'AvgSA' or im_type[im] == 'SA':
					P = imt.SA(period=T_CS[i])
				S = [const.StdDev.TOTAL]
				mu0,sigma0 = bgmpe.get_mean_and_stddevs(sctx,rctx,dctx,P,S)
				mu_im[i] = mu0[0]
				sigma_im[i] = sigma0[0][0]

				# Compute the correlations between each T and Tstar
				if im_type[im] == 'AvgSA':
					rho = compute_rho_avgSA(T_CS[i],avg_periods,sctx,rctx,dctx,sigma_im_cond, bgmpe, corr_type)
					rho = rho[0]
				elif im_type[im] == 'SA' or im_type[im] == 'PGA':
					if corr_type == 'baker_jayaram':
						rho = baker_jayaram_correlation(T_CS[i],Tstar)
					elif corr_type == 'akkar':
						rho = akkar_correlation(T_CS[i],Tstar)
				rho_T_Tstar[i] = rho

				# Get the value of the CMS
				mu_im_im_cond[i] = mu_im[i] + rho_T_Tstar[i]*epsilon[0]*sigma_im[i]

#                TgtMean=[]
#                rho=[]
#                mean=[]
#                sigma=[]
#                for per in TgtPer:
#                    if(per==0):
#                        P = imt.PGA()
#                    else:
#                        P=imt.SA(period=per)
#                    S=[const.StdDev.TOTAL]
#                    bMean_SA, bStDev_SA = bgmpe.get_mean_and_stddevs(sctx, rctx, dctx, P, S)
#                    mean.append(bMean_SA)
#                    sigma.append(bStDev_SA[0])
#                    if(intensity_measures[im]=='AvgSA'):
#                        rho_per=compute_rho_avgSA(per,avg_periods,sctx,rctx,dctx,stddvs_SaTcond, bgmpe, corr_type)
#                        rho.append(rho_per[0])
#                    else:
#                        if(corr_type=='baker_jayaram'):
#                            rho_per = baker_jayaram_correlation(per,Tstar[im])
#                        if(corr_type=='akkar'):
#                            rho_per = akkar_correlation(per,Tstar[im])
#                        rho.append(rho_per)
#                    spectrum=bMean_SA+rho_per*bStDev_SA[0]*epsilon
#                    # (Log) Response Spectrum Mean: TgtMean
#                    TgtMean.append(spectrum[0])
#
#                TgtMean=np.array(TgtMean)

			# Compute the Covariance
			Cov = np.zeros((len(T_CS),len(T_CS)))

			for i in range(len(T_CS)):
				for j in range(len(T_CS)):
					var1 = sigma_im[i]**2
					var2 = sigma_im[j]**2
					varTstar = sigma_im_cond**2

					if corr_type == 'baker_jayaram':
						sigma_Corr = baker_jayaram_correlation(T_CS[i],T_CS[j])*np.sqrt(var1*var2)
					elif corr_type == 'akkar':
						sigma_Corr = akkar_correlation(T_CS[i],T_CS[j])*np.sqrt(var1*var2)

					Sigma11 = np.matrix([[var1, sigma_Corr], [sigma_Corr, var2]])
					Sigma22 = np.array([varTstar])
					Sigma12 = np.array([rho_T_Tstar[i]*np.sqrt(var1*varTstar), rho_T_Tstar[j]*np.sqrt(varTstar*var2)])
					Sigma_cond = Sigma11 - Sigma12*1./(Sigma22)*Sigma12.T
					Cov[i,j]=Sigma_cond[0,1]

			Cov[np.absolute(Cov) < 1e-10] = 1e-10 # find covariance values of zero and set them to a small number so that random number generation can be performed
			Sigma_im_im_cond = np.sqrt(np.diagonal(Cov))

#                ### Compute covariances and correlations at all periods
#                TgtCovs = np.zeros((len(TgtPer),len(TgtPer)))
#                for i in np.arange(len(TgtPer)):
#                    for j in np.arange(len(TgtPer)):
#                        Ti = TgtPer[i]
#                        Tj = TgtPer[j]
#                        varT = stddvs_SaTcond**2
#                        varT=varT[0]
#                        var1 = sigma[i]**2
#                        var1=var1[0]
#                        var2 = sigma[j]**2
#                        var2=var2[0]
#
#                        if(corr_type=='baker_jayaram'):
#                            sigmaCorr=baker_jayaram_correlation(Ti, Tj)*np.sqrt(var1*var2)
#                        if(corr_type=='akkar'):
#                            sigmaCorr=akkar_correlation(Ti, Tj)*np.sqrt(var1*var2)
#                        sigma11 = np.matrix([[var1, sigmaCorr], [sigmaCorr, var2]])
#                        sigma22 = np.array(varT)
#                        sigma12 =np.array([[rho[i]*np.sqrt(var1*varT)],[rho[j]*np.sqrt(var2*varT)]])
#                        sigmaCond = sigma11 - sigma12*1./(sigma22)*sigma12.T
#                        TgtCovs[i,j] = sigmaCond.item((0,1))

#%% -----------------------------------------------------------------------------
			# Simulate and find suitable spectra to match the target CS defined above
			meanReq=mu_im_im_cond
			covReq=Cov
			stdevs=np.sqrt(np.diagonal(covReq))

			# Simulate response spectra using the target CS
			simulated_spectra=simulate_spectra(nTrials,meanReq,covReq,stdevs,nGM,weights)
			sampleBig =np.log(SaKnown[:,indPer])

			id_sel=[]
			if im_type[im] == 'AvgSA':
				id_sel_bool=np.isin(TgtPer,avg_periods)
				for i in np.arange(len(TgtPer)):
					if(id_sel_bool[i] == True):
						id_sel.append(i)
				id_sel=np.array(id_sel)
			elif im_type[im] == 'SA' or im_type[im] == 'PGA':
				id_sel=np.where(TgtPer==Tstar[im])
			lnSa1=np.mean(meanReq[id_sel])

			recID = np.zeros(nGM,dtype=int)
			sampleSmall = []
			IMScaleFac = np.ones(nGM)

			# Find database spectra most similar to each simulated spectrum
			for i in np.arange(nGM): # for each simulated spectrum
				err=np.zeros(nBig)*1000000 #initialize error matrix
				scaleFac = np.ones(nBig) # initialize scale factors to 1
				# compute scale factors and errors for each candidate ground motion
				for j in np.arange(nBig):
					rec_value=np.exp(sum(sampleBig[j,id_sel])/len(id_sel))
					if type(rec_value) == 'list':
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
#%% -----------------------------------------------------------------------------
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

#%% -----------------------------------------------------------------------------
			# Create the output information

			# Collect information of the final record set
			finalRecords = recID
			recIdx = [allowedIndex[i] for i in finalRecords]
			finalScaleFactors = IMScaleFac
			meanrecorded=np.mean(np.exp(sampleSmall),axis=0)
			meanrecorded_p2sigma = np.percentile(np.exp(sampleSmall),50+34.1+13.6,axis=0)
			meanrecorded_n2sigma = np.percentile(np.exp(sampleSmall),50-34.1-13.6,axis=0)
			meanrecorded_eps =  (np.log(meanrecorded_p2sigma)-np.log(meanrecorded_n2sigma))/(2*1.96)

			# Create the outputs folder
			folder=output_folder+'/'+name
			if not os.path.exists(folder):
				os.makedirs(folder)

			# Plot the figure
			plot_final_selection(name,im_type_lbl[im],nGM,TgtPer,T_CS,sampleSmall,meanReq,stdevs,meanrecorded,meanrecorded_p2sigma,meanrecorded_n2sigma,meanrecorded_eps,output_folder)

			# Output a summary of the results to a text file
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
				scale_acc(nGM,recIdx,record_sequence_number_NGA,path_NGA_folder,path_ESM_folder,source,event_id,station_code,name,output_folder,finalScaleFactors)
>>>>>>> 40c4f7a2b40068964caee3b4841dff76c5af59f3
