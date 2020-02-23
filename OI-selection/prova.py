#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:38:53 2020

@author: gor

Plot the spectra of the ground motions identified (provided spectrum)
Compare them with the spectra computed manually from record
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np

#%% Import the spectra used in selection from pickle file
SaKnown,indPer,TgtPer,nBig,allowedIndex,event_id,station_code,source,record_sequence_number_NGA,source,event_mw,event_mag,acc_distance = pickle.load( open( 'provaSA(0.5)-site_1-poe-0.pkl' , 'rb'))
ngms = 2 # Just look a few ground motion

#%% Plot the spectra used in selection
clr = ['b', 'g', 'r']
plt.figure()
for i in range(ngms):
	plt.plot(TgtPer,SaKnown[i][indPer],color=clr[i],linewidth=1,label='RotD50 (OpenInsel)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Period [s]')
plt.ylabel('Spectral Acceleration [g]')

#%% Define function to read the AT2 files directly
def read_NGA(filepath):
	with open(filepath,'r') as f:
		content = f.readlines()
	counter = 0
	desc, row4Val, acc_data = "","",[]
	for x in content:
		if counter == 1:
			desc = x
		elif counter == 3:
			row4Val = x
			if row4Val[0][0] == 'N':
				val = row4Val.split()
				npts = float(val[(val.index('NPTS='))+1].rstrip(','))
				dt = float(val[(val.index('DT='))+1])
			else:
				val = row4Val.split()
				npts = float(val[0])
				dt = float(val[1])
		elif counter > 3:
			data = str(x).split()
			for value in data:
				a = float(value)
				acc_data.append(a)
			inp_acc = np.asarray(acc_data)
			time = []
			for i in range (0,len(acc_data)):
				t = i * dt
				time.append(t)
		counter = counter + 1
	return desc, npts, dt, time, inp_acc


#%% Load the actual time series
a1, a2 = [], []
t1, t2 = [], []
dts = []
for i in range(ngms):
	file1 = '../GM-Records-Database/NGA-West2-Reduced/RSN'+str(record_sequence_number_NGA[allowedIndex[i]])+'_1.AT2'
	file2 = '../GM-Records-Database/NGA-West2-Reduced/RSN'+str(record_sequence_number_NGA[allowedIndex[i]])+'_2.AT2'
	_, npts1, dt, time1, acc1 = read_NGA(file1)
	_, npts2, _, time2, acc2 = read_NGA(file2)
	dts.append(dt)
	npts = int(min(npts1,npts2))
	a1.append(acc1[0:npts])
	a2.append(acc2[0:npts])
	t1.append(time1[0:npts])
	t2.append(time2[0:npts])

#%% Compute the response spectra of the two orthogonal components
import eqsig
num_T = 50
Trange = np.logspace(start=-2, stop=1, num=num_T)

Sa1, Sa2 = [], []
for i in range(ngms):
	_, _, Sa = eqsig.sdof.pseudo_response_spectra(a1[i], dts[i], Trange, xi=0.05)
	Sa1.append(Sa)
	_, _, Sa = eqsig.sdof.pseudo_response_spectra(a2[i], dts[i], Trange, xi=0.05)
	Sa2.append(Sa)

#%% Plot the two directions' actual spectra
#for i in range(ngms):
#	plt.plot(Trange,Sa1[i],color=clr[i],linestyle=':',label='Dir 1 (Computed)')
#	plt.plot(Trange,Sa2[i],color=clr[i],linestyle='-.',label='Dir 2 (Computed)')

#%% Define a function to compute RotDxx
def get_RotDxx(acc1, acc2, dt, damp, xx, num_theta=100, num_T=100):
	theta = np.linspace(start=0, stop=np.pi, num=num_theta)
	Trange = np.logspace(start=-2, stop=1, num=num_T)
	
	Rot = np.zeros((num_theta,num_T))
	RotDxx = np.zeros((num_T))
	
	# Just do for one ground motion for now
	for j in range(num_theta):
		a = acc1*np.cos(theta[j])+acc2*np.sin(theta[j])
		_, _, Rot[j][:] = eqsig.sdof.pseudo_response_spectra(a, dt, Trange, xi=damp)
	RotDxx = np.percentile(Rot,xx,axis=0)
	
	return RotDxx, Trange

#%% Compute the RotD50 directly
RotD50_1 = np.zeros((ngms,100))
RotD50_2 = np.zeros((ngms,100))
for i in range(ngms):
	RotD50_1[i], Trange_RotD50 = get_RotDxx(a1[i], a2[i], dts[i], 0.005, 50)
	RotD50_2[i], Trange_RotD50 = get_RotDxx(a1[i], a2[i], dts[i], 0.05, 50)


#%% Plot the computed RotD50
for i in range(ngms):
	plt.plot(Trange_RotD50,RotD50_1[i],color=clr[i],linestyle='--',linewidth=1,label='RotD50 0.5% (Computed)')
	plt.plot(Trange_RotD50,RotD50_2[i],color=clr[i],linestyle='-.',linewidth=1,label='RotD50 5% (Computed)')
plt.legend()
plt.savefig('Spectra.pdf', bbox_inches='tight')
plt.close()
