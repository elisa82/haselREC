def plot_final_selection(name,lbl,nGM,TgtPer,T_CS,sampleSmall,meanReq,stdevs,meanrecorded,meanrecorded_p2sigma,meanrecorded_n2sigma,meanrecorded_eps,output_folder):
	# Import libraries
	import numpy as np
	import matplotlib.pyplot as plt

#	plt.figure()
#	for i in np.arange(nGM):
#		plt.loglog(TgtPer,np.exp(sampleSmall[i,:]),'g')
#	plt.loglog(TgtPer,np.exp(meanReq),'r',label='target')
#	plt.loglog(TgtPer,np.exp(meanReq+2*stdevs),'--r',label='target+2*sigma')
#	plt.loglog(TgtPer,np.exp(meanReq-2*stdevs),'--r',label='target-2*sigma')
#	plt.loglog(TgtPer,meanrecorded,'b',label='mean recorded')
#	plt.xlabel('Period (s)')
#	plt.ylabel('Spectral acceleration (g)')
#	plt.xlim((0.01, np.max(TgtPer)))
#	plt.legend()
#	plt.grid(True)
#	plt.savefig(output_folder+'/'+name+'/'+name+'_selection.png', bbox_inches='tight')
#	plt.close()

	# Spectra with ground motions
	plt.figure(figsize=(1.5*2.36,2.36))
	plt.rcParams.update({'font.size': 8})
	for i in np.arange(nGM):
		plt.loglog(TgtPer,np.exp(sampleSmall[i,:]),'g', linewidth=.5)
	plt.plot(T_CS,np.exp(meanReq),'r',label='CMS', linewidth=1.0)
	plt.plot(T_CS,np.exp(meanReq+2*stdevs),'--r',label=r'CMS $\pm 2\sigma$', linewidth=1.0)
	plt.plot(T_CS,np.exp(meanReq-2*stdevs),'--r', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded,'k',label='Selected', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded_p2sigma,'--k',label=r'Selected $\pm 2\sigma$', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded_n2sigma,'--k', linewidth=1.0)
	plt.xlabel('Period [s]')
	plt.ylabel('Intensity Measure - '+lbl)
	plt.xlim(min(T_CS),max(T_CS))
	plt.ylim(1e-2,1e1)
	plt.yscale('log')
	plt.xscale('log')
	plt.grid(True)
	plt.legend()
	plt.savefig(output_folder+'/'+name+'/'+name+'_spectra_gms.pdf', bbox_inches='tight')
	plt.close()

	# Spectra
	plt.figure(figsize=(1.5*2.36,2.36))
	plt.rcParams.update({'font.size': 8})
	plt.plot(T_CS,np.exp(meanReq),'r',label='CMS', linewidth=1.0)
	plt.plot(T_CS,np.exp(meanReq+2*stdevs),'--r',label=r'CMS $\pm 2\sigma$', linewidth=1.0)
	plt.plot(T_CS,np.exp(meanReq-2*stdevs),'--r', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded,'k',label='Selected', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded_p2sigma,'--k',label=r'Selected $\pm 2\sigma$', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded_n2sigma,'--k', linewidth=1.0)
	plt.xlabel('Period [s]')
	plt.ylabel('Intensity Measure - '+lbl)
	plt.xlim(min(T_CS),max(T_CS))
	plt.ylim(1e-2,1e1)
	plt.yscale('log')
	plt.xscale('log')
	plt.grid(True)
	plt.legend()
	plt.savefig(output_folder+'/'+name+'/'+name+'_spectra.pdf', bbox_inches='tight')
	plt.close()

	# Epsilon
	plt.figure(figsize=(1.5*2.36,2.36))
	plt.rcParams.update({'font.size': 8})
	plt.plot(T_CS,stdevs,'r',label='CMS', linewidth=1.0)
	plt.plot(TgtPer,meanrecorded_eps,'k',label='Selected', linewidth=1.0)
	plt.xlabel('Period [s]')
	plt.ylabel('Dispersion')
	plt.xlim(min(T_CS),max(T_CS))
	plt.ylim(0,1)
	plt.xscale('log')
	plt.grid(True)
	plt.legend()
	plt.savefig(output_folder+'/'+name+'/'+name+'_dispersion.pdf', bbox_inches='tight')
	plt.close()
