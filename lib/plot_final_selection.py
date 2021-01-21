def plot_final_selection(name,lbl,nGM,T_CS,sampleSmall,meanReq,stdevs,output_folder):
	# Import libraries
    import numpy as np
    import matplotlib.pyplot as plt

    meanrecorded=np.mean(np.exp(sampleSmall),axis=0)
    meanrecorded_p2sigma = np.percentile(np.exp(sampleSmall),50+34.1+13.6,axis=0)
    meanrecorded_n2sigma = np.percentile(np.exp(sampleSmall),50-34.1-13.6,axis=0)
    meanrecorded_eps =  (np.log(meanrecorded_p2sigma)-np.log(meanrecorded_n2sigma))/(2*1.96)
    indexes = [i for i,x in enumerate(T_CS) if x > 0]
    stdevs_gt0=[]
    T_CS_gt0=[]
    meanReq_gt0=[]
    meanrecorded_eps_gt0=[]
    meanrecorded_p2sigma_gt0=[]
    meanrecorded_n2sigma_gt0=[]
    meanrecorded_gt0=[]
    for i in indexes:
        T_CS_gt0.append(T_CS[i])
        meanrecorded_eps_gt0.append(meanrecorded_eps[i])
        meanrecorded_p2sigma_gt0.append(meanrecorded_p2sigma[i])
        meanrecorded_n2sigma_gt0.append(meanrecorded_n2sigma[i])
        meanrecorded_gt0.append(meanrecorded[i])
        stdevs_gt0.append(stdevs[i])
        meanReq_gt0.append(meanReq[i])
    stdevs_gt0=np.asarray(stdevs_gt0)
    T_CS_gt0=np.asarray(T_CS_gt0)
    meanReq_gt0=np.asarray(meanReq_gt0)
    meanrecorded_eps_gt0=np.asarray(meanrecorded_eps_gt0)
    meanrecorded_p2sigma_gt0=np.asarray(meanrecorded_p2sigma_gt0)
    meanrecorded_n2sigma_gt0=np.asarray(meanrecorded_n2sigma_gt0)
    meanrecorded_gt0=np.asarray(meanrecorded_gt0)

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
        sampleSmall_gt0=[]
        for j in indexes:
            sampleSmall_gt0.append(sampleSmall[i,j])
        plt.loglog(T_CS_gt0,np.exp(sampleSmall_gt0),'g', linewidth=.5)
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0),'r',label='CMS', linewidth=1.0)
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0+2*stdevs_gt0),'--r',label=r'CMS $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0-2*stdevs_gt0),'--r', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_gt0,'k',label='Selected', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_p2sigma_gt0,'--k',label=r'Selected $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_n2sigma_gt0,'--k', linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel(lbl+' [g]')
    plt.xlim(min(T_CS_gt0),max(T_CS_gt0))
    plt.ylim(1e-2,1e1)
    plt.yscale('log')
    plt.xscale('log')
    #number=int(name[9])+1
    #plt.title('site '+str(number))
    plt.grid(True)
    plt.legend()
    plt.savefig(output_folder+'/'+name+'/'+name+'_spectra_gms.pdf', bbox_inches='tight')
    plt.close()

        # Spectra
    plt.figure(figsize=(1.5*2.36,2.36))
    plt.rcParams.update({'font.size': 8})
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0),'r',label='CMS', linewidth=1.0)
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0+2*stdevs_gt0),'--r',label=r'CMS $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(T_CS_gt0,np.exp(meanReq_gt0-2*stdevs_gt0),'--r', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_gt0,'k',label='Selected', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_p2sigma_gt0,'--k',label=r'Selected $\pm 2\sigma$', linewidth=1.0)
    plt.loglog(T_CS_gt0,meanrecorded_n2sigma_gt0,'--k', linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel(lbl+' [g]')
    plt.xlim(min(T_CS_gt0),max(T_CS_gt0))
    plt.ylim(1e-2,1e1)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_folder+'/'+name+'/'+name+'_spectra.pdf', bbox_inches='tight')
    plt.close()

	# Dispersion
    plt.figure(figsize=(1.5*2.36,2.36))
    plt.rcParams.update({'font.size': 8})
    plt.plot(T_CS_gt0,stdevs_gt0,'r',label='CMS', linewidth=1.0)
    plt.plot(T_CS_gt0,meanrecorded_eps_gt0,'k',label='Selected', linewidth=1.0)
    plt.xlabel('Period [s]')
    plt.ylabel('Dispersion')
    plt.xlim(min(T_CS_gt0),max(T_CS_gt0))
    plt.ylim(0,1)
    plt.xscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_folder+'/'+name+'/'+name+'_dispersion.pdf', bbox_inches='tight')
    plt.close()
