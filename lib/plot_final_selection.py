def plot_final_selection(name,nGM,TgtPer,sampleSmall,meanReq,stdevs,meanrecorded,output_folder):
	# Import libraries
	import numpy as np
	import matplotlib.pyplot as plt
	
	plt.figure()
	for i in np.arange(nGM):
		plt.loglog(TgtPer,np.exp(sampleSmall[i,:]),'g')
	plt.loglog(TgtPer,np.exp(meanReq),'r',label='target')
	plt.loglog(TgtPer,np.exp(meanReq+2*stdevs),'--r',label='target+2*sigma')
	plt.loglog(TgtPer,np.exp(meanReq-2*stdevs),'--r',label='target-2*sigma')
	plt.loglog(TgtPer,meanrecorded,'b',label='mean recorded')
	plt.xlabel('Period (s)')
	plt.ylabel('Spectral acceleration (g)')
	plt.xlim((0.1, np.max(TgtPer)))
	plt.legend()

	plt.grid(True)

	plt.savefig(output_folder+'/'+name+'/'+name+'_selection.png', bbox_inches='tight')

	plt.close()
