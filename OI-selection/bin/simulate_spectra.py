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