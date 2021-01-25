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

def simulate_spectra(random, nTrials,meanReq,covReq,stdevs,nGM,weights):
    # Import libraries
    import numpy as np
    from scipy.stats import skew
    
    # simulate response spectra from the target mean and covariance matrix
    devTotalSim =[]
    spettri=[]
    for j in np.arange(nTrials):
        SpectraSample=np.exp(random.multivariate_normal(meanReq,covReq,nGM))
        spettri.append(SpectraSample)
        # evaluate simulation
        sampleMeanErr=np.mean(np.log(SpectraSample),axis=0)-meanReq
        sampleStdErr=np.std(np.log(SpectraSample),axis=0)-stdevs
        sampleSkewnessErr=skew(np.log(SpectraSample),axis=0,bias=True)
        devTotalSim.append(weights[0]*sum(sampleMeanErr**2)+weights[1]**sum(sampleStdErr**2)+weights[2]*sum(sampleSkewnessErr**2))

    bestSample=devTotalSim.index(min(devTotalSim)) # find the simulated spectra that best match the target
    simulated_spectra = spettri[bestSample] #return the best set of simulations

    return simulated_spectra
