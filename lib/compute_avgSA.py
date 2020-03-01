def compute_avgSA(avg_periods,sctx, rctx, dctx, bgmpe, corr_type):
    # Import libraries
    from openquake.hazardlib import imt, const
    import numpy as np
    from lib.im_correlation import baker_jayaram_correlation
    from lib.im_correlation import akkar_correlation
    
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

def compute_rho_avgSA(per,avg_periods,sctx,rctx,dctx,stddvs_avgsa, bgmpe, corr_type):
    # Import libraries
    from openquake.hazardlib import imt, const
    from im_correlation import baker_jayaram_correlation
    from im_correlation import akkar_correlation
    
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
