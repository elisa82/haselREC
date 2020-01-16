def akkar_correlation(t1, t2):
    """
    Read the period-dependent correlation coefficient matrix rho_H(T0,T) given in
    Akkar S, Sandikkaya MA, Ay BO (2014) Compatible ground-motion prediction equations
    for damping scaling factors and vertical-to-horizontal spectral amplitude ratios
    for the broader Europe region, Bull Earthquake Eng 12:517-547
    """
    # Import libraries
    from openquake.hazardlib.gsim.mgmpe import akkar_coeff_table as act
    from scipy import interpolate

    x=act.periods
    y=act.periods
    z=act.coeff_table
    f=interpolate.interp2d(x, y, z, kind='linear')
    interpolated_value=t1,t2,f(t1,t2)

    return interpolated_value[2][0]


def baker_jayaram_correlation(t1, t2):
    """
    NOTE: subroutine taken from: https://usgs.github.io/shakemap/shakelib

    Produce inter-period correlation for any two spectral periods.

    Based upon:
    Baker, J.W. and Jayaram, N., "Correlation of spectral acceleration
    values from NGA ground motion models," Earthquake Spectra, (2007).

    Args:
        t1, t2 (float):
            The two periods of interest.

    Returns:
        rho (float): The predicted correlation coefficient

    """
    # Import libraries
    import numpy as np
    
    t_min = min(t1, t2)
    t_max = max(t1, t2)

    c1 = 1.0 - np.cos(np.pi / 2.0 - np.log(t_max / max(t_min, 0.109)) * 0.366)

    if t_max < 0.2:
        c2 = 1.0 - 0.105 * (1.0 - 1.0 / (1.0 + np.exp(100.0 * t_max - 5.0)))*(t_max - t_min) / (t_max - 0.0099)
    else:
        c2 = 0

    if t_max < 0.109:
        c3 = c2
    else:
        c3 = c1

    c4 = c1 + 0.5 * (np.sqrt(c3) - c3) * (1.0 + np.cos(np.pi * t_min / 0.109))

    if t_max <= 0.109:
        rho = c2
    elif t_min > 0.109:
        rho = c1
    elif t_max < 0.2:
        rho = min(c2, c4)
    else:
        rho = c4

    return rho
