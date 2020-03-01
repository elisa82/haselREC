def akkar_correlation(t1, t2):
    """
    Read the period-dependent correlation coefficient matrix rho_H(T0,T) given in
    Akkar S, Sandikkaya MA, Ay BO (2014) Compatible ground-motion prediction equations
    for damping scaling factors and vertical-to-horizontal spectral amplitude ratios
    for the broader Europe region, Bull Earthquake Eng 12:517-547
    """

    import pathlib
    this = pathlib.Path(__file__)
    from scipy import interpolate

    coeff_table = []
    with open(this.parent / 'akkar_coeff_table.csv') as f:
        for row in f:
            coeff_table.append([float(col) for col in row.split(',')])

    periods = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.11, 0.12, 0.13, 0.14,
            0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3,
            0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.55, 0.6,
            0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
            1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]

    x=periods
    y=periods
    z=coeff_table
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
