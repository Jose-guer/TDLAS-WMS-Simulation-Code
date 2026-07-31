# Spectrum_Plot_tools.py

"""
Numerical tools for the simulation of and fitting to molecular absorption spectra

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu
"""

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

# Physical constants
h   = 6.6260680e-34      # J·s
kB  = 1.3806503e-23      # J/K
c   = 2.99792458e10      # cm/s
Tref= 296.0              # K

############################################################################

def HITRAN_Spectrum_Plot(T, P, Xrad, Lcm, dwn, wavenumLow, wavenumHigh, species, dat_path, struct_name='htd'):

    """
    This function simulates the absorbance spectrum and returns a wavenumber
    vector and the absorbance spectrum.

    Inputs
    T - Temperature [K]
    P - Pressure  [atm]
    Xrad - Radiating species mole fraction
    L - Path length [cm]
    dwn - wavenumber spacing [cm^-1]
    wavenumLow - Lower wavenumber bound [cm^-1]
    wavenumHigh - Upper wavenumber bound [cm^-1]
    species - Molecule name
    """

        # ------- Species selection -------
    if species == 'H2O':
        qfilename = dat_path + species + '/q_h2o_2020.txt'
        matfile   = dat_path + species + '/h2o_hitran2020.mat'
        MW_rad    = 18.01528
    elif species == 'CO':
        qfilename = dat_path + species + '/q_co_2020.txt'
        matfile   = dat_path + species + '/co_hitran2020.mat'
        MW_rad    = 27.994915
    elif species == 'CO2':
        qfilename = dat_path + species + '/q_co2_2020.txt'
        matfile   = dat_path + species + '/co2_hitran2020.mat'
        MW_rad    = 43.989830
    elif species == '13CO2':
        qfilename = dat_path + 'CO2/isotopes/q_13co2_2020.txt'
        matfile   = dat_path + 'CO2/isotopes/13co2_hitran2020.mat'
        MW_rad    = 44.993185
    elif species == 'CH4':
        qfilename = dat_path + species + '/q_ch4_2020.txt'
        matfile   = dat_path + species + '/ch4_hitran2020.mat'
        MW_rad    = 16.031300
    elif species == '13CH4':
        qfilename = dat_path + 'CH4/isotopes/q_13ch4_2020.txt'
        matfile   = dat_path + 'CH4/isotopes/13ch4_hitran2020.mat'
        MW_rad    = 17.034655
    elif species == 'NO':
        qfilename = dat_path + species + '/q_no_2020.txt'
        matfile   = dat_path + species + '/no_hitran2020.mat'
        MW_rad    = 29.9979895
    elif species == 'NO2':
        qfilename = dat_path + species + '/q_no2_2020.txt'
        matfile   = dat_path + species + '/no2_hitran2020.mat'
        MW_rad    = 45.992904
    elif species == 'H2O_HiTemp':
        qfilename = dat_path +  'H2O/q_h2o_2020.txt'
        matfile   = dat_path + 'H2O/H2O_2000_2015_HITEMP2010.mat'
        MW_rad    = 18.01528
    elif species == 'CO_HiTemp':
        qfilename = dat_path + 'CO/q_co_2020.txt'
        matfile   = dat_path + 'CO/CO_2000_2015_HITEMP2019.mat'
        MW_rad    = 27.994915
    elif species == 'CO2_HiTemp':
        qfilename = dat_path + 'CO2/q_co2_2020.txt'
        matfile   = dat_path + 'CO2/CO2_2000_2015_HITEMP2010.mat'
        MW_rad    = 43.989830
    else:
        raise ValueError(f"Unsupported species: {species}")
    # -------------------------------------------------------------

    # Load HITRAN data
    htd = load_hitran_struct(matfile, struct_name=struct_name)
    vo = htd['transitionWavenumber']

    # Window padding (±0.5 cm^-1) and selection
    wnL = wavenumLow - 0.5
    wnH = wavenumHigh + 0.5
    jj = (vo >= wnL) & (vo <= wnH)

    vo_vac     = htd['transitionWavenumber'][jj]
    STo        = htd['lineIntensity'][jj]
    Es         = htd['lowerStateEnergy'][jj]
    gamma_self = htd['selfBroadenedWidth'][jj]
    gamma_air  = htd['airBroadenedWidth'][jj]
    n_air      = htd['temperatureDependence'][jj]
    d_air      = htd['pressureShift'][jj]

    # n_self defaults and H2O special cases
    n_self = np.full_like(gamma_self, 0.75)
    if species == 'H2O' and vo_vac.size > 0:
        if vo_vac[0] > 7000:
            if (7185.596 > vo_vac.min()) and (7185.596 < vo_vac.max()):
                ind = np.argsort(np.abs(vo_vac - 7185.596))[:2]
                n_self[ind] = 0.53
        elif vo_vac[0] < 7000:
            if (6806.03 > vo_vac.min()) and (6806.03 < vo_vac.max()):
                ind = np.argsort(np.abs(vo_vac - 6806.03))[:2]
                n_self[ind] = 0.86

    # Partition function Q(T) 
    Q = fit_partition_q(qfilename)

    Xair = 1.0 - Xrad

    # Line strength with stimulated emission factor 
    ST = STo * (Q(Tref)/Q(T)) * np.exp(-h*c*Es/kB * (1.0/T - 1.0/Tref)) * \
         ((1.0 - np.exp(-h*c*vo_vac/(kB*T))) / (1.0 - np.exp(-h*c*vo_vac/(kB*Tref))))

    # Pressure-shifted line center
    vo_vac = vo_vac + P * Xair * d_air * (Tref / T) ** 0.96

    # Doppler and collisional HWHM 
    wG = (3.58115e-7) * vo_vac * np.sqrt(T / MW_rad)
    wL = P * (Xair * gamma_air * (Tref / T) ** n_air + Xrad * gamma_self * (Tref / T) ** n_self)

    # Wavenumber grid and accumulation
    X = np.arange(wnL, wnH + dwn/2.0, dwn)
    SVP = np.zeros_like(X)
    for j in range(vo_vac.size):
        VP = Voigt(X, vo_vac[j], wG[j], wL[j])
        SVP += ST[j] * VP

    # Absorbance
    absorption = SVP * ((P * Xrad) / (T * kB)) * Lcm * (0.101325)

    # Trim back to requested plotting range with ±0.02 cm^-1 pad (scalar bounds)
    ind = (X >= (wavenumLow - 0.02)) & (X <= (wavenumHigh + 0.02))
    return X[ind], absorption[ind]

############################################################################

def HITRAN_Spectrum_Plot_Alt(T, Pi, wLi, Lcm, dwn, wavenumLow, wavenumHigh, species, dat_path, struct_name='htd'):

    """
    This function simulates the absorbance spectrum and returns a wavenumber
    vector and the absorbance spectrum.

    Inputs
    T - Temperature [K]
    Pi - Partial Pressure  [atm]
    wLi - Radiating species collisional HWHM [cm^-1]
    L - Path length [cm]
    dwn - wavenumber spacing [cm^-1]
    wavenumLow - Lower wavenumber bound [cm^-1]
    wavenumHigh - Upper wavenumber bound [cm^-1]
    species - Molecule name
    """

    if species == 'H2O':
        qfilename = dat_path + species + '/q_h2o_2020.txt'
        matfile   = dat_path + species + '/h2o_hitran2020.mat'
        MW_rad    = 18.01528
    elif species == 'CO':
        qfilename = dat_path + species + '/q_co_2020.txt'
        matfile   = dat_path + species + '/co_hitran2020.mat'
        MW_rad    = 27.994915
    elif species == 'CO2':
        qfilename = dat_path + species + '/q_co2_2020.txt'
        matfile   = dat_path + species + '/co2_hitran2020.mat'
        MW_rad    = 43.989830
    elif species == '13CO2':
        qfilename = dat_path + 'CO2/isotopes/q_13co2_2020.txt'
        matfile   = dat_path + 'CO2/isotopes/13co2_hitran2020.mat'
        MW_rad    = 44.993185
    elif species == 'CH4':
        qfilename = dat_path + species + '/q_ch4_2020.txt'
        matfile   = dat_path + species + '/ch4_hitran2020.mat'
        MW_rad    = 16.031300
    elif species == '13CH4':
        qfilename = dat_path +  'CH4/isotopes/q_13ch4_2020.txt'
        matfile   = dat_path +  'CH4/isotopes/13ch4_hitran2020.mat'
        MW_rad    = 17.034655
    elif species == 'NO':
        qfilename = dat_path + species + '/q_no_2020.txt'
        matfile   = dat_path + species + '/no_hitran2020.mat'
        MW_rad    = 29.9979895
    elif species == 'NO2':
        qfilename = dat_path + species + '/q_no2_2020.txt'
        matfile   = dat_path + species + '/no2_hitran2020.mat'
        MW_rad    = 45.992904
    elif species == 'H2O_HiTemp':
        qfilename = dat_path +  'H2O/q_h2o_2020.txt'
        matfile   = dat_path + 'H2O/H2O_2000_2015_HITEMP2010.mat'
        MW_rad    = 18.01528
    elif species == 'CO_HiTemp':
        qfilename = dat_path + 'CO/q_co_2020.txt'
        matfile   = dat_path + 'CO/CO_2000_2015_HITEMP2019.mat'
        MW_rad    = 27.994915
    elif species == 'CO2_HiTemp':
        qfilename = dat_path + 'CO2/q_co2_2020.txt'
        matfile   = dat_path + 'CO2/CO2_2000_2015_HITEMP2010.mat'
        MW_rad    = 43.989830    
    else:
        raise ValueError(f"Unsupported species: {species}")

    htd = load_hitran_struct(matfile, struct_name=struct_name)
    vo = htd['transitionWavenumber']

    wnL = wavenumLow - 0.5
    wnH = wavenumHigh + 0.5
    jj = (vo >= wnL) & (vo <= wnH)

    vo_vac = htd['transitionWavenumber'][jj]
    STo    = htd['lineIntensity'][jj]
    Es     = htd['lowerStateEnergy'][jj]

    Q = fit_partition_q(qfilename)

    ST = STo * (Q(Tref)/Q(T)) * np.exp(-h*c*Es/kB * (1.0/T - 1.0/Tref)) * \
         ((1.0 - np.exp(-h*c*vo_vac/(kB*T))) / (1.0 - np.exp(-h*c*vo_vac/(kB*Tref))))

    wG = (3.58115e-7) * vo_vac * np.sqrt(T / MW_rad)
    wL = wLi * np.ones_like(wG)

    X = np.arange(wnL, wnH + dwn/2.0, dwn)
    SVP = np.zeros_like(X)
    for j in range(vo_vac.size):
        VP = Voigt(X, vo_vac[j], wG[j], wL[j])
        SVP += ST[j] * VP

    absorption = SVP * (Pi / (T * kB)) * Lcm * 0.101325

    ind = (X >= (wavenumLow - 0.02)) & (X <= (wavenumHigh + 0.02))
    return X[ind], absorption[ind]

############################################################################

def SpectrumPlot(Trange, Prange, Xrange, Lcm, wavenum, species, dat_path, struct_name='htd'):
    
    """
    This function loops over the range of temperature, pressures, and mole fraction, and make a plot.

    Inputs
    Trange - Temperature range [K]
    Prange - Pressure range [atm]
    Xrange - Radiating species mole fraction range
    L - Path length [cm]
    wavenum - wavenumber vector [cm^-1]
    species - Molecule name
    """

    if species == 'H2O':
        qfilename = dat_path + species + '/q_h2o_2020.txt'
        matfile   = dat_path + species + '/h2o_hitran2020.mat'
        MW_rad    = 18.01528
    elif species == 'CO':
        qfilename = dat_path + species + '/q_co_2020.txt'
        matfile   = dat_path + species + '/co_hitran2020.mat'
        MW_rad    = 27.994915
    elif species == 'CO2':
        qfilename = dat_path + species + '/q_co2_2020.txt'
        matfile   = dat_path + species + '/co2_hitran2020.mat'
        MW_rad    = 43.989830
    elif species == '13CO2':
        qfilename = dat_path + 'CO2/isotopes/q_13co2_2020.txt'
        matfile   = dat_path + 'CO2/isotopes/13co2_hitran2020.mat'
        MW_rad    = 44.993185
    elif species == 'CH4':
        qfilename = dat_path + species + '/q_ch4_2020.txt'
        matfile   = dat_path + species + '/ch4_hitran2020.mat'
        MW_rad    = 16.031300
    elif species == '13CH4':
        qfilename = dat_path + 'CH4/isotopes/q_13ch4_2020.txt'
        matfile   = dat_path + 'CH4/isotopes/13ch4_hitran2020.mat'
        MW_rad    = 17.034655
    elif species == 'NO':
        qfilename = dat_path + species + '/q_no_2020.txt'
        matfile   = dat_path + species + '/no_hitran2020.mat'
        MW_rad    = 29.9979895
    elif species == 'NO2':
        qfilename = dat_path + species + '/q_no2_2020.txt'
        matfile   = dat_path + species + '/no2_hitran2020.mat'
        MW_rad    = 45.992904
    elif species == 'H2O_HiTemp':
        qfilename = dat_path +  'H2O/q_h2o_2020.txt'
        matfile   = dat_path + 'H2O/H2O_2000_2015_HITEMP2010.mat'
        MW_rad    = 18.01528
    elif species == 'CO_HiTemp':
        qfilename = dat_path + 'CO/q_co_2020.txt'
        matfile   = dat_path + 'CO/CO_2000_2015_HITEMP2019.mat'
        MW_rad    = 27.994915
    elif species == 'CO2_HiTemp':
        qfilename = dat_path + 'CO2/q_co2_2020.txt'
        matfile   = dat_path + 'CO2/CO2_2000_2015_HITEMP2010.mat'
        MW_rad    = 43.989830    
    else:
        raise ValueError(f"Unsupported species: {species}")

    htd = load_hitran_struct(matfile, struct_name=struct_name)
    vo_all = htd['transitionWavenumber']

    wnL = float(np.min(wavenum)) - 0.5
    wnH = float(np.max(wavenum)) + 0.5
    dwn = float(wavenum[1] - wavenum[0])

    jj = (vo_all >= wnL) & (vo_all <= wnH)

    vo_base    = htd['transitionWavenumber'][jj]
    STo        = htd['lineIntensity'][jj]
    Es         = htd['lowerStateEnergy'][jj]
    gamma_self = htd['selfBroadenedWidth'][jj]
    gamma_air  = htd['airBroadenedWidth'][jj]
    n_air      = htd['temperatureDependence'][jj]
    d_air      = htd['pressureShift'][jj]

    n_self = np.full_like(gamma_self, 0.75)
    if species == 'H2O' and vo_base.size > 0:
        if vo_base[0] > 7000:
            if (7185.596 > vo_base.min()) and (7185.596 < vo_base.max()):
                ind = np.argsort(np.abs(vo_base - 7185.596))[:2]
                n_self[ind] = 0.53
        elif vo_base[0] < 7000:
            if (6806.03 > vo_base.min()) and (6806.03 < vo_base.max()):
                ind = np.argsort(np.abs(vo_base - 6806.03))[:2]
                n_self[ind] = 0.86

    Q = fit_partition_q(qfilename)

    plt.figure()
    for Xrad in np.atleast_1d(Xrange):
        Xair = 1.0 - Xrad
        for T in np.atleast_1d(Trange):
            for P in np.atleast_1d(Prange):
                ST = STo * (Q(Tref)/Q(T)) * np.exp(-h*c*Es/kB * (1.0/T - 1.0/Tref)) * \
                     ((1.0 - np.exp(-h*c*vo_base/(kB*T))) / (1.0 - np.exp(-h*c*vo_base/(kB*Tref))))

                vo = vo_base + P * Xair * d_air * (Tref / T) ** 0.96
                wG = (3.58115e-7) * vo * np.sqrt(T / MW_rad)
                wL = P * (Xair * gamma_air * (Tref / T) ** n_air + Xrad * gamma_self * (Tref / T) ** n_self)

                X = np.arange(wnL, wnH + dwn/2.0, dwn)
                SVP = np.zeros_like(X)
                for j in range(vo.size):
                    VP = Voigt(X, vo[j], wG[j], wL[j])
                    SVP += ST[j] * VP

                absorption = SVP * ((P * Xrad) / (T * kB)) * Lcm * 0.101325
                ind = (X >= np.min(wavenum)) & (X <= np.max(wavenum))
                plt.plot(X[ind], absorption[ind])

    z = 1
    legStr = []
    for Xrad in np.atleast_1d(Xrange):
        for T in np.atleast_1d(Trange):
            for P in np.atleast_1d(Prange):
                legStr.append(f"T = {T:.0f}K, P = {np.round(P,2)}atm, X = {np.round(Xrad,2)}")
                z += 1

    plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
    plt.ylabel(r'$\alpha(\nu)$')
    plt.legend(legStr, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1)
    plt.gcf().set_size_inches(7.32, 5.95)
    plt.tight_layout()
    plt.show()

############################################################################

def load_hitran_struct(mat_path: str, struct_name: str = 'htd'):

    d = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    s = d[struct_name] 
    out = {}

    fields = [
        'transitionWavenumber',
        'lineIntensity',
        'lowerStateEnergy',
        'selfBroadenedWidth',
        'airBroadenedWidth',
        'temperatureDependence',
        'pressureShift',
    ]
    for f in fields:
        out[f] = np.atleast_1d(np.array(getattr(s, f), dtype=float)).squeeze()
    return out

############################################################################

def fit_partition_q(poly_path: str):

    dat = np.loadtxt(poly_path)
    Tcol = dat[:, 0]
    Qcol = dat[:, 1]
    coeffs = np.polyfit(Tcol, Qcol, 6)  # degree-6 fit
    def Q(T):
        return np.polyval(coeffs, T)
    return Q

############################################################################

def Gaussian(v, vo, wG):
    """
    v  : wavenumber array
    vo : linecenter 
    wG : FWHM
    """
    v = np.asarray(v, dtype=np.float64)
    y = 2.0 / wG * np.sqrt(np.log(2.0) / np.pi) * np.exp(-4.0 * np.log(2.0) * (v - vo)**2 / wG**2)
    return y

############################################################################

def Lorentzian(v, vo, wL):
    """
    v  : wavenumber array
    vo : linecenter
    wL : FWHM
    """
    v = np.asarray(v, dtype=np.float64)
    y = (1.0 / (2.0 * np.pi)) * wL / ((v - vo)**2 + (wL / 2.0)**2)
    return y

############################################################################

def Voigt(wavenumber_array, center_line, width_gauss, width_lorentz):

    """
    Voigt profile calculation

    INPUT ARGUMENTS
    wavenumberArray - array 1*N of wavenumbers 
    centerLine - transition line center
    widthGauss - Gaussian / Doppler HWHM
    widthLorentz - Lorentzian / Collisional HWHM

	OUTPUT
    Voigt profile - array 1*N of intensities

    For more details on algorithm see the publication:
    F. Schreier: Optimized Implementations of Rational Approximations for the Voigt ane Complex Error Function. 
    J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 10101025, doi 10.1016/j.jqsrt.2010.12.010, 2011. 

    """
    wavenumber_array = np.asarray(wavenumber_array, dtype=np.float64)
    width_gauss = float(width_gauss)
    width_lorentz = float(width_lorentz)

    sqrt_log2 = np.sqrt(np.log(2.0))
    x = sqrt_log2 * (wavenumber_array - center_line) / width_gauss
    y = sqrt_log2 * (width_lorentz / width_gauss)

    w = complexErrorFunction(x, y)
    V = np.sqrt(np.log(2.0)/np.pi) / width_gauss * np.real(w)
    return V

############################################################################

def complexErrorFunction(x, y):

    x = np.asarray(x, dtype=np.float64).ravel()   # column-like vector
    y = np.asarray(y, dtype=np.float64)           # allow scalar or array
    if y.ndim == 0:
        y = np.full_like(x, y)
    else:
        y = np.asarray(y, dtype=np.float64).ravel()
        if y.size == 1 and x.size > 1:
            y = np.full_like(x, y.item())
        elif y.size != x.size:
            raise ValueError("x and y must be broadcastable to the same length")

    # constants
    half = 0.5
    recSqrtPi = 1.0/np.sqrt(np.pi)
    s15 = 15.0
    l = 4.1195342878142354  # sqrt(n/sqrt(2)) with n=24

    a = np.array([
        -1.5137461654527820e-10,  4.9048215867870488e-09,  1.3310461806370372e-09,
        -3.0082822811202271e-08, -1.9122258522976932e-08,  1.8738343486619108e-07,
         2.5682641346701115e-07, -1.0856475790698251e-06, -3.0388931839840047e-06,
         4.1394617248575527e-06,  3.0471066083243790e-05,  2.4331415462641969e-05,
        -2.0748431511424456e-04, -7.8166429956142650e-04, -4.9364269012806686e-04,
         6.2150063629501763e-03,  3.3723366855316413e-02,  1.0838723484566792e-01,
         2.6549639598807689e-01,  5.3611395357291292e-01,  9.2570871385886788e-01,
         1.3948196733791203e+00,  1.8562864992055408e+00,  2.1978589365315417e+00
    ], dtype=np.float64)

    w = np.zeros(x.shape, dtype=np.complex128)

    # region test like MATLAB: s = |x| + y
    s = np.abs(x) + y
    region_I = ( (y > s15) | (s > s15) )
    # Humliček region I
    if np.any(region_I):
        t = y[region_I] - 1j*x[region_I]
        w[region_I] = (recSqrtPi * t) / (half + t*t)

    # Weideman region (complement)
    ind = ~region_I
    if np.any(ind):
        xm = x[ind]; ym = y[ind]
        recLmZ = 1.0 / (l + ym - 1j*xm)             # 1/(l+y - i x)
        t = (l - ym + 1j*xm) * recLmZ               # (l - y + i x)/(l + y - i x)
        # Horner evaluation of the polynomial in t with coeffs a[0..23]
        poly = a[0]
        for k in range(1, a.size):
            poly = poly * t + a[k]
        w[ind] = recLmZ * (recSqrtPi + 2.0*recLmZ*poly)

    return w

############################################################################

def multi_line_voigt_fit(vo, A, wL, T, wavenum, abs_meas):
    """
    Reconstructs the absorption spectrum using the input parameters and
    returns the residual.

    vo       : vector of line centers
    A        : areas of Voigt profiles for lines in vo
    wL       : collisional broadening HWHM of lines
    T        : temperature
    wavenum  : wavenumber vector for measured absorption spectrum
    abs_meas : measured absorption spectrum
    """

    MW_rad = 18.01528  # H2O
    wG = (3.58115e-7) * vo * np.sqrt(T / MW_rad)

    VP = np.zeros((len(wavenum), len(vo)))

    for j in range(len(vo)):
        VP[:, j] = Voigt(wavenum, vo[j], wG[j], wL[j])

    abs_rec = A * VP

    if len(vo) == 1:
        abs_rec = abs_rec[:, 0]
    else:
        abs_rec = np.sum(abs_rec, axis=1)

    Resid = abs_meas - abs_rec

    return Resid

############################################################################

def cosinefit_lsq(t, dat, a, phi, b, f):
    t = np.asarray(t, dtype=float)
    dat = np.asarray(dat, dtype=float)
    return dat - (a * np.cos(2.0 * np.pi * f * t + phi) + b)

######################################################################


def getHITRAN(T, Xrad, wavenumLow, wavenumHigh, threshold, species, wn_units, S_units, dat_path):

    # ---- physical constants ----
    h  = 6.6260680e-34    # J·s
    kB = 1.3806503e-23    # J/K
    c  = 2.99792458e10    # cm/s
    Tref = 296.0

    # ---- file selection (paths are appended to dat_path) ----
    if species == 'H2O':
        qfilename = dat_path + 'H2O/q_h2o_2020.txt'
        matfile   = dat_path + 'H2O/h2o_hitran2020.mat'
    elif species == 'CO':
        qfilename = dat_path + 'CO/q_co_2020.txt'
        matfile   = dat_path + 'CO/co_hitran2020.mat'
    elif species == 'CO2':
        qfilename = dat_path + 'CO2/q_co2_2020.txt'
        matfile   = dat_path + 'CO2/co2_hitran2020.mat'
    elif species == '13CO2':
        qfilename = dat_path + 'CO2/isotopes/q_13co2_2020.txt'
        matfile   = dat_path + 'CO2/isotopes/13co2_hitran2020.mat'
    elif species == 'CH4':
        qfilename = dat_path + 'CH4/q_ch4_2020.txt'
        matfile   = dat_path + 'CH4/ch4_hitran2020.mat'
    elif species == '13CH4':
        qfilename = dat_path + 'CH4/isotopes/q_13ch4_2020.txt'
        matfile   = dat_path + 'CH4/isotopes/13ch4_hitran2020.mat'
    elif species == 'NO':
        qfilename = dat_path + 'NO/q_no_2020.txt'
        matfile   = dat_path + 'NO/no_hitran2020.mat'
    elif species == 'NO2':
        qfilename = dat_path + 'NO2/q_no2_2020.txt'
        matfile   = dat_path + 'NO2/no2_hitran2020.mat'
    elif species == 'H2O_HiTemp':
        qfilename = dat_path +  'H2O/q_h2o_2020.txt'
        matfile   = dat_path + 'H2O/H2O_2000_2015_HITEMP2010.mat'
        MW_rad    = 18.01528
    elif species == 'CO_HiTemp':
        qfilename = dat_path + 'CO/q_co_2020.txt'
        matfile   = dat_path + 'CO/CO_2000_2015_HITEMP2019.mat'
        MW_rad    = 27.994915
    elif species == 'CO2_HiTemp':
        qfilename = dat_path + 'CO2/q_co2_2020.txt'
        matfile   = dat_path + 'CO2/CO2_2000_2015_HITEMP2010.mat'
        MW_rad    = 43.989830    
    else:
        raise ValueError(f"Unsupported species: {species}")

    # ---- load HITRAN struct ----
    htd = load_hitran_struct(matfile, struct_name='htd')
    vo_all = htd['transitionWavenumber']

    # selection in requested range
    ind = (vo_all >= wavenumLow) & (vo_all <= wavenumHigh)

    vo_vac     = htd['transitionWavenumber'][ind]
    So         = htd['lineIntensity'][ind]
    gamma_self = htd['selfBroadenedWidth'][ind]
    gamma_air  = htd['airBroadenedWidth'][ind]
    n_air      = htd['temperatureDependence'][ind]
    Es         = htd['lowerStateEnergy'][ind]
    d_air      = htd['pressureShift'][ind]

    # ---- partition function Q(T) ----
    pf = np.loadtxt(qfilename)
    coeffs = np.polyfit(pf[:, 0], pf[:, 1], 6)
    def Q(Tin):
        return np.polyval(coeffs, Tin)

    # ---- line strength S(T) ----
    if np.isclose(T, Tref):
        S = So.copy()
    else:
        S = So * (Q(Tref) / Q(T)) * np.exp(-h * c * Es / kB * (1.0 / T - 1.0 / Tref)) * \
            ((1.0 - np.exp(-h * c * vo_vac / (kB * T))) /
             (1.0 - np.exp(-h * c * vo_vac / (kB * Tref))))

    # S units conversion to cm^-2/atm 
    if S_units == 1:
        S = S * (7.34e21) / T
        So = So * (7.34e21) / Tref

    # Convert frequency micrometer units    
    if wn_units == 0:
        vo_vac = 10**4 / vo_vac

    XS = Xrad * S

    ind = XS >= threshold

    vo_vac = vo_vac[ind]
    Es = Es[ind]
    XS = XS[ind]
    So = So[ind]
    gamma_air = gamma_air[ind]
    gamma_self = gamma_self[ind]
    n_air = n_air[ind]
    d_air = d_air[ind]

    return vo_vac, Es, XS, So, gamma_air, gamma_self, n_air, d_air
