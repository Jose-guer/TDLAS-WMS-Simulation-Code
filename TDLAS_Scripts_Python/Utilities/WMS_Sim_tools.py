# WMS_Sim_tools.py

"""
Numerical tools for the simulation and extraction of WMS Harmonics

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu
"""

import numpy as np
from scipy.signal import butter, lfilter

############################################################################

def wms_calibration_free(t, tau, i0, psi1, i2, psi2):
    """
    Returns X and Y components for WMS-1f, -2f, and -4f signals and Fourier
    coefficients Hk using the calibration-free WMS model.

    Inputs
    t - time vector
    tau - transmittance = exp(-absorption)
    i0, i2, - first and second order DC normalized intensity modulation
    amplitudes
    psi1, psi2, - first and second order relative phase-shifts

    Rieker, Gregory B., Jay B. Jeffries, and Ronald K. Hanson. 
    Calibration-free wavelength-modulation spectroscopy for measurements of 
    gas temperature and concentration in harsh environments.
    Applied optics 48.29 (2009): 5546-5560.

    """
    H0 = wms_Hk(t, tau, 0)
    H1 = wms_Hk(t, tau, 1)
    H2 = wms_Hk(t, tau, 2)
    H3 = wms_Hk(t, tau, 3)
    H4 = wms_Hk(t, tau, 4)
    H5 = wms_Hk(t, tau, 5)
    H6 = wms_Hk(t, tau, 6)

    H = np.array([H1, H2, H3, H4, H5, H6])

    # 1f
    X1f = H1 + i0 * (H0 + H2/2.0) * np.cos(psi1) + (i2/2.0) * (H1 + H3) * np.cos(psi2)
    Y1f = -i0 * (H0 - H2/2.0) * np.sin(psi1) + (i2/2.0) * (H1 - H3) * np.sin(psi2)

    # 2f
    X2f = H2 + (i0/2.0) * (H1 + H3) * np.cos(psi1) + i2 * (H0 + H4/2.0) * np.cos(psi2)
    Y2f = -(i0/2.0) * (H1 - H3) * np.sin(psi1) + i2 * (H0 - H4/2.0) * np.sin(psi2)

    # 4f
    X4f = H4 + (i0/2.0) * (H5 + H3) * np.cos(psi1) + (i2/2.0) * (H6 + H2) * np.cos(psi2)
    Y4f = (i0/2.0) * (H5 - H3) * np.sin(psi1) + (i2/2.0) * (H6 - H2) * np.sin(psi2)

    return X1f, Y1f, X2f, Y2f, X4f, Y4f, H

############################################################################

def wms_Hk(t, tau, k):
    """
    Returns the k-th Fourier coefficient.
    t : time vector (seconds)
    tau : transmittance = exp(-absorption)
    k : harmonic number (integer >= 0)
    """
    t = np.asarray(t, dtype=float)
    tau = np.asarray(tau, dtype=float)

    T = t[-1] - t[0]
    w = 2.0 * np.pi / T

    if k == 0:
        Hk = (1.0 / T) * np.trapz(tau, t)
    else:
        Hk = (2.0 / T) * np.trapz(tau * np.cos(k * w * t), t)

    return Hk

############################################################################

def filter_butter(inputSignal, fs, fc, N):
    """
    Filters inputSignal using an Nth-order low-pass Butterworth filter.
    fs : sampling frequency (Hz)
    fc : cutoff frequency (Hz)
    N  : filter order (integer)
    """
    inputSignal = np.asarray(inputSignal, dtype=float)
    Wn = fc / (fs / 2.0)  # normalize by Nyquist
    b, a = butter(N, Wn, btype='low', analog=False, output='ba')
    filteredSignal = lfilter(b, a, inputSignal)
    return filteredSignal

############################################################################

def WMS_harmonic(S, t, fs, fc, fm, harmonic, filterOrder):
    """
    Digital lock-in to extract the n-th harmonic of signal S.
    S : detector signal
    t : time vector (seconds)
    fs : sampling frequency (Hz)
    fc : LPF cutoff frequency (Hz)
    fm : modulation frequency (Hz)
    harmonic : n (1, 2, 3, ...)
    filterOrder : Butterworth filter order
    """
    S = np.asarray(S, dtype=float)
    t = np.asarray(t, dtype=float)

    n = harmonic
    w = 2.0 * np.pi * fm

    # X-component
    X = S * np.cos(n * w * t)
    X_nf = filter_butter(X, fs, fc, filterOrder)

    # Y-component
    Y = S * np.sin(n * w * t)
    Y_nf = filter_butter(Y, fs, fc, filterOrder)

    # n-th harmonic magnitude
    WMS_nf = np.sqrt(X_nf**2 + Y_nf**2)

    return WMS_nf, X_nf, Y_nf

############################################################################
