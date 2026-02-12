
"""
This script simulates the two-color WMS-2f/1f ratio and WMS-2f/1f signals 
at a specified pressure and laser specific tuning parameters

Hardcoded values are taken from [1] and Fig. 5 is reproduced.

[1] Schwartz, Charles J., et al. "Near-MHz temperature and H2O measurements 
in post-detonation fireballs of 25 g hemispherical explosives using 
scanned-wavelength-modulation spectroscopy." 
Applied Optics 62.6 (2023): 1598-1609.

Linecenter pressure-shift can be commented out in
"HITRAN_Spectrum_Plot.m" line 128, when generating scanned-WMS tables. 

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu

"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt

from Utilities.Spectrum_Plot_tools import HITRAN_Spectrum_Plot
from Utilities.WMS_Sim_tools import wms_calibration_free

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"

# ************ SIMULATION INPUTS ***************

M = 10**6
G = 10**9

fs = 3 * G  # Sampling frequency [Hz]

# Test gas details
species = 'H2O'

P = 1.0  # [atm]
T = np.arange(800.0, 2400.0 + 100.0, 200.0)  # 800:100:2400
H2O = np.array([0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31])  
L = 12.6  # [cm]

# Background details
Tbg = 296.0
H2O_bg = 0.0
Lbg = 0.0  # [cm]

# Laser 1 - Absorption and WMS-2f/1f Signal

wn0 = 7185.596
fmod = 35.0 * M  # Hz
a = 0.113        # cm^-1
i0 = 0.160
i2 = 3.18e-3
psi1 = 6.08
psi2 = 6.92

t = np.arange(-0.5 / fmod, 0.5 / fmod + 1.0 / fs, 1.0 / fs)
w = 2.0 * np.pi * fmod
wn_sim = wn0 + a * np.cos(w * t)

delta_wn = 0.001
wnL = np.min(wn_sim)
wnH = np.max(wn_sim)

# Background absorption
X, bg_abs = HITRAN_Spectrum_Plot(Tbg, P, H2O_bg, Lbg, delta_wn, wnL, wnH, species, dat_path)
tau_bg_copy = np.exp(-bg_abs)
# 1D linear interpolation
tau_bg = np.interp(wn_sim, X, tau_bg_copy)

# Background WMS signal
X1f_bg, Y1f_bg, X2f_bg, Y2f_bg, X4f_bg, Y4f_bg, _ = wms_calibration_free(t, tau_bg, i0, psi1, i2, psi2)

S2f1f_L1 = np.zeros((H2O.size, T.size))

for i in range(H2O.size):
    for j in range(T.size):
        # combustion absorption
        _, absorption = HITRAN_Spectrum_Plot(T[j], P, H2O[i], L, delta_wn, wnL, wnH, species, dat_path)
        
        tau = np.exp(-(absorption + bg_abs))
        tau_sim = np.interp(wn_sim, X, tau)

        # Absorption Signal
        X1f, Y1f, X2f, Y2f, X4f, Y4f, _ = wms_calibration_free(t, tau_sim, i0, psi1, i2, psi2)

        R1f = np.sqrt(X1f**2 + Y1f**2)
        R1f_bg = np.sqrt(X1f_bg**2 + Y1f_bg**2)

        S2f1f_L1[i, j] = np.sqrt(((X2f / R1f) - (X2f_bg / R1f_bg))**2 + ((Y2f / R1f) - (Y2f_bg / R1f_bg))**2)

# ---------------------------------------------
# Laser 2 - Absorption and WMS-2f/1f Signal
# ---------------------------------------------

wn0 = 6806.03
fmod = 45.5 * M  # Hz

a = 0.166       # cm^-1
i0 = 0.175
i2 = 8.394e-3
psi1 = 5.88
psi2 = 4.22

t = np.arange(-0.5 / fmod, 0.5 / fmod + 1.0 / fs, 1.0 / fs)
w = 2.0 * np.pi * fmod
wn_sim = wn0 + a * np.cos(w * t)

delta_wn = 0.001
wnL = np.min(wn_sim)
wnH = np.max(wn_sim)

# Background absorption
X, bg_abs = HITRAN_Spectrum_Plot(Tbg, P, H2O_bg, Lbg, delta_wn, wnL, wnH, species, dat_path)
tau_bg_copy = np.exp(-bg_abs)
tau_bg = np.interp(wn_sim, X, tau_bg_copy)

# Background WMS signal
X1f_bg, Y1f_bg, X2f_bg, Y2f_bg, X4f_bg, Y4f_bg, _ = wms_calibration_free(t, tau_bg, i0, psi1, i2, psi2)

S2f1f_L2 = np.zeros((H2O.size, T.size))

for i in range(H2O.size):
    for j in range(T.size):
        # combustion absorption
        _, absorption = HITRAN_Spectrum_Plot(T[j], P, H2O[i], L, delta_wn, wnL, wnH, species, dat_path)
        
        tau = np.exp(-(absorption + bg_abs))
        tau_sim = np.interp(wn_sim, X, tau)

        # Absorption Signal
        X1f, Y1f, X2f, Y2f, X4f, Y4f, _ = wms_calibration_free(t, tau_sim, i0, psi1, i2, psi2)

        R1f = np.sqrt(X1f**2 + Y1f**2)
        R1f_bg = np.sqrt(X1f_bg**2 + Y1f_bg**2)

        S2f1f_L2[i, j] = np.sqrt(((X2f / R1f) - (X2f_bg / R1f_bg))**2 + ((Y2f / R1f) - (Y2f_bg / R1f_bg))**2)

# ---------------------------------------------
# Plots
# ---------------------------------------------

legStr = [rf"$X_{{\mathrm{{H_2O}}}}$ = {x:.2f}" for x in H2O]

plt.figure()
for ii in range(H2O.size):
    plt.plot(T, S2f1f_L2[ii, :] / S2f1f_L1[ii, :], label=legStr[ii])

plt.xlabel('T [K]')
plt.ylabel('Two-Color WMS-2f/1f Ratio')
plt.legend(loc='upper left', frameon=True)
plt.xlim([np.min(T), np.max(T)])
plt.gcf().set_size_inches(7.65, 6.28)
plt.tight_layout()

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.show()