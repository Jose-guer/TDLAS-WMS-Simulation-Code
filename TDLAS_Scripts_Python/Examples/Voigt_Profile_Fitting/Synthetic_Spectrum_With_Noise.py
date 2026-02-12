# Synthetic_Spectrum_With_Noise.py

"""
%This script is used to simulate synthetic absorption spectra with noise.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import matplotlib.pyplot as plt
from Utilities.Spectrum_Plot_tools import HITRAN_Spectrum_Plot


# -------- SIMULATION INPUTS --------
species = 'H2O'
T = 1800.0      # K
P = 1.0         # atm
Xrad = 0.15
L = 10.0        # cm
dwn = 0.001     # cm^-1
dat_path = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"

# ----- Window 1 -----
wnL = 7185.3
wnH = 7185.8
x1 = np.linspace(wnL, wnH, 500)
X1, absor1 = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

# ----- Window 2 -----
wnL = 6806.03 - 0.3
wnH = 6806.03 + 0.2
x2 = np.linspace(wnL, wnH, 500)
X2, absor2 = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

# Interpolate and add noise
y1 = np.interp(x1, X1, absor1)
y2 = np.interp(x2, X2, absor2)

A = 0.005  # noise amplitude
noise1 = A * (np.random.rand(y1.size) - 0.5)
noise2 = A * (np.random.rand(y2.size) - 0.5)

y1 = y1 + noise1
y2 = y2 + noise2

dat = {
    "x1": x1, "absor1": y1,
    "x2": x2, "absor2": y2,
    "P": P, "T": T, "X": Xrad, "L": L
}

# Plots
plt.figure(figsize=(7.32, 5.95))
plt.plot(x1, y1)
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.tight_layout()

plt.figure(figsize=(7.32, 5.95))
plt.plot(x2, y2)
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.tight_layout()

plt.show()
