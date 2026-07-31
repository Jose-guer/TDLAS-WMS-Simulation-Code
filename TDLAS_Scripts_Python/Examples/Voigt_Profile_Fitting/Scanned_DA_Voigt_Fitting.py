"""
This script implements a multi-line Voigt profile fitting routine to infer
gas properties from noisy synthetic absorption spectra. The algorithm
implements a nonlinear fitting routine where A, wL, and vo are varied.
The fitting routine is initiated by guessing an initial temperature. The
fitting routine will loop until the inferred temperature from the best-fit
parameters is within a specified tolerance of the current iteration's
guess.

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu

"""

import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # project root
os.chdir(Path(__file__).resolve().parent)                    # use script folder for data files

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import least_squares

from Utilities.Spectrum_Plot_tools import Voigt, multi_line_voigt_fit


def fit_residual(X, n_lines, T, wavenum, abs_meas):
    vo = X[:n_lines]
    A = X[n_lines:2*n_lines]
    wL = X[2*n_lines:3*n_lines]

    return multi_line_voigt_fit(vo, A, wL, T, wavenum, abs_meas)


# -------------------- Load synthetic data --------------------
mat = loadmat("data_noise.mat", squeeze_me=True, struct_as_record=False)
# mat = loadmat("data_noise_free.mat", squeeze_me=True, struct_as_record=False)
dat = mat["dat"]

# Unwrap the MATLAB struct
if isinstance(dat, np.ndarray):
    dat = dat.item()

x1 = np.asarray(dat.x1).squeeze()
y1 = np.asarray(dat.absor1).squeeze()
x2 = np.asarray(dat.x2).squeeze()
y2 = np.asarray(dat.absor2).squeeze()

P = float(dat.P)
T = float(dat.T)
X = float(dat.X)
L = float(dat.L)


# -------------------- Plot raw inputs --------------------
plt.figure()

plt.subplot(2, 1, 1)
plt.plot(x1, y1)
plt.ylabel(r'$\alpha(\nu)$')

plt.subplot(2, 1, 2)
plt.plot(x2, y2, 'r')
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')

plt.tight_layout()


# -------------------- Infer gas properties --------------------

# Using the sum of the doublet linestrength values in the HITRAN database
Es1 = 1045.058     # cm^-1
STo1 = 0.01962     # cm^-2/atm

Es2 = 3291.15      # cm^-1
STo2 = 6.55e-07    # cm^-2/atm

Tsim = np.arange(1000.0, 2500.0 + 100.0, 100.0)  # K

qfilename = Path("/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/H2O/q_h2o_2020.txt")

h = 6.6260680e-34      # Planck constant, J-s
kB = 1.3806503e-23     # Boltzmann constant, J/K
c = 2.99792458e10      # Speed of light, cm/s
Tref = 296.0           # K
MW_rad = 18.01528

# Partition function
dat = np.loadtxt(qfilename)
pn = np.polyfit(dat[:, 0], dat[:, 1], 6)
Q = lambda T: np.polyval(pn, T)

# Pressure-normalized definition of the line-strength function
ST1 = STo1 * (Q(Tref) / Q(Tsim)) * (Tref / Tsim) * np.exp(-h * c * Es1 / kB * (1.0 / Tsim - 1.0 / Tref))
ST2 = STo2 * (Q(Tref) / Q(Tsim)) * (Tref / Tsim) * np.exp(-h * c * Es2 / kB * (1.0 / Tsim - 1.0 / Tref))
R = ST1 / ST2


# -------------------- Initial guesses --------------------

Tguess = 1000.0

vo1 = np.array([7185.4, 7185.596])
A1 = np.array([0.0015, 0.015])
wG1 = (3.58115e-7) * vo1 * np.sqrt(Tguess / MW_rad)
wL1 = 2.0 * wG1

vo2 = np.array([6805.81, 6806.033, 6806.12])
A2 = np.array([0.001, 0.0047, 0.0008])
wG2 = (3.58115e-7) * vo2 * np.sqrt(Tguess / MW_rad)
wL2 = 2.0 * wG2


errT = 100.0

while errT > 0.1:

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X0 = np.array([vo1[0], vo1[1],
                   A1[0], A1[1],
                   wL1[0], wL1[1]])

    lb = np.array([vo1[0] - 0.005, vo1[1] - 0.005,
                   A1[0] * 0.5, A1[1] * 0.5,
                   wL1[0] * 0.5, wL1[1] * 0.5])

    ub = np.array([vo1[0] + 0.005, vo1[1] + 0.005,
                   A1[0] * 2.0, A1[1] * 2.0,
                   wL1[0] * 2.0, wL1[1] * 2.0])

    fit_param1 = least_squares(fit_residual,X0,bounds=(lb, ub),args=(2, Tguess, x1, y1),method="trf").x

    vo1 = np.array([fit_param1[0], fit_param1[1]])
    A1 = np.array([fit_param1[2], fit_param1[3]])
    Area1 = A1[1]
    wL1 = np.array([fit_param1[4], fit_param1[5]])

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X0 = np.array([vo2[0], vo2[1], vo2[2],
                   A2[0], A2[1], A2[2],
                   wL2[0], wL2[1], wL2[2]])

    lb = np.array([vo2[0] - 0.005, vo2[1] - 0.005, vo2[2] - 0.005,
                   A2[0] * 0.5, A2[1] * 0.5, A2[2] * 0.5,
                   wL2[0] * 0.5, wL2[1] * 0.5, wL2[2] * 0.5])

    ub = np.array([vo2[0] + 0.005, vo2[1] + 0.005, vo2[2] + 0.005,
                   A2[0] * 2.0, A2[1] * 2.0, A2[2] * 2.0,
                   wL2[0] * 2.0, wL2[1] * 2.0, wL2[2] * 2.0])

    fit_param2 = least_squares(fit_residual,X0,bounds=(lb, ub),args=(3, Tguess, x2, y2),method="trf").x

    vo2 = np.array([fit_param2[0], fit_param2[1], fit_param2[2]])
    A2 = np.array([fit_param2[3], fit_param2[4], fit_param2[5]])
    Area2 = A2[1]
    wL2 = np.array([fit_param2[6], fit_param2[7], fit_param2[8]])

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Tmeas = np.interp(Area1 / Area2, R, Tsim)
    Tmeas = np.interp(Area1 / Area2, R[::-1], Tsim[::-1])

    errT = abs(Tmeas - Tguess) / Tguess * 100.0

    wG1 = (3.58115e-7) * vo1 * np.sqrt(Tmeas / MW_rad)
    wG2 = (3.58115e-7) * vo2 * np.sqrt(Tmeas / MW_rad)

    # Update guessed temperature
    Tguess = Tmeas


# -------------------- Species mole fraction --------------------

S1 = np.interp(Tmeas, Tsim, ST1)
XH2O = Area1 / (S1 * P * L)

# -------------------- Reconstruct absorption spectra -----------

# 7185.596 cm^-1
VP = np.zeros((len(x1), len(vo1)))

for j in range(len(vo1)):
    VP[:, j] = Voigt(x1, vo1[j], wG1[j], wL1[j])

abs_fit1 = A1 * VP
abs_fit1 = np.sum(abs_fit1, axis=1)


# 6806.03 cm^-1
VP = np.zeros((len(x2), len(vo2)))

for j in range(len(vo2)):
    VP[:, j] = Voigt(x2, vo2[j], wG2[j], wL2[j])

abs_fit2 = A2 * VP
abs_fit2 = np.sum(abs_fit2, axis=1)


# -------------------- Plot results --------------------

plt.figure()

plt.subplot(2, 1, 1)
plt.plot(x1, y1, 'k*', markersize=1)
plt.plot(x1, abs_fit1, 'b', linewidth=2)
plt.ylabel(r'$\alpha(\nu)$')

plt.subplot(2, 1, 2)
plt.plot(x2, y2, 'k*', markersize=1)
plt.plot(x2, abs_fit2, 'r', linewidth=2)
plt.ylabel(r'$\alpha(\nu)$')
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')

plt.tight_layout()


# -------------------- Print results --------------------

errT = (T - Tmeas) / T * 100.0
errX = (X - XH2O) / X * 100.0

A = np.array([Tmeas, errT, XH2O, errX])

print("\n\nResults:")
print(f"T_meas [K] : {Tmeas:8.2f}    err T [%] : {errT:7.3f}")
print(f"X_meas     : {XH2O:8.5f}    err X [%] : {errX:7.3f}")
print("\n\n")

plt.show()
