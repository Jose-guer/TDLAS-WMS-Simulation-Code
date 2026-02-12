
"""
This script implements a multi-line Voigt profile fitting routine to infer
gas properties from noisy synthetic absorption spectra. The algorithm
implements a nonlinear fitting routine where A, wL, and vo are varried.
The fitting routine is initiated by guessing an inital temperature. The
fitting routine will loop until the inferred temperature from the best-fit
parameters is within a specified tolerance of the current iterations
guess.

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu

"""
import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # project root
os.chdir(Path(__file__).resolve().parent)                      # use script folder for data files

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import least_squares

from Utilities.Spectrum_Plot_tools import Voigt, multi_line_voigt_fit

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
P  = float(dat.P)
T  = float(dat.T)
X  = float(dat.X)
L  = float(dat.L)

# -------------------- Plot raw inputs --------------------
plt.figure()
ax1 = plt.subplot(2, 1, 1)
plt.plot(x1, y1, 'b', markersize=1)
plt.ylabel(r'$\alpha(\nu)$')

ax2 = plt.subplot(2, 1, 2)
plt.plot(x2, y2, 'r', markersize=1)
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.tight_layout()

# -------------------- Line-strength ratio model --------------------
# using the sum of the doublets linestrength values in the HITRAN database
Es1  = 1045.058   # cm^-1
STo1 = 0.01962    # cm^-2/atm  (pressure-normalized reference strength)

Es2  = 3291.15    # cm^-1
STo2 = 6.55e-07   # cm^-2/atm

Tsim = np.arange(1000.0, 2500.0 + 100.0, 100.0)

# constants
h   = 6.6260680e-34      # J·s
kB  = 1.3806503e-23      # J/K
c   = 2.99792458e10      # cm/s
Tref= 296.0              # K
MW_rad = 18.01528        # g/mol (H2O)

# Partition function Q(T) from file
dat_path = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"
qfilename = Path(dat_path) / "H2O" / "q_h2o_2020.txt"
pq = np.loadtxt(qfilename)
coeffs = np.polyfit(pq[:, 0], pq[:, 1], 6)
Q = lambda TT: np.polyval(coeffs, TT)

# Pressure-normalized line strengths S(T) [cm^-2/atm]
ST1 = STo1 * (Q(Tref) / Q(Tsim)) * (Tref / Tsim) * np.exp(-h * c * Es1 / kB * (1.0 / Tsim - 1.0 / Tref))
ST2 = STo2 * (Q(Tref) / Q(Tsim)) * (Tref / Tsim) * np.exp(-h * c * Es2 / kB * (1.0 / Tsim - 1.0 / Tref))
R    = ST1 / ST2

# -------------------- Initial guesses --------------------
Tguess = 1000.0

vo1 = np.array([7185.4, 7185.596])
A1  = np.array([0.0015, 0.0150])
wG1 = (3.58115e-7) * vo1 * np.sqrt(Tguess / MW_rad)
wL1 = 2.0 * wG1

vo2 = np.array([6805.81, 6806.033, 6806.12])
A2  = np.array([0.0010, 0.0047, 0.0008])
wG2 = (3.58115e-7) * vo2 * np.sqrt(Tguess / MW_rad)
wL2 = 2.0 * wG2

# -------------------- Iterative fit loop --------------------
errT = 100.0  # percent

while errT > 0.1:
    # --- Fit window 1 (two lines) ---
    X0 = np.hstack([vo1, A1, wL1])
    lb = np.hstack([vo1 - 0.005, 0.5 * np.hstack([A1, wL1])])
    ub = np.hstack([vo1 + 0.005, 2.0 * np.hstack([A1, wL1])])

    def resid1(X):
        v_centers = np.array([X[0], X[1]])
        areas     = np.array([X[2], X[3]])
        widthsL   = np.array([X[4], X[5]])
        return multi_line_voigt_fit(v_centers, areas, widthsL, Tguess, x1, y1)

    sol1 = least_squares(resid1, X0, bounds=(lb, ub), method="trf")
    vo1 = np.array([sol1.x[0], sol1.x[1]])
    A1  = np.array([sol1.x[2], sol1.x[3]])
    wL1 = np.array([sol1.x[4], sol1.x[5]])
    Area1 = A1[1]

    # --- Fit window 2 (three lines) ---
    X0 = np.hstack([vo2, A2, wL2])
    lb = np.hstack([vo2 - 0.005, 0.5 * np.hstack([A2, wL2])])
    ub = np.hstack([vo2 + 0.005, 2.0 * np.hstack([A2, wL2])])

    def resid2(X):
        v_centers = np.array([X[0], X[1], X[2]])
        areas     = np.array([X[3], X[4], X[5]])
        widthsL   = np.array([X[6], X[7], X[8]])
        return multi_line_voigt_fit(v_centers, areas, widthsL, Tguess, x2, y2)

    sol2 = least_squares(resid2, X0, bounds=(lb, ub), method="trf")
    vo2 = np.array([sol2.x[0], sol2.x[1], sol2.x[2]])
    A2  = np.array([sol2.x[3], sol2.x[4], sol2.x[5]])
    wL2 = np.array([sol2.x[6], sol2.x[7], sol2.x[8]])
    Area2 = A2[1]

    # --- Update temperature from area ratio ---
    ratio = Area1 / Area2
    Tmeas = np.interp(ratio, R, Tsim)
    errT  = abs(Tmeas - Tguess) / Tguess * 100.0

    # update Doppler widths with new temperature
    wG1 = (3.58115e-7) * vo1 * np.sqrt(Tmeas / MW_rad)
    wG2 = (3.58115e-7) * vo2 * np.sqrt(Tmeas / MW_rad)

    Tguess = Tmeas

# -------------------- Partial pressure / mole fraction --------------------
S1   = np.interp(Tmeas, Tsim, ST1)
XH2O = Area1 / (S1 * P * L)

# -------------------- Reconstruct spectra from best-fit --------------------
# Window 1
VP1 = np.zeros((x1.size, vo1.size))
for j in range(vo1.size):
    VP1[:, j] = Voigt(x1, vo1[j], wG1[j], wL1[j])
abs_fit1 = (A1 * VP1).sum(axis=1)

# Window 2
VP2 = np.zeros((x2.size, vo2.size))
for j in range(vo2.size):
    VP2[:, j] = Voigt(x2, vo2[j], wG2[j], wL2[j])
abs_fit2 = (A2 * VP2).sum(axis=1)

# -------------------- Plots --------------------
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(x1, y1, 'k*', markersize=1)
plt.plot(x1, abs_fit1, 'b', linewidth=2)
plt.ylabel(r'$\alpha(\nu)$')

plt.subplot(2, 1, 2)
plt.plot(x2, y2, 'k*', markersize=1)
plt.plot(x2, abs_fit2, 'r', linewidth=2)
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.tight_layout()

# -------------------- Report --------------------

errT_pct = (T - Tmeas) / T * 100.0
errX_pct = (X - XH2O) / X * 100.0
print("\n\nResults:")
print(f"T_meas [K] : {Tmeas:8.2f}    err T [%] : {errT_pct:7.3f}")
print(f"X_meas     : {XH2O:8.5f}    err X [%] : {errX_pct:7.3f}")
print("\n\n")

plt.show()

