"""
This script loads in the sample data and determines the intensity
characterization parameters for WMS. First a cosine wave at 1f is fit to
the baseline and then a cosine wave at 2f is fit to the residual. The DC
normalized intensity modulation amplitudes are printed at the end.
When the nonlinear solver fails to converge, generally chaning the inital
phase guess solves the problem. The phase of the 1st and 2nd order fits
are inputs into the etalon_post_processing.m script where the actual
relative phase shifts will be determined.

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

from Utilities.Spectrum_Plot_tools import cosinefit_lsq

# ---------------- Load data ----------------
mat = loadmat("etalon_dat.mat", squeeze_me=True, struct_as_record=False)
dat = mat["dat"]
if isinstance(dat, np.ndarray):
    dat = dat.item()

t        = np.asarray(dat.t).squeeze()
baseline = np.asarray(dat.baseline).squeeze()
etalon   = np.asarray(dat.etalon).squeeze()
f        = float(dat.f)

rgb = (0.0, 0.5, 1.0)

plt.figure()
plt.plot(t, baseline, 'k')
plt.plot(t, etalon, color=rgb)
plt.legend(['Baseline', 'Etalon'], frameon=True)
plt.xlabel('time [s]')
plt.ylabel('Amplitude [V]')

# ---------------- First-order (1f) fit ----------------
i0_init = 0.5 * (baseline.max() - baseline.min())
phi_init = 1.4 * np.pi
b_init = float(np.mean(baseline))
x0 = np.array([i0_init, phi_init, b_init], dtype=float)

lb = np.array([0.8 * i0_init, 0.0, 0.8 * b_init])
ub = np.array([1.2 * i0_init, 2.0 * np.pi, 1.2 * b_init])

def resid_1f(x):
    i0, phi, b = x
    return cosinefit_lsq(t, baseline, i0, phi, b, f)

sol1 = least_squares(resid_1f, x0, bounds=(lb, ub), method="trf")
i0, phi, b = sol1.x

plt.figure()
plt.plot(t, baseline, 'k')
plt.plot(t, i0 * np.cos(2.0 * np.pi * f * t + phi) + b, 'r--', linewidth=2)
plt.xlabel('time (s)')
plt.ylabel('Amplitude [V]')
plt.legend(['Baseline', '1st-order fit'], frameon=True)

# ---------------- Second-order (2f) fit on residual ----------------
res = baseline - (i0 * np.cos(2.0 * np.pi * f * t + phi) + b)

i2_init = 0.5 * (res.max() - res.min())
phi2_init = 0.2 * np.pi
b2 = 0.0  # fixed to 0 in original routine
x0_2 = np.array([i2_init, phi2_init], dtype=float)

lb2 = np.array([0.25 * i2_init, 0.0])
ub2 = np.array([1.5 * i2_init, 2.0 * np.pi])

def resid_2f(x):
    i2, phi2 = x
    return cosinefit_lsq(t, res, i2, phi2, b2, 2.0 * f)

sol2 = least_squares(resid_2f, x0_2, bounds=(lb2, ub2), method="trf")
i2, phi2 = sol2.x

plt.figure()
plt.plot(t, res, 'k')
plt.plot(t, i2 * np.cos(2.0 * np.pi * 2.0 * f * t + phi2), 'r-', linewidth=2)
plt.xlabel('time (s)')
plt.ylabel('Amplitude [V]')
plt.legend(['Residual', '2nd-order fit'], frameon=True)

# ---------------- Summary plots ----------------
Idc = b
scan_fit = Idc + i0 * np.cos(2.0 * np.pi * f * t + phi) + i2 * np.cos(2.0 * np.pi * (2.0 * f) * t + phi2)

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].plot(t, baseline)
axs[0].plot(t, scan_fit, '--')
axs[0].legend(['Baseline', 'Fit'], frameon=True)
axs[0].set_ylabel('Amplitude [V]')

axs[1].plot(t, (scan_fit - baseline) * 100.0 / baseline)
axs[1].set_ylabel(r'$\%$')
axs[1].legend(['Residual'], frameon=True)
axs[1].set_xlabel('time (s)')

plt.tight_layout()
plt.show()

# ---------------- Results ----------------
i0_dc = i0 / Idc
i2_dc = i2 / Idc
phi1_over_pi = phi / np.pi
phi2_over_pi = phi2 / np.pi

print("\nIntensity characterization:")
print(f"i0        : {i0_dc:.6g} (DC-normalized)")
print(f"phi1/pi   : {phi1_over_pi:.6g}")
print(f"i2        : {i2_dc:.6g} (DC-normalized)")
print(f"phi2/pi   : {phi2_over_pi:.6g}")
print(f"Idc       : {Idc:.6g}")
