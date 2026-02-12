"""
% This script loads in the etalon signal and labelled peaks and determines
% the modulation depth and the relative-phase shifts if the phase of the
% 1st and 2nd order fits to the intensity signal are specified. 
%
% The etalon peaks here have already been found and labeled depending on
% whether the wavenumber is increasing or decreasing. Each peak is evenly
% spaced by the FSR of the etalon so the peaks only need to be "counted"
% adding 1 when the wavenumber increases, subtracting 1 when the wavenumber
% decreases, or adding 0 after a frequency reversal event. The data here 
% illistrates the general process.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

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

t       = np.asarray(dat.t).squeeze()
etalon  = np.asarray(dat.etalon).squeeze()
f       = float(dat.f)

tpeaks      = np.asarray(dat.tpeaks).squeeze()
ypeaks      = np.asarray(dat.ypeaks).squeeze()
peak_label  = np.asarray(dat.peak_label).squeeze()
fsr         = float(dat.fsr)

baseline = np.asarray(dat.baseline).squeeze()

# From intensity post-processing
phi1 = 0.9715 * np.pi
phi2 = 0.8672 * np.pi

# ---------------- Compute relative frequency from labels ----------------
rel_freq = peak_label * fsr

# Initial guess and bounds: [a, phi, b]
a0   = 0.5 * (rel_freq.max() - rel_freq.min())
phi0 = np.pi
b0   = float(np.mean(rel_freq))
x0   = np.array([a0, phi0, b0], dtype=float)

lb = np.array([0.0, 0.0, -1.0])
ub = np.array([0.5, 2.0 * np.pi, 1.0])

def resid_cosine(x):
    a, phi, b = x
    return rel_freq - (a * np.cos(2.0 * np.pi * f * tpeaks + phi) + b)

sol = least_squares(resid_cosine, x0, bounds=(lb, ub), method="trf")
a, phi, b = sol.x
a = float(np.round(a, 3))

# ---------------- Relative phase shifts ----------------
psi_1 = (phi1 - phi) / np.pi
psi_2 = (phi2 - 2.0 * phi) / np.pi

psi_1 = (psi_1 + 10.0) % 2.0
psi_2 = (psi_2 + 10.0) % 2.0

# ---------------- Plots ----------------
rgb = (0.0, 0.5, 1.0)

baseline_n = baseline - np.mean(baseline)
baseline_n = baseline_n / np.max(np.abs(baseline_n))
m = 10

plt.figure()
plt.plot(t[::m], baseline_n[::m], 'k', linewidth=2)
plt.plot(t, np.cos(2.0 * np.pi * f * t + phi), color=rgb, linewidth=2)
plt.title('FM/IM Phase shift')
plt.xlabel('time (s)')
plt.legend(['Intensity', 'Wavenumber'], frameon=True)

fig, axs = plt.subplots(3, 1, sharex=True)
axs[0].plot(t, etalon, color=rgb)
axs[0].plot(tpeaks, ypeaks, 'k*')
axs[0].set_ylabel('Amplitude [V]')
axs[0].legend(['Etalon', 'Peaks'], frameon=True)

axs[1].plot(tpeaks, peak_label, 'k*')
axs[1].grid(True, which='both')
axs[1].set_ylabel('Label')
axs[1].legend(['Peaks'], frameon=True)
axs[1].set_ylim([-10, 10])

axs[2].plot(tpeaks, rel_freq, 'k*')
t_fit = np.linspace(tpeaks.min(), tpeaks.max(), 2000)
axs[2].plot(t_fit, a * np.cos(2.0 * np.pi * f * t_fit + phi) + b, 'r-', linewidth=2)
axs[2].set_xlabel('time (s)')
axs[2].set_ylabel(r'rel. freq. $[\mathrm{cm^{-1}}]$')
axs[2].legend(['Peaks', 'Best Fit'], frameon=True)
axs[2].set_ylim([-0.2, 0.15])

plt.tight_layout()
plt.show()

# ---------------- Results ----------------
print(f"modulation_depth = {a} cm^-1")
print(f"psi1 = {psi_1}*pi")
print(f"psi2 = {psi_2}*pi")
