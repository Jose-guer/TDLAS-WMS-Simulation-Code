# scanned_wms_sim.py

"""
This script demonstrates how to simulate scanned-WMS signals using the
calibration-free WMS model.

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

# constants
k = 10**3
M = 10**6
G = 10**9

# Simulation inputs
species = 'H2O'
wn0 = 7185.596  # or 6806.03

T = 1800.0
P = 1.0
Xrad = 0.15
L = 10.0

T_bg = 293.0
P_bg = 1.0
Xrad_bg = 0.001
L_bg = 10.0

fm = 35.0 * M
fscan = 500 * k
fs = 3 * G

a = 0.116
i0 = 0.1671
i2 = 2.48e-03
psi1 = 1.9356 * np.pi
psi2 = 4.4138 * np.pi
Idc = 1.2957

as_ = 0.06
psi_s = 1.411 * np.pi
i_s = 0.1661

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"

# Incident light intensity simulation
wm = 2.0 * np.pi * fm
ws = 2.0 * np.pi * fscan
t = np.arange(0.0, 1.0 / fscan + 1.0 / fs, 1.0 / fs)

wn_scan = wn0 + as_ * np.cos(ws * t)
I0_scan = 1.0 + i_s * np.cos(ws * t + psi_s)

wn_sim = wn_scan + a * np.cos(wm * t)
I0_sim = Idc * (1.0 + i_s * np.cos(ws * t + psi_s) + i0 * np.cos(wm * t + psi1) + i2 * np.cos(2.0 * wm * t + psi2))

i0 = i0 / I0_scan
i2 = i2 / I0_scan

# WMS signal simulation
delta_wn = 0.0005
wnL = float(np.min(wn_sim))
wnH = float(np.max(wn_sim))

X, bg_abs = HITRAN_Spectrum_Plot(T_bg, P_bg, Xrad_bg, L_bg, delta_wn, wnL, wnH, species, dat_path)
bg_abs_scan = np.interp(wn_scan, X, bg_abs)
bg_abs_sim = np.interp(wn_sim, X, bg_abs)
tau_bg_copy = np.exp(-bg_abs)

X, absorption = HITRAN_Spectrum_Plot(T, P, Xrad, L, delta_wn, wnL, wnH, species, dat_path)
absor_scan = np.interp(wn_scan, X, absorption)
absor_sim = np.interp(wn_sim, X, absorption)
tau = np.exp(-(absorption + bg_abs))

It_sim = I0_sim * np.exp(-(absor_sim + bg_abs_sim))

S_2f1f = np.zeros_like(wn_scan)
S_2f = np.zeros_like(wn_scan)
S_1f = np.zeros_like(wn_scan)
S_4f = np.zeros_like(wn_scan)
S_4f2f = np.zeros_like(wn_scan)
tsim = np.arange(-0.5 / fm, 0.5 / fm + 1.0 / fs, 1.0 / fs)

for j in range(wn_scan.size):
    wn0_j = wn_scan[j]
    wn_sim_new = wn0_j + a * np.cos(wm * tsim)
    tau_sim = np.interp(wn_sim_new, X, tau)
    tau_bg = np.interp(wn_sim_new, X, tau_bg_copy)

    X1f, Y1f, X2f, Y2f, X4f, Y4f, _ = wms_calibration_free(tsim, tau_sim, i0[j], psi1, i2[j], psi2)
    X1f_bg, Y1f_bg, X2f_bg, Y2f_bg, X4f_bg, Y4f_bg, _ = wms_calibration_free(tsim, tau_bg, i0[j], psi1, i2[j], psi2)

    R1f = np.sqrt(X1f**2 + Y1f**2)
    R1f_bg = np.sqrt(X1f_bg**2 + Y1f_bg**2)

    S_2f1f[j] = np.sqrt(((X2f / R1f) - (X2f_bg / R1f_bg))**2 + ((Y2f / R1f) - (Y2f_bg / R1f_bg))**2)
    S_2f[j] = np.sqrt((X2f - X2f_bg)**2 + (Y2f - Y2f_bg)**2)
    S_1f[j] = np.sqrt((X1f - X1f_bg)**2 + (Y1f - Y1f_bg)**2)
    S_4f[j] = np.sqrt((X4f - X4f_bg)**2 + (Y4f - Y4f_bg)**2)
    S_4f2f[j] = S_4f[j] / S_2f[j]

# Plots
fig, axs = plt.subplots(4, 1, figsize=(11.74, 7.55), sharex=True)

axs[0].plot(np.tile(It_sim, (1, 5)).T)
axs[0].legend([r'$I_t(t)$'], frameon=True)
axs[0].spines['top'].set_visible(False);  axs[0].spines['right'].set_visible(False)

axs[1].plot(np.tile(absor_scan, (1, 5)).T)
axs[1].legend([r'$\alpha(\nu)$'], frameon=True)
axs[1].spines['top'].set_visible(False);  axs[1].spines['right'].set_visible(False)

axs[2].plot(np.tile(S_2f1f, (1, 5)).T)
axs[2].plot(np.tile(S_4f * 3.0, (1, 5)).T)
axs[2].legend([r'WMS-2$f$/1$f$', r'WMS-4$f$ $\times 3$'], frameon=True)
axs[2].spines['top'].set_visible(False);  axs[2].spines['right'].set_visible(False)

axs[3].plot(np.tile(i0, (1, 5)).T)
axs[3].legend([r'$i_0(t)$'], frameon=True)
axs[3].spines['top'].set_visible(False);  axs[3].spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
