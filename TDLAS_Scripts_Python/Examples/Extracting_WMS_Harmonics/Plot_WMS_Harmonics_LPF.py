# plot_WMS_Harmonics_LPF.py

"""
This script loads in the aligned background and absorption intensity
signals and extracts the background subtracted WMS harmonics. A power
spectral density plot highlighting the hamonics and side-bands is also
computed.

written by Jose Guerrero, University of Michigan - Aerospace Department 
joseguer@umich.edu

"""
import sys, os
from pathlib import Path

# make Utilities import work (project root)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
# make relative files load from THIS script’s folder
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from Utilities.WMS_Sim_tools import WMS_harmonic

k = 10**3
M = 10**6

raw = loadmat('raw_data.mat', squeeze_me=True, struct_as_record=False)['raw']
fscan = float(raw.fscan)
fm1 = float(raw.fm1)
fm2 = float(raw.fm2)

dat = loadmat('aligned_data.mat', squeeze_me=True)['dat']
t = dat[:, 0]
sig = dat[:, 1]
bg_sig = dat[:, 2]

fs = 1.0 / (t[1] - t[0])
N = sig.size

filterOrder = 5
fc = 2.2 * fscan

# -------- Power Spectral Density of Raw Signal --------
plt.figure()
plt.plot(t * 1e6, bg_sig)
plt.plot(t * 1e6, sig)
ax = plt.gca()
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.xlabel('time [$\\mu$s]')
plt.ylabel('Signal [V]')
plt.legend(['Background', 'Absorption Sig.'], frameon=True)
plt.tight_layout()

PSD_signal = np.abs(np.fft.fft(sig)) ** 2
f_PSD = np.arange(N) / N * fs / 1e6  # MHz

ind = (f_PSD > (fm1 - 2.2 * fscan) / M) & (f_PSD < (fm1 + 2.2 * fscan) / M)
PSD_1f = PSD_signal[ind]
freq_1f = f_PSD[ind]

ind = (f_PSD > (2 * fm1 - 2.2 * fscan) / M) & (f_PSD < (2 * fm1 + 2.2 * fscan) / M)
PSD_2f = PSD_signal[ind]
freq_2f = f_PSD[ind]

ind = (f_PSD > (fm2 - 2.2 * fscan) / M) & (f_PSD < (fm2 + 2.2 * fscan) / M)
PSD_1f_2 = PSD_signal[ind]
freq_1f_2 = f_PSD[ind]

ind = (f_PSD > (2 * fm2 - 2.2 * fscan) / M) & (f_PSD < (2 * fm2 + 2.2 * fscan) / M)
PSD_2f_2 = PSD_signal[ind]
freq_2f_2 = f_PSD[ind]

plt.figure()
plt.semilogy(f_PSD, PSD_signal, 'k')
plt.semilogy(freq_1f, PSD_1f, 'b')
plt.semilogy(freq_2f, PSD_2f, 'b')
plt.semilogy(freq_1f_2, PSD_1f_2, 'r')
plt.semilogy(freq_2f_2, PSD_2f_2, 'r')
ax = plt.gca()
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.xlabel('Frequency (MHz)')
plt.ylabel('PSD ($\\mathrm{V^2/Hz}$)')
plt.xlim([30, 93])
plt.legend(['', '$\\nu_o$ = 1391.7nm', '', '$\\nu_o$ = 1469.3nm'], frameon=True)
plt.tight_layout()

# ------------------ Extract WMS Harmonics ------------------
# Laser 1
_, X1f, Y1f = WMS_harmonic(sig, t, fs, fc, fm1, 1, filterOrder)
_, X1f_bg, Y1f_bg = WMS_harmonic(bg_sig, t, fs, fc, fm1, 1, filterOrder)

_, X2f, Y2f = WMS_harmonic(sig, t, fs, fc, fm1, 2, filterOrder)
_, X4f, Y4f = WMS_harmonic(sig, t, fs, fc, fm1, 4, filterOrder)

_, X2f_bg, Y2f_bg = WMS_harmonic(bg_sig, t, fs, fc, fm1, 2, filterOrder)
_, X4f_bg, Y4f_bg = WMS_harmonic(bg_sig, t, fs, fc, fm1, 4, filterOrder)

R1f = np.sqrt(X1f**2 + Y1f**2)
R1f_bg = np.sqrt(X1f_bg**2 + Y1f_bg**2)
WMS_2f1f_L1 = np.sqrt(((X2f / R1f) - (X2f_bg / R1f_bg))**2 + ((Y2f / R1f) - (Y2f_bg / R1f_bg))**2)

WMS_2f_L1 = np.sqrt((X2f - X2f_bg)**2 + (Y2f - Y2f_bg)**2)
WMS_4f_L1 = np.sqrt((X4f - X4f_bg)**2 + (Y4f - Y4f_bg)**2)
WMS_4f2f_L1 = WMS_4f_L1 / WMS_2f_L1

# Laser 2 (4f outside detector bandwidth)
_, X1f, Y1f = WMS_harmonic(sig, t, fs, fc, fm2, 1, filterOrder)
_, X1f_bg, Y1f_bg = WMS_harmonic(bg_sig, t, fs, fc, fm2, 1, filterOrder)

_, X2f, Y2f = WMS_harmonic(sig, t, fs, fc, fm2, 2, filterOrder)
_, X2f_bg, Y2f_bg = WMS_harmonic(bg_sig, t, fs, fc, fm2, 2, filterOrder)

R1f = np.sqrt(X1f**2 + Y1f**2)
R1f_bg = np.sqrt(X1f_bg**2 + Y1f_bg**2)
WMS_2f1f_L2 = np.sqrt(((X2f / R1f) - (X2f_bg / R1f_bg))**2 + ((Y2f / R1f) - (Y2f_bg / R1f_bg))**2)

WMS_2f_L2 = np.sqrt((X2f - X2f_bg)**2 + (Y2f - Y2f_bg)**2)

# --------------------------- Plots --------------------------
m = 20  # downsample
Nsubplot = 2

fig, axs = plt.subplots(Nsubplot, 1, sharex=True)
axs[0].plot(t[::m] * 1e6, WMS_2f1f_L1[::m])
axs[0].spines['top'].set_visible(False); axs[0].spines['right'].set_visible(False)
axs[0].set_ylabel('1391.7 nm')

axs[1].plot(t[::m] * 1e6, WMS_2f1f_L2[::m], 'r')
axs[1].spines['top'].set_visible(False); axs[1].spines['right'].set_visible(False)
axs[1].set_ylabel('1469.3 nm')
axs[1].set_xlabel('time [$\\mu$s]')
plt.tight_layout()

fig, axs = plt.subplots(Nsubplot, 1, sharex=True)
axs[0].plot(t[::m] * 1e6, WMS_2f_L1[::m])
axs[0].spines['top'].set_visible(False); axs[0].spines['right'].set_visible(False)
axs[0].set_ylabel('WMS-2$f$')

axs[1].plot(t[::m] * 1e6, WMS_4f_L1[::m], 'r')
axs[1].set_ylim([0, 0.005])
axs[1].spines['top'].set_visible(False); axs[1].spines['right'].set_visible(False)
axs[1].set_ylabel('WMS-4$f$')
axs[1].set_xlabel('time [$\\mu$s]')
plt.tight_layout()

plt.show()
