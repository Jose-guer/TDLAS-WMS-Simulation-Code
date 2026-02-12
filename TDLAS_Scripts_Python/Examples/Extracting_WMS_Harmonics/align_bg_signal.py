
"""
This script implements an algorithm based on translating the background
intensity signal and computes the root-mean-squared-deviation. The 
signals are aligned at the minimum RMSD. The script assumes that both
modulation frequencies are divisible by fscan such every scan is periodic 
and therefore only one period of the scanning frequency needs to be 
searched to align the signals. 

This script loads in the sample data "raw_data.mat" and saves the
"aligned_data.mat" file which is used by "plot_WMS_Harmonics_LPF.m" to
extract the background subtracted WMS harmonics.

The script crops the background singal here to 26 scanning periods but
this number can be chosen to be any length greater than 2 scanning
periods upto the total length of the absorption signal. 

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
from scipy.io import loadmat, savemat

# Load sample data (expects raw_data.mat in the working directory)
d = loadmat('raw_data.mat', squeeze_me=True, struct_as_record=False)
raw = d['raw']
t = np.asarray(raw.t).ravel()
sig = np.asarray(raw.sig).ravel()
bg_sig = np.asarray(raw.bgsig).ravel()
fscan = float(raw.fscan)
fs = float(raw.fs)

Nsamples = int(round(fs / fscan))          # samples per scan period
bg_sig = bg_sig[: 100 * Nsamples]          # crop background signal

# Copies for saving / scaling
sig_copy = sig.copy()
bg_sig_copy = bg_sig.copy()

# Normalize by the peak over one scan (first period)
sig = sig / np.max(sig[:Nsamples])
bg_sig = bg_sig / np.max(bg_sig[:Nsamples])

# Plot initial normalized, unaligned signals (first period)
fig1, ax1 = plt.subplots()
ax1.plot(sig[:Nsamples])
ax1.plot(bg_sig[:Nsamples])
ax1.set_title('Normalized')
ax1.set_xlabel('Sample')
ax1.set_ylabel('Amplitude')
ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
ax1.legend(['sig', 'bg'])

# Find minimum residual alignment for each block
Nbg = len(bg_sig) - Nsamples         # effective background length per block
Nsig = len(sig)

P = Nsig // Nbg                      # number of full blocks
R = Nsig - P * Nbg                   # remainder length
N = P if R == 0 else P + 1

sig_bg = np.zeros(Nsig)
RMSD = None
Nshift = None

if R != 0:
    print(f'R = {R}')

for i in range(1, N + 1):
    if i < N:
        iL = Nbg * (i - 1)
        iR = Nbg * i
        this_Nbg = Nbg
    else:
        iL = Nbg * (i - 1)
        iR = Nsig
        this_Nbg = R if R != 0 else Nbg
        if R != 0:
            print('Length of last section is R')

    sig_seg = sig[iL:iR]
    RMSD_vals = np.zeros(Nsamples)
    # j = 1..Nsamples in MATLAB -> 0..Nsamples-1 in Python
    for j in range(Nsamples):
        bg_seg = bg_sig[j : j + this_Nbg]
        SSR = np.sum((sig_seg - bg_seg) ** 2)
        RMSD_vals[j] = np.sqrt(SSR / this_Nbg)

    # Track last block's RMSD/Nshift for plotting
    RMSD = RMSD_vals
    Nshift = int(np.argmin(RMSD_vals))  # 0-based index

    # Scale intensity of the aligned background to match the signal segment
    bg_seg_copy = bg_sig_copy[Nshift : Nshift + this_Nbg]
    ratio = np.max(sig_copy[iL:iR]) / np.max(bg_seg_copy)
    sig_bg[iL:iR] = bg_seg_copy * ratio

# Plot RMSD vs shift for the last-aligned block
fig2, ax2 = plt.subplots()
ax2.plot(RMSD)
ax2.plot(Nshift, RMSD[Nshift], 'ro', linewidth=2)
ax2.set_xlabel('Sample Offset')
ax2.set_ylabel('RMSD')
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)

# Plot a few scan periods of aligned signals (normalized)
Nsig_plot = int(round(10 * fs / fscan))
fig3, ax3 = plt.subplots()
ax3.plot(t[:Nsig_plot], sig_copy[:Nsig_plot] / np.max(sig_copy[:Nsig_plot]))
ax3.plot(t[:Nsig_plot], sig_bg[:Nsig_plot] / np.max(sig_bg[:Nsig_plot]))
ax3.legend(['sig', 'bg sig'])
ax3.set_xlabel('time [s]')
ax3.spines['top'].set_visible(False); ax3.spines['right'].set_visible(False)

# Save aligned data
dat = np.column_stack((t, sig_copy, sig_bg))
#savemat('aligned_data.mat', {'dat': dat})

plt.tight_layout()
plt.show()
