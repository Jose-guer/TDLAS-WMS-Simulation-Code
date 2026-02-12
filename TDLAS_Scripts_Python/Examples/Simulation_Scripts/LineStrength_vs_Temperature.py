# line_strength_ratio.py

"""
This script plots the normalized linestrength function S(T)/S(T0) over a
range of temperatures for specified lower-state energies (Eg)

written by Jose Guerrero, University of Michigan - Aerospace Department 
joseguer@umich.edu
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt

# Simulation details
dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/H2O/"
qfilename = "q_h2o_2020.txt"
Eg = np.array([1045.0, 3291.0])
Trange = np.arange(300.0, 3000.0 + 10.0, 10.0)
ymax = 0         # 1 to scale by peak amplitude
logplot = 1      # 1 for semilogy
S_P_units = 1    # 0: cm^-1/(molecule-cm^-2), 1: cm^-2/atm

# Constants
h = 6.6260680e-34
k = 1.3806503e-23
c = 2.99792458e10
Tref = 296.0

# Partition function fit
dat = np.loadtxt(Path(dat_path) / qfilename)
pn = np.polyfit(dat[:, 0], dat[:, 1], 6)
def Q(T): return np.polyval(pn, T)

# Line strength ratio
LineStrengthRatio = np.zeros((Eg.size, Trange.size))
for i in range(Eg.size):
    if S_P_units == 1:
        LineStrengthRatio[i, :] = (Q(Tref)/Q(Trange)) * (Tref/Trange) * np.exp(-h*c*Eg[i]/k * (1.0/Trange - 1.0/Tref))
    else:
        LineStrengthRatio[i, :] = (Q(Tref)/Q(Trange)) * np.exp(-h*c*Eg[i]/k * (1.0/Trange - 1.0/Tref))

legStr = [rf"$E''$ = {int(E)} $\mathrm{{cm^{{-1}}}}$" for E in Eg]

# Plot
fig, ax = plt.subplots(figsize=(6.9, 5.67))

if ymax == 1:
    Smax = LineStrengthRatio.max(axis=1, keepdims=True)
    Y = LineStrengthRatio / Smax
else:
    Y = LineStrengthRatio

if logplot == 1:
    for i in range(Eg.size): ax.semilogy(Trange, Y[i, :])
else:
    for i in range(Eg.size): ax.plot(Trange, Y[i, :])

ax.axvline(1000.0, color='k', linestyle='--', linewidth=1.0)
ax.axvline(2500.0, color='k', linestyle='--', linewidth=1.0)

ax.set_xlabel('T [K]')
ax.legend(legStr, frameon=True, loc='lower left')
ax.grid(True, which='both')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_xlim(Trange.min(), Trange.max())

if S_P_units == 1:
    ax.set_ylabel(r'$S(T)/S(T_0)$')
else:
    ax.set_ylabel(r'$S^*(T)/S^*(T_0)$')

plt.tight_layout()
plt.show()
