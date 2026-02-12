
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt
from Utilities.Spectrum_Plot_tools import HITRAN_Spectrum_Plot, getHITRAN

# This example plots the cluster of CO lines near 5 microns commonly used in thermometry applications in the mid-IR

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/HiTemp Data/"
species = 'CO_HiTemp'  # absorbing molecule
T = 2000.0       # K
Xrad = 0.05       # radiating species mole fraction
L = 2.0          # cm path length

dwn = 0.001      # cm^-1 spacing
wnL = 2008.0     # cm^-1 low
wnH = 2009.0     # cm^-1 high

threshold = 10**-4
wn_units = 1
S_units = 1

P = 2.0          # atm
X1, abs1 = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

P = 5.0          # atm
X2, abs2 = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

Xrad = 1
vo,_,XS,_,_,_,_,_ = getHITRAN(T, Xrad, wnL, wnH, threshold, species, wn_units, S_units, dat_path)

plt.figure(figsize=(7.5, 6))
plt.plot(X1, abs1) 
plt.plot(X2, abs2) 
plt.stem(vo, XS, linefmt='k-', markerfmt=' ', basefmt=' ')
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.xlim([wnL, wnH])
plt.tight_layout()
plt.show()

