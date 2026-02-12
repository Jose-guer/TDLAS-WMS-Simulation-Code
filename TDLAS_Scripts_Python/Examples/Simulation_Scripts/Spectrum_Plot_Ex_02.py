
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt
from Utilities.Spectrum_Plot_tools import HITRAN_Spectrum_Plot

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"
species = 'H2O'  # absorbing molecule
T = 1200.0       # K
P = 1.0          # atm
Xrad = 0.1       # radiating species mole fraction
L = 6.0          # cm path length

dwn = 0.001      # cm^-1 spacing
wnL = 7185.2     # cm^-1 low
wnH = 7186.0     # cm^-1 high

X, absorption = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

#np.savetxt("spectrum.csv", np.column_stack((X, absorption)),
#            delimiter=",", header="wavenumber,absorption", comments="")

plt.figure(figsize=(7.5, 6))
plt.plot(X, absorption) 
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.xlim([wnL, wnH])
plt.tight_layout()
plt.show()

