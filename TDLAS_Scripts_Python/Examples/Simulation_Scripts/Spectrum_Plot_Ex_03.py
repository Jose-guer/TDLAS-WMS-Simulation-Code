
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt
from Utilities.Spectrum_Plot_tools import HITRAN_Spectrum_Plot_Alt, HITRAN_Spectrum_Plot

# This example simulates the absorption spectrum in two ways. The first
# uses the specified T, P, Xi and the second using T, Pi, and wLi.

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"
species = 'H2O'   # molecule
T = 1800.0        # Temperature [K]
P = 1.0           # Pressure [atm]
Xrad = 0.1        # radiating species mole fraction
Pi = P * Xrad
L = 6.0           # Path length [cm]
dwn = 0.001       # wavenumber spacing [cm^-1]

############################################################
# wnL = 7185.0
# wnH = 7186.0
# gamma_self = 0.198  # 7185.596
# n_self = 0.53
# gamma_air = 0.0403
# n_air = 0.587
############################################################

wnL = 6805.5
wnH = 6806.5

gamma_self = 0.205  # 6806.03
n_self = 0.86
gamma_air = 0.0104
n_air = -0.164

############################################################

wn = np.arange(wnL, wnH + dwn / 2.0, dwn)

# estimate collisional width (HWHM) for comparison
wLi = P * (Xrad * gamma_self * (296.0 / T)**n_self + (1.0 - Xrad) * gamma_air * (296.0 / T)**n_air)

X1, absor1 = HITRAN_Spectrum_Plot_Alt(T, Pi, wLi, L, dwn, wnL, wnH, species, dat_path)

X2, absor2 = HITRAN_Spectrum_Plot(T, P, Xrad, L, dwn, wnL, wnH, species, dat_path)

plt.figure(figsize=(7.32, 5.95))
plt.plot(X1, absor1)
plt.plot(X2, absor2)
plt.legend([r'$\alpha(\nu,T,P_i,\Delta\nu_c)$', r'$\alpha(\nu,T,P,X_i)$'], frameon=True)
plt.xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
plt.ylabel(r'$\alpha(\nu)$')
plt.xlim([wnL, wnH])
plt.tight_layout()
plt.show()
