# Discrete_Intensity_Spectrum.py
"""
This script simulates the discrete intensity spectrum of modulated light 
following the theory outlined in :

Mathews, Garrett, and Christopher Goldenstein. 
"Near-GHz scanned-wavelength-modulation spectroscopy for MHz thermometry 
and H2O measurements in aluminized fireballs of energetic materials." 
Applied Physics B 126.11 (2020): 189.

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu

"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv  # Bessel J

# 1469 nm Laser
lambda_nm = 1469.28
fscan = 500e3
fmod = 45.5e6

# Laser characterization
a = 0.0795
aHz = 2.385e9
i0 = 0.175
i2 = 8.394e-3
psi = 5.88
psi2 = 4.22

c = 3e10           # cm/s
vc = 1e7 / lambda_nm  # cm^-1

M = i0 / 2.0
beta = aHz / fmod

K = 500
k = np.arange(-K, K + 1)

v = vc + k * (fmod / c)

Jk  = jv(k, beta)
Jkm = jv(k - 1, beta)
Jkp = jv(k + 1, beta)

E = Jk + (M / (2j)) * np.exp(1j * psi) * Jkm - (M / (2j)) * np.exp(-1j * psi) * Jkp

F = c * v / 1e9
F = F - F.mean()

intensity = np.abs(E) ** 2
Im = intensity.max()

plt.figure()
(markerline, stemlines, baseline) = plt.stem(F, intensity / Im, linefmt='k-', markerfmt=' ', basefmt=' ')
plt.xlabel('Relative frequency [GHz]')
plt.ylabel('Discrete spectrum')
plt.xlim([-3, 3])
plt.tight_layout()
plt.show()
