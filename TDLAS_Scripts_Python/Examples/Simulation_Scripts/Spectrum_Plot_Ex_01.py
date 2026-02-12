
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
from Utilities.Spectrum_Plot_tools import SpectrumPlot

dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"
species = 'H2O'  # molecule
T = np.array([1200.0, 1800.0, 2200.0])  # Temperature [K]
P = 1.0  # Pressure [atm]
Xrad = 0.1  # Radiating species mole fraction
L = 6.0  # Path length [cm]
dwn = 0.001  # wavenumber spacing [cm^-1]

# wnL = 7185.0  # [cm^-1] - min wavenumber
wnH = 7186.0  # [cm^-1] - max wavenumber

wnL = 6805.5   # [cm^-1] - min wavenumber
wnH = 6806.5   # [cm^-1] - max wavenumber

wn = np.arange(wnL, wnH + dwn / 2.0, dwn)

# This function is meant as a tool for quickly plotting the absorbance
# spectrum. For more versatility see the next example.
SpectrumPlot(T, P, Xrad, L, wn, species, dat_path)
