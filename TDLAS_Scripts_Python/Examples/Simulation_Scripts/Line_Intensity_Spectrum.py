# Line_Intensity_Spectrum.py

"""
This script simulates the line intensity spectrum at a specified
temperature and species mole fraction

written by Jose Guerrero, University of Michigan - Aerospace Department
joseguer@umich.edu

"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # adds project root

import numpy as np
import matplotlib.pyplot as plt
from Utilities.Spectrum_Plot_tools import getHITRAN

# ************ SIMULATION INPUTS ***************

T = 1800.0  # K
X = 1.0
threshold = 1e-4
wnL = 1.0e4 / 6.0  # lower wavenumber
wnH = 1.0e4 / 1.0  # upper wavenumber
wn_units = 0       # 0 for μm, 1 for cm⁻¹
S_units = 1
Xrad = 1

# ************ H2O *****************************

species = 'H2O'
dat_path  = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"

plt.figure(figsize=(11.94, 4.98))  # same physical dimensions as MATLAB [inches]
vo, _,XS,_,_,_,_,_ = getHITRAN(T, Xrad, wnL, wnH, threshold, species, wn_units, S_units, dat_path)

plt.stem(vo, XS, linefmt='b-', markerfmt=' ', basefmt=' ')
plt.yscale('log')
plt.legend([species], loc='upper left', frameon=True)
plt.ylabel(r'$X_i \cdot S(T)~\mathrm{[cm^{-2}/atm]}$')

if wn_units == 1:
    plt.xlim([wnL, wnH])
    plt.xlabel(r'wavenumber $\mathrm{[cm^{-1}]}$')
else:
    plt.xlim([1, 6])
    plt.xlabel(r'wavelength $(\mu\mathrm{m})$')

plt.tight_layout()

# ************ CO2 AND 13CO2 *******************

plt.figure(figsize=(11.94, 4.98))

species = 'CO2'

vo, _,XS,_,_,_,_,_ = getHITRAN(T, Xrad, wnL, wnH, threshold, species, wn_units, S_units, dat_path)
plt.stem(vo, XS, linefmt='b-', markerfmt=' ', basefmt=' ')

species = '13CO2'
vo, _,XS,_,_,_,_,_ = getHITRAN(T, Xrad, wnL, wnH, threshold, species, wn_units, S_units, dat_path)
plt.stem(vo, XS, linefmt='r-', markerfmt=' ', basefmt=' ')

plt.yscale('log')
plt.legend(['CO2', '13CO2'], loc='upper left', frameon=True)
plt.ylabel(r'$X_i \cdot S(T)~\mathrm{[cm^{-2}/atm]}$')

if wn_units == 1:
    plt.xlim([wnL, wnH])
    plt.xlabel(r'wavenumber $\mathrm{[cm^{-1}]}$')
else:
    plt.xlim([1, 6])
    plt.xlabel(r'wavelength $(\mu\mathrm{m})$')

plt.tight_layout()
plt.show()
