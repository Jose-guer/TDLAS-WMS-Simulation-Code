# TDLAS-WMS-Simulation-Code
MATALB and Python code for the simulation of molecular absorption spectra using the HITRAN and HITEMP databases is made available here as supplementary material to our review paper cited below [[1](https://doi.org/10.48550/arXiv.2512.18201)]. Code for the simulation of WMS harmonics using the calibration-free WMS model, and laser characterization examples with data are also included. 

## System requirements
The MATALB implementation of the code requires the signal processing toolbox. 

## Description of Code
The code is organized into the following folders:
  - Examples. - This folder contains the main simulation scripts.
  - Functions. - This folder contains the custom functions used in the examples and are further divided into scripts used to simulated the molecular absorption spectra, and those used to simulate WMS            harmonics, including a digital lock-in amplifier.
  - Hitran Data. - This folder contains hitran data for various molecules of interest to combustion diagnostics (H2O, CO, CO2, NO, etc.) and is formated to be used with the provided simulation scripts.
  - HiTemp Data. - This folder contains mid-IR data for CO, CO2, and H2O.

## Please cite this article when using code provided here in your work
Guerrero, Jose, and Mirko Gamba. "A Review of Theory and Practical Considerations of Tunable Diode Laser Absorption Spectroscopy Diagnostics." arXiv:2512.18201 (2025) 
https://doi.org/10.48550/arXiv.2512.18201

## Reference
[1] Guerrero, Jose, and Mirko Gamba. "A Review of Theory and Practical Considerations of Tunable Diode Laser Absorption Spectroscopy Diagnostics." arXiv:2512.18201 (2025) 
https://doi.org/10.48550/arXiv.2512.18201
