# Error_code_HITRAN.py

import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # project root
os.chdir(Path(__file__).resolve().parent)                     # use script folder for data files

import numpy as np
import pandas as pd
from scipy.io import loadmat
from Utilities.Spectrum_Plot_tools import load_hitran_struct

"""
% Code   Relative Uncertainty range
% 0      Not Reported
% 1      Default or constant
% 2      Average or estimate
% 3      >= 20%
% 4      >= 10% and < 20%
% 5      >= 5% and < 10%
% 6      >= 2% and < 5%
% 7      >= 1% and < 2%
% 8      < 1%
"""

# ----------------------- Inputs -----------------------
wnL = 7185.59
wnH = 7185.60
species = 'H2O'

# Update this to your HITRAN data root
dat_path = "/Users/joseguerrero/Downloads/TDLAS_Scripts_Python/Hitran Data/"

if species == 'H2O':
        qfilename = dat_path + species + '/q_h2o_2020.txt'
        matfile   = dat_path + species + '/h2o_hitran2020.mat'
        MW_rad    = 18.01528
elif species == 'CO':
        qfilename = dat_path + species + '/q_co_2020.txt'
        matfile   = dat_path + species + '/co_hitran2020.mat'
        MW_rad    = 27.994915
elif species == 'CO2':
        qfilename = dat_path + species + '/q_co2_2020.txt'
        matfile   = dat_path + species + '/co2_hitran2020.mat'
        MW_rad    = 43.989830
else:
    raise ValueError(f"Unsupported species: {species}")

# ---------------- Load HITRAN struct ----------------
m = loadmat(str(matfile), squeeze_me=True, struct_as_record=False)
htd = m["htd"]

# Scalars/vectors
vo_all = np.ravel(np.array(htd.transitionWavenumber, dtype=float))  # wavenumber range
ST0_all = np.ravel(np.array(htd.lineIntensity, dtype=float))        # line intensity
E_all = np.ravel(np.array(htd.lowerStateEnergy, dtype=float))       # lower state energy

# Error codes (MATLAB char array). We’ll normalize to a list of 6-char strings, one per line
errs = htd.errorCodes
errs_arr = np.array(errs)

if errs_arr.ndim == 2:
    # join characters across columns for each row
    err_rows = [''.join(row.astype(str).tolist()) for row in errs_arr]
elif errs_arr.ndim == 1:
    err_rows = [str(s) for s in errs_arr.tolist()]
else:
    raise ValueError("Unexpected shape for errorCodes.")

# ---------------- Window select ----------------
jj = (vo_all >= wnL) & (vo_all <= wnH)
vo_vac = vo_all[jj]
ST0 = ST0_all[jj]
E = E_all[jj]

# match error rows to selection by index
idx = np.nonzero(jj)[0]
err_rows_sel = [err_rows[i] for i in idx]

# ---------------- Build table ----------------
# MATLAB used digits 2..5 (1-based); here take positions 1..4 (0-based)
tab_rows = []
for v, s0, e, code in zip(vo_vac, ST0, E, err_rows_sel):
    # guard length; pad if needed
    code = (code + "000000")[:6]
    err_block = [int(code[1]), int(code[2]), int(code[3]), int(code[4])]
    tab_rows.append([v, s0, e, *err_block])

df = pd.DataFrame(
    tab_rows,
    columns=["vo", "S0", "E", "err-S(To)", "err-g-air", "err-g-self", "err-n-air"]
)

# ---------------- Display ----------------
pd.set_option("display.precision", 4)
print("\n\n")
print(df.to_string(index=False))
print("\n\n")