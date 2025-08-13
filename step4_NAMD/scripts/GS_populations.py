import os, glob, time, h5py, warnings
import multiprocessing as mp
import matplotlib.pyplot as plt   # plots
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit

from liblibra_core import *
import util.libutil as comn

import matplotlib.ticker as ticker

import libra_py
from libra_py import units, data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import libra_py.data_visualize

warnings.filterwarnings('ignore')
# Define functions
def exp_functfssh(t, _x):
    return 1 - np.exp(-np.power(t / _x, 2.0))

def exp_functdish(t, _x):
    return 1 - np.exp(-np.power(t / _x, 1.0))

# Define separate ICONDS lists for MSDM and FSSH
ICONDS = list(range(1, 3000, 100))  # FSSH: From 1 to 3000, step 100

# Create a single figure
fig, ax = plt.subplots(figsize=(12, 11))  # Single subplot for combined plot

methods = ['FSSH',  'DISH']
colors = ['red', 'green']  # Unique colors for each method

# Store timescale strings for both methods
timescale_labels = []

# Loop over methods and apply different exponential functions based on the method
for idx, method in enumerate(methods):
    taus = []

    # Select the appropriate exponential function for the method
    if method == 'DISH':
        exp_func = exp_functdish
    elif method == 'FSSH':
        exp_func = exp_functfssh
    else:
        raise ValueError(f"Unknown method: {method}")

    for icond in ICONDS:
        try:
            F = h5py.File(f'{method}_icond_{icond}/mem_data.hdf')
            sh_pop = np.array(F['../key_outputs/sh_pop_adi/data'][:, 0])
            md_time = np.array(F['../key_outputs/time/data'][:]) * units.au2fs
#            ax.plot(md_time, sh_pop, alpha=0.1, color=colors[idx])  # Light background curves
            F.close()

            # Fit data using the selected exponential function
            popt, pcov = curve_fit(exp_func, md_time, sh_pop, bounds=([0.0], [np.inf]))
            _tau = popt

            # Compute R-squared
            residuals = sh_pop - exp_func(md_time, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sh_pop - np.mean(sh_pop))**2)
            r_squared = 1.0 - (ss_res / ss_tot)
            print(f"Method: {method}, IC {icond}, R2 = {r_squared}")

            if r_squared > 0.0:
                taus.append(_tau)
        except Exception as e:
            print(f"Error processing {method}, IC {icond}: {e}")
            continue

    taus = np.array(taus)
    ave_tau = np.average(taus)
    s = np.std(taus)
    Z = 1.96
    N = taus.shape[0]
    error_bar = Z * s / np.sqrt(N)

    # Timescale in ps with 3 decimal places
    ave_tau_ps = ave_tau / 1000
    error_bar_ps = error_bar / 1000
    print(f'Timescales for {method}: {ave_tau_ps:.3f} ± {error_bar_ps:.3f}')

    # Append timescale label for the upper-left region
    timescale_labels.append(f"{method}: {ave_tau_ps:.3f} ± {error_bar_ps:.3f} ps")

    # Plot fitted curve and error bars
    ax.plot(md_time, exp_func(md_time, ave_tau-error_bar), ls='--', linewidth=2, color=colors[idx])
    ax.plot(md_time, exp_func(md_time, ave_tau), ls='-', linewidth=6, color=colors[idx], label=f"{method}")
    ax.plot(md_time, exp_func(md_time, ave_tau+error_bar), ls='--', linewidth=2, color=colors[idx])
#    ax.set_xlabel("Time (fs)", fontsize=30)
ax.set_ylabel("Ground State Population", fontsize=45)

# Add a title (optional)
ax.set_title("(TiO$_{2})_{2}$ 100 K", fontsize=50)

# Add a legend
# Customize tick parameters for both axes
ax.tick_params(axis='both', which='major', labelsize=40)  # Larger font for major ticks
# Adjust timescale box for center-left positioning
timescale_text = "\n".join(timescale_labels)
props = dict(boxstyle='round', facecolor='white', alpha=0.8)
#ax.text(0.05, 0.5, timescale_text, transform=ax.transAxes, fontsize=20,
#        verticalalignment='center', horizontalalignment='left', bbox=props)

# Show grid for better readability
plt.savefig('TiO2_2_100K_GS_pop_evolution.pdf', format='pdf')

# Show the plot
plt.tight_layout()
plt.show()



