import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods

# %matplotlib notebook  # Uncomment if running interactively

# === Path for only the 02_100K case ===
paths_100K_only = {
    2: "02_100K/step2/step3/res-sd-9-occ/"
}

system_labels = {
    2: "(TiO$_{2}$)$_{2}$"
}

def analyze_energy_gap_0_1(path_dict, temperature):
    plt.figure()

    cmap = cm.get_cmap("coolwarm")
    norm = mcolors.Normalize(vmin=min(path_dict.keys()), vmax=max(path_dict.keys()))

    for np_units, folder in path_dict.items():
        energy_gaps = []

        params = {
            "path_to_energy_files": folder,
            "dt": 1.0,
            "prefix": "Hvib_ci_",
            "suffix": "_re",
            "istep": 1000,
            "fstep": 3990
        }

        try:
            md_time, energies = CP2K_methods.extract_energies_sparse(params)
        except Exception as e:
            print(f"[Error] Couldn't extract energies for {folder}: {e}")
            continue

        energies = energies * units.au2ev  # Convert to eV

        if energies.shape[1] < 2:
            print(f"[Skip] Not enough states in {folder}")
            continue

        # Only gap between states 0 and 1
        gap = np.abs(energies[:, 1] - energies[:, 0])
        for val in gap:
            x_mb = MATRIX(1, 1)
            x_mb.set(0, 0, val)
            energy_gaps.append(x_mb)

        if energy_gaps:
            bin_supp, dens, cum = data_stat.cmat_distrib(energy_gaps, 0, 0, 0, 0, 50, 0.01)
            label = system_labels.get(np_units, f"System {np_units}")
            color = cmap(norm(np_units))
            plt.plot(bin_supp, dens, label=label, color=color)

    plt.xlabel('Energy gap (S$_1$ − S$_0$) (eV)', fontsize=20)
    plt.ylabel('PD (1/eV)', fontsize=20)
    plt.yscale('log')
    plt.xlim(2, 5)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'energy_gap_0_1_{temperature}K.jpg', dpi=600)

# === Run only for 100K, 02 case ===
analyze_energy_gap_0_1(paths_100K_only, 100)