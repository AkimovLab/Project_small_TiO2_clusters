import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods
import numpy as np

# %matplotlib notebook  # Uncomment if running interactively

# === Path for only the 02_100K case ===
paths_100K_only = {
    2: "results/"
}

system_labels = {
    2: "(TiO$_{2}$)$_{2}$"
}

def analyze_energy_gap_fluctuations(path_dict, temperature):
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

        energies = energies * units.au2ev  # convert to eV

        # Adjacent gaps
        for i in range(1, energies.shape[1] - 1):
            gap_fwd = energies[:, i+1] - energies[:, i]
            gap_bwd = energies[:, i] - energies[:, i-1]
            for val in np.concatenate((gap_fwd, gap_bwd)):
                x_mb = MATRIX(1, 1)
                x_mb.set(0, 0, abs(val))
                energy_gaps.append(x_mb)

        # Forward gaps again (matches original structure)
        for i in range(1, energies.shape[1] - 1):
            gap = energies[:, i+1] - energies[:, i]
            for val in gap:
                x_mb = MATRIX(1, 1)
                x_mb.set(0, 0, abs(val))
                energy_gaps.append(x_mb)

        if energy_gaps:
            bin_supp, dens, cum = data_stat.cmat_distrib(energy_gaps, 0, 0, 0, 0, 50, 0.01)
            label = system_labels.get(np_units, f"System {np_units}")
            color = cmap(norm(np_units))
            plt.plot(bin_supp, dens, label=label, color=color)

    plt.xlabel('Energy gap (S$_{\\mathrm{i}}$, S$_{\\mathrm{i}\\pm1}$) (eV)', fontsize=20)
    plt.ylabel('PD (1/eV)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'energy_gap_fluct_{temperature}K.jpg', dpi=600)

# === Run only for 100K, 02 case ===
analyze_energy_gap_fluctuations(paths_100K_only, 100)