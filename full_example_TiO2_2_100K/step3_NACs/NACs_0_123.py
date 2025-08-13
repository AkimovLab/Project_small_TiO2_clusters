
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat

# %matplotlib notebook  # Uncomment if running interactively

# === Path for only the 02_100K case ===
paths_100K_only = {
    2: "results/"
}

system_labels = {
    2: "(TiO$_{2}$)$_{2}$"
}

def analyze_nac_01(path_dict, temperature):
    plt.figure()
    
    cmap = cm.get_cmap("coolwarm")
    norm = mcolors.Normalize(vmin=min(path_dict.keys()), vmax=max(path_dict.keys()))

    for np_units, folder in path_dict.items():
        nac_values = []

        nac_files = glob.glob(f"{folder}/Hvib_ci*im*")
        if not nac_files:
            print(f"[Warning] No NAC file found in {folder}")
            continue

        for nac_file in nac_files:
            hvib = sp.load_npz(nac_file)
            hvib_dense = hvib.todense().real

            if hvib.shape[0] < 4:
                print(f"[Skip] File {nac_file} has fewer than 4 states")
                continue

            for j in [1, 2, 3]:
                if j < hvib.shape[0]:
                    nac_0j = np.abs(hvib_dense[0, j]) * 1000.0 * units.au2ev  # in meV
                    x_mb = MATRIX(1, 1)
                    x_mb.set(0, 0, nac_0j)
                    nac_values.append(x_mb)

        if nac_values:
            bin_supp, dens, cum = data_stat.cmat_distrib(nac_values, 0, 0, 0, 0, 50, 0.1)
            label = system_labels.get(np_units, f"System {np_units}")
            color = cmap(norm(np_units))
            plt.plot(bin_supp, dens, label=label, color=color)

    plt.xlabel('|NAC (S$_0$, S$_1$–S$_3$)| (meV)', fontsize=20)
    plt.ylabel('PD (1/meV)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0, 2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'nac_0j_dist_{temperature}K.jpg', dpi=600, bbox_inches='tight')

# === Run only for 100K, 02 case ===
analyze_nac_01(paths_100K_only, 100)