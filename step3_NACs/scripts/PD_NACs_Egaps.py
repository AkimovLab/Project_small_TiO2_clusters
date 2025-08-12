#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods 

get_ipython().run_line_magic('matplotlib', 'notebook')

# === Folder paths for each system ===
paths_100K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/100K/newAS_8_11/",
    2: "02_100K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    5: "0500_100K/isomer/step3/res-sd-30-occ/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    6: "0600_100K/step3/res-sd-36-occ/",
    7: "07_100/0700/100K/notaligned/step3/res-sd-42-occ/",
    8: "../100K/100K/step3/res-sd-48-occ/"
}

paths_300K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/300K/newAS_8_11/",
    2: "02_300K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_24_25/",
    5: "0500_300K/isomer2/step3/AS_30_39/",
    6: "0600_300K/step3/res-sd-33-occ/",
    7: "0700/300K/step3/res-sd-42-occ/",
    8: "../08300/300K/300K/step3/res-sd-48-occ/"
}

system_labels = {
    1: "(TiO$_{2}$)$_{1}$",
    2: "(TiO$_{2}$)$_{2}$",
    3: "(TiO$_{2}$)$_{3}$",
    4: "(TiO$_{2}$)$_{4}$",
    5: "(TiO$_{2}$)$_{5}$",
    6: "(TiO$_{2}$)$_{6}$",
    7: "(TiO$_{2}$)$_{7}$",
    8: "(TiO$_{2}$)$_{8}$",
}

def analyze_nac_01(path_dict, temperature):
    plt.figure()
    
    # Setup colormap from blue to red for indices 1 through 8
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
            color = cmap(norm(np_units))  # Get color for this system
            plt.plot(bin_supp, dens, label=label, color=color)

    plt.xlabel('|NAC (S$_0$, S$_1$-S$_3$)| (meV)', fontsize=20)
    plt.ylabel('PD (1/meV)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0, 2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    # 4 rows × 2 columns legend layout
    plt.legend(fontsize=15, ncol=2, loc='lower left')
    
    plt.tight_layout()
    plt.savefig(f'nac_0j_dist_{temperature}K.jpg', dpi=600, bbox_inches='tight')
    
# === Run for both temperatures ===
analyze_nac_01(paths_100K, 100)
analyze_nac_01(paths_300K, 300)


# In[8]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods 

get_ipython().run_line_magic('matplotlib', 'notebook')

# === Folder paths for each system ===
paths_100K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/100K/newAS_8_11/",
    2: "02_100K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    5: "0500_100K/isomer/step3/res-sd-30-occ/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    6: "0600_100K/step3/res-sd-36-occ/",
    7: "07_100/0700/100K/notaligned/step3/res-sd-42-occ/",
    8: "../100K/100K/step3/res-sd-48-occ/"
}

paths_300K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/300K/newAS_8_11/",
    2: "02_300K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_24_25/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_24_25/",
    5: "0500_300K/isomer2/step3/AS_30_39/",
    6: "0600_300K/step3/res-sd-33-occ/",
    7: "0700/300K/step3/res-sd-42-occ/",
    8: "../08300/300K/300K/step3/res-sd-48-occ/"
}

system_labels = {
    1: "(TiO$_{2}$)$_{1}$",
    2: "(TiO$_{2}$)$_{2}$",
    3: "(TiO$_{2}$)$_{3}$",
    4: "(TiO$_{2}$)$_{4}$",
    5: "(TiO$_{2}$)$_{5}$",
    6: "(TiO$_{2}$)$_{6}$",
    7: "(TiO$_{2}$)$_{7}$",
    8: "(TiO$_{2}$)$_{8}$",
}

def analyze_nac_adjacent(path_dict, temperature):
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
            nstates = hvib.shape[0]

            if nstates < 3:
                print(f"[Skip] File {nac_file} has fewer than 3 states")
                continue

            # Go over adjacent pairs: 1-2, 2-3, ..., n-1 - n
            for i in range(1, nstates):  # start at 1 to skip state 0
    # i-1 coupling
                if i - 1 >= 1:
                   nac_bwd = np.abs(hvib_dense[i, i - 1]) * 1000.0 * units.au2ev
                   x_mb = MATRIX(1, 1)
                   x_mb.set(0, 0, nac_bwd)
                   nac_values.append(x_mb)

            # i+1 coupling
                if i + 1 < nstates:
                   nac_fwd = np.abs(hvib_dense[i, i + 1]) * 1000.0 * units.au2ev
                   x_mb = MATRIX(1, 1)
                   x_mb.set(0, 0, nac_fwd)
                   nac_values.append(x_mb)
            for i in range(1, nstates - 1):
                j = i + 1
                nac_ij = np.abs(hvib_dense[i, j]) * 1000.0 * units.au2ev  # meV
                x_mb = MATRIX(1, 1)
                x_mb.set(0, 0, nac_ij)
                nac_values.append(x_mb)

        if nac_values:
            bin_supp, dens, cum = data_stat.cmat_distrib(nac_values, 0, 0, 0, 0, 50, 0.1)
            label = system_labels.get(np_units, f"System {np_units}")
            color = cmap(norm(np_units))
            plt.plot(bin_supp, dens, label=label, color=color)

    plt.xlabel('|NAC (S$_{\\mathrm{i}}$, S$_{\\mathrm{i}\\pm1}$)| (meV)', fontsize=20)
    plt.ylabel('PD (1/meV)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(20, 50)
    plt.ylim(0.0001, 0.05)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
  #  plt.legend(fontsize=14, ncol=2, loc='lower left')
    plt.tight_layout()
    plt.savefig(f'nac_adjacent_dist_{temperature}K_ZOOM.jpg', dpi=600, bbox_inches='tight')

# === Run for both temperatures ===
analyze_nac_adjacent(paths_100K, 100)
analyze_nac_adjacent(paths_300K, 300)


# In[4]:


import matplotlib.cm as cm
import matplotlib.colors as mcolors

get_ipython().run_line_magic('matplotlib', 'notebook')

# === Folder paths for each system ===
paths_100K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/100K/newAS_8_11/",
    2: "02_100K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_19_19/",
    5: "0500_100K/isomer/step3/res-sd-30-occ/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    6: "0600_100K/step3/res-sd-36-occ/",
    7: "07_100/0700/100K/notaligned/step3/res-sd-42-occ/",
    8: "../100K/100K/step3/res-sd-48-occ/"
}

paths_300K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/300K/newAS_8_11/",
    2: "02_300K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_24_25/",
    5: "0500_300K/isomer2/step3/AS_30_39/",
    6: "0600_300K/step3/res-sd-33-occ/",
    7: "0700/300K/step3/res-sd-42-occ/",
    8: "../08300/300K/300K/step3/res-sd-48-occ/"
}

system_labels = {
    1: "(TiO$_{2}$)$_{1}$",
    2: "(TiO$_{2}$)$_{2}$",
    3: "(TiO$_{2}$)$_{3}$",
    4: "(TiO$_{2}$)$_{4}$",
    5: "(TiO$_{2}$)$_{5}$",
    6: "(TiO$_{2}$)$_{6}$",
    7: "(TiO$_{2}$)$_{7}$",
    8: "(TiO$_{2}$)$_{8}$",
}

def analyze_energy_gap_fluctuations(path_dict, temperature):
    plt.figure()

    # Colormap setup: coolwarm from blue (low index) to red (high index)
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

        # Only adjacent states, skipping 0→1
        for i in range(1, energies.shape[1] - 1):  # avoid state 0 and last index for i+1
        # i to i+1 gap
          gap_fwd = energies[:, i+1] - energies[:, i]
        # i to i-1 gap (i starts at 1, so i-1 ≥ 0)
          gap_bwd = energies[:, i] - energies[:, i-1]
    
          for val in np.concatenate((gap_fwd, gap_bwd)):
              x_mb = MATRIX(1, 1)
              x_mb.set(0, 0, abs(val))
              energy_gaps.append(x_mb)
        for i in range(1, energies.shape[1] - 1):  # skip state 0
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
   # plt.legend(fontsize=14, ncol=2, loc='lower left')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'energy_gap_fluct_{temperature}K.jpg', dpi=600)
# === Run for both temperatures ===
analyze_energy_gap_fluctuations(paths_100K, 100)
analyze_energy_gap_fluctuations(paths_300K, 300)


# In[5]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods 

get_ipython().run_line_magic('matplotlib', 'notebook')


# === Folder paths for each system ===
paths_100K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/100K/newAS_8_11/",
    2: "02_100K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_19_19/",
    5: "0500_100K/isomer/step3/res-sd-30-occ/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/100K/AS_24_25/",
    6: "0600_100K/step3/res-sd-36-occ/",
    7: "07_100/0700/100K/notaligned/step3/res-sd-42-occ/",
    8: "../100K/100K/step3/res-sd-48-occ/"
}

paths_300K = {
    1: "TiO2_1/H2O_0/step3_NACs_CAM_B3LYP/300K/newAS_8_11/",
    2: "02_300K/step2/step3/res-sd-9-occ/",
    3: "TiO2_3/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_19_19/",
    4: "TiO2_4/H2O_0/step3_NACs_CAM_B3LYP/300K/AS_24_25/",
    5: "0500_300K/isomer2/step3/AS_30_39/",
    6: "0600_300K/step3/res-sd-33-occ/",
    7: "0700/300K/step3/res-sd-42-occ/",
    8: "../08300/300K/300K/step3/res-sd-48-occ/"
}

system_labels = {
    1: "(TiO$_{2}$)$_{1}$",
    2: "(TiO$_{2}$)$_{2}$",
    3: "(TiO$_{2}$)$_{3}$",
    4: "(TiO$_{2}$)$_{4}$",
    5: "(TiO$_{2}$)$_{5}$",
    6: "(TiO$_{2}$)$_{6}$",
    7: "(TiO$_{2}$)$_{7}$",
    8: "(TiO$_{2}$)$_{8}$",
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
   # plt.legend(fontsize=14, ncol=2, loc='upper right')
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlim(2, 5)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'energy_gap_0_1_{temperature}K.jpg', dpi=600)
analyze_energy_gap_0_1(paths_100K, 100)
analyze_energy_gap_0_1(paths_300K, 300)


# In[ ]:




