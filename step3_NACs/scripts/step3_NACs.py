import os, glob, time, h5py, warnings
import numpy as np
import scipy.sparse as sp
from libra_py import units, data_stat, influence_spectrum, data_conv
import matplotlib.pyplot as plt
from liblibra_core import *
from libra_py.workflows.nbra import step3
import libra_py.packages.cp2k.methods as CP2K_methods
import multiprocessing as mp

import util.libutil as comn
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times

params_active_space = {
    'lowest_orbital': 12-10, 'highest_orbital': 13+10, 'num_occ_orbitals': 6, 'num_unocc_orbitals': 8,
    'path_to_npz_files': os.getcwd()+'/../../step2_Exci_CAM_B3LYP/300K/longer/res', 'logfile_directory': os.getcwd()+'/../../step2_Exci_CAM_B3LPYP/300K/longer/all_logfiles',
    'path_to_save_npz_files': os.getcwd()+'/../../step2_Exci_CAM_B3LYP/300K/longer/new_res'
}
new_lowest_orbital, new_highest_orbital = step3.limit_active_space(params_active_space)

params_mb_sd = {
          'lowest_orbital': 12-10, 'highest_orbital': 13+10, 'num_occ_states': 6, 'num_unocc_states': 8,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.0, 'verbosity': 0, 'use_multiprocessing': False, 'nprocs': 4,
          'is_many_body': True, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../../step2_Exci_CAM_B3LYP/300K/longer/res',
          'logfile_directory': os.getcwd()+'/../../step2_Exci_CAM_B3LYP/300K/longer/all_logfiles/',
          'path_to_save_sd_Hvibs': os.getcwd()+'/newAS_6_8-new',
          'outdir': os.getcwd()+'/newAS_6_8-new','start_time': 1000, 'finish_time': 3998, 'sorting_type': 'identity',
         }

step3.run_step3_sd_nacs_libint(params_mb_sd)

