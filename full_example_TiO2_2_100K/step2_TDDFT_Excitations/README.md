# Project_relaxation_recombination_TiO2_nanoclusters

**inputs and scripts to be in the same directory to start the calculation**

**inputs**:

Aligned AIMD trajectory for this (TiO2)2 system at 100 Ke (`dynamcis_AIMD-aligned_100K.xyz`).

CP2K input where TD-DFT options are detailed (`es_diag_temp.inp`).

**scripts**:

`distribute_jobs.py` is the file that you need modify by specifying initial and final step of the trajectory and the number of jobs. Libra will split the trajectory based on these values and will submit them by creating the specific folders for each job.

`run_template.py`: A file that contains the data to run the MO overlap calculations. It is important to properly define the Active Space. In this case we have chosen orbitals going from 24-15 to 25+15. 24th and 25th orbitals are respectively the HOMO and the LUMO.

The `submit_template.slm` file is the one which will be submitted and runs the python `run.py` where the `run.py` file is the copy of the `run_template.py` file but with the initial and final steps filled based on the requested number of jobs and the initial and final step for the trajectory.

`vmd_input_template`: The VMD input template. This input will be modified by Libra but only the mol load cube and rendering line, render Tachyon ... will be changed. You can make your own input template by opening the VMD on your computer and from File choose Log Tcl commands to file... and it will append all the commands for reproducing that image into a file.

To start the TD-DFT calculation simply run: **python distribute_jobs.py**



**outputs**:

All logfiles with TD-DFT CP2K output information are included in `logfiles_TiO2_2_100K.tar.bz2`.
`res_TiO2_2_100K.tar.bz2`contains all KS energies (E_ks_*); MO overalps (S_ks_*) and time-overlaps (St_ks_*).

To unpack, decompress each file separately:

    tar -xvjf logfiles_TiO2_2_100K.tar.bz2
    tar -xvjf res_TiO2_2_100K.tar.bz2


