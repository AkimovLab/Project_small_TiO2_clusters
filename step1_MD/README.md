	# Project_relaxation_recombination_TiO2_nanoclusters

**inputs and scripts to be in the same directory to start the calculation**

**inputs**:

AIMD is performed with CP2K. The folder contains the (TiO2)1 files as case
example: geometry (`01_00.xyz`), and cp2k input (`md.inp`).


**scripts**:

The `submit.slm` file for performing step 1 calculation in a HPC.

**outputs**:

compressed file with the 3000 geometries' 16 trajectories (8 nanosystems x 2
temperatures)
To unpack:
`unzip trajectory_files.zip`

