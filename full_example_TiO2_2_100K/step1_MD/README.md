# Project_relaxation_recombination_TiO2_nanoclusters

**inputs and scripts to be in the same directory to start the calculation**

**inputs**:

AIMD is performed with CP2K. The folder contains the relaxed (TiO2)2 geometry at 100 K as case
example (`pos.xyz`), and cp2k input (`md.inp`) for performing the MD run.


**scripts**:

The `submit.slm` file for performing step 1 calculation in a HPC.

**outputs**:

3000 steps AIMD generated trajectory: `trajectory_TiO2_2_100K.xyz`

-------------

To prepare this step, combine all input and script files in the same folder and run:

    sbatch submit.slm
