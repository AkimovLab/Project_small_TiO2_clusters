# Project_relaxation_recombination_TiO2_nanoclusters

**scripts and selected outputs from computing the nonadiabatic couplings (NACs) between pairs of TD-DFT states.**

**scripts**:

The inputs in this case are part of the outcome from step2, so the required files are invoked
in the `step3_NACs.py` file. It is not needed to move such files to this directory.
There is not a "input" folder in this case.


To prepare this step, make sure the previous setps have been properly done, check the path from step 2 outputs is well-written in `step3_NACs.py`and run:

**sbatch `submit_template.sh`**

We add the script `PD_NACs_Egaps.py` that (after conducting all step3
calculations) generate the plots for (at both 100 K and 300 K) probability distribution (PD) of 
all NACs between ground state (S0) and each of the first three excited states (S1-S3); NACs between 
all adjacent (i, i+1) states excluding S0; energy gaps between ground and first excited state; and 
energy gaps between all adjacent (i, i+1) states excluding S0.

We add the `Egaps_NACs_PD_visualization.py` which use is explained below.

**key_outputs**:

The outputs from this step (the NACs) are large files and are therefore not included here.
You can generate them by running step3_NACs.py with Python or by using the submission script (`submit_template.sh`), as described above.

For convenience, we have included .txt files containing the minimal data required to plot the PD function for:

The energy gap evolution between adjacent states (besides 0-1), and between states 0 and 1.

The NACs for the same cases but the latter goes from 0 to 1-2-3.

The files are named (1 per temperature value):

gap_adjacent_*, gap_01_*, nac_adjacent_*, and nac_0_123_*. 

To create the corresponding plots, go to the **scripts** folder and run:

**python `Egaps_NACs_PD_visualization.py`**

Before running, edit the script to select the desired property for the PD function.


