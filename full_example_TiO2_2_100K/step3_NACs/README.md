**inputs and scripts to be in the same directory to start the calculation**

**inputs**:

`step3_NACs.py` is the key file in which again the Active Space has to be defined (consistent with the previous step), as well as the number of states.

The `submit_template.slm` file is the one which will be submitted and runs the python `step3_NACs.py` 


**outputs**:

This directory contains multiple compressed files with nonadiabatic couplings and energy gap data, generated after ths step is properly run.
Also includes scripts for plotting the energy gap and NAC probability distributions for two cases: (1) between the ground state and the first excited state (related to ground-state recombination) and (2) between all excited states (related to relaxation from a higher-energy initial state).

All output files are included in the various `step3_TiO2_100K_part*.tar.bz2`.

To unpack, decompress each file separately:

    tar -xvjf `step3_TiO2...`
