**inputs and scripts to be in the same directory to start the calculation**

**inputs**:

`NAMD_FSSH.py`; `NAMD_DISH.py`; `NAMD_MSDM.py`: scripts for performinig the NA-MD runs following different TSH schemes.

`recipes`: Various combinations of parameters, chosen to represent particular known method. In principle, don't modify those.

The `submit_*.slm` files are used to run the `NAMD_*` files in a HPC. 

**outputs**:

All output files are included in the two compressed files..

To unpack, decompress each file separately:

    tar -xvjf `step4_FSSH_TiO2_2_100K.tar.bz2`
    tar -xvjf `step4_DISH_TiO2_2_100K.tar.bz2`


This directory also contains:

`GS_populations.py`— script for plotting the ground-state population evolution for both FSSH and DISH results.
`single_exponential_decay.py`; `bi_exponential_decay.py`; `stretched_exponential_decay.py`- 3 scripts that read data from Step 3 and Step 4, plots the average excitation energy decay to S1 and 
fits the decay with three functions (Stretched exponential, Bi-exponential, Single exponential). These scripts plot the raw average energy (cyan) with ±1σ standard deviation across trajectories (gray), smoothed average (red), and fit (black dashed). Horizontal lines mark the initial energy (green) and the population-weighted S₁ energy (blue), toward which the 
system decays. They are prepared for invoking just the DISH results but can be easily modified to account for FSSH evolution as well.
Also includes scripts for plotting the energy gap and NAC probability distributions for two cases: (1) between the ground state and the first excited state (related to ground-state recombination) and (2) between all excited states (related to relaxation from a higher-energy initial state).

