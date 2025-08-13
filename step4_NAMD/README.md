# Project_relaxation_recombination_TiO2_nanoclusters

**scripts and selected outputs from performing the nonadiabatic molecular dynamics (NA-MD) runs.**

**scripts**

Again, the inputs are part of the outcome from the previous steps. No "inputs" folder here.

`NAMD_*.py` the script for starting the NA-MD using one trajectory surface hoping
(TSH) method. In our case, we have applied FSSH, and the 2023 revised version
of DISH.

recipes: This folder contains the pre-defined recipes for the different TSH schemes. Don't modify
this folder, just leave in the same directory where the NA-MD runs are conducted.


To prepare this step, make sure the previous setps have been properly done, check the path from step 3 outputs is well-written in `NAMD_*.py`and run:

    sbatch submit_*

After finishing the NA-MD runs, the script `GS_populations.py` is used to plot the ground state population evolution, as deailed above. The `average_excitation_energy_decay_3functions.py` is sed for fitting the average excitation energy decay to S1, but it requires a lot of heavy files which are not included here (see below).

**key_outputs**

The outputs from this step (evolution of excited-state populations) are provided as compressed files, organized by system, trajectory surface hopping (TSH) scheme, and temperature.

The ground-state population evolution (associated with recombination) can be plotted by running the script `GS_populations.py`. Stretched and Gaussian fitting functions are used for the FSSH and DISH methods, respectively.

The decay from the initial state to S1 is harder to simulate because it requires the results from step 3, which are too large to include here. However, we have included one complete case—with all outputs and inputs—as a test example in the **full_example_TiO2_2_100K folder**. Please see that case for reference.
