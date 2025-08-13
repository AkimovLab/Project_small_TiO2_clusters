# Project_relaxation_recombination_TiO2_nanoclusters

**scripts and selected outputs from performing the nonadiabatic molecular dynamics (NA-MD) runs.**

**scripts**

Again, the inputs are part of the outcome from the previous steps. No "inputs" folder here.

`NAMD_*.py` the script for starting the NA-MD using one trajectory surface hoping
(TSH) method. In our case, we have applied FSSH, and the 2023 revised version
of DISH.

recipes: This folder contains the pre-defined recipes for the different TSH schemes. Don't modify
this folder, just leave in the same directory where the NA-MD runs are conducted.

To start the calculations, simply run the `NAMD_*py` file with Python or call it with
the submit (`submit_*`) file.

After finishing the NA-MD runs, the script `GS_populations.py` is used to plot the Ground State population
evolution (5-molecules TiO2 system as example case). Besides, the included `average_excitation_energy_decay_3functions.py`
plots the average excitation energy decay to S1 fitting it with 3 different functions (stretched exponential, bi-exponential,
and single exponential). (TiO2)5 used again as example case.

**key_outputs**

The outputs from this step (evolution of excited state populations) are heavy files so are not included here. They
are easy to obtain by running the `NAMDs_*.py` files with Python or calling it with
the proper submit (`submit_template.sh`) file. We have included PNG files with the GS population evolutions for all
systems and both 100 and 300 K temperatures; and the same for the average excitatione energy decay.


