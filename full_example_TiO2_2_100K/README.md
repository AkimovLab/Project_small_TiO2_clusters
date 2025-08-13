	#full_example_TiO2_2_100K

**scripts and outputs from performing the complete four-step workflow for the (TiO₂)₂ nanocluster case.**


The provided inputs and scripts correspond to the four sequential steps of the workflow. 
The scripts are prepared to be run like:

**python `name_of_the_script`**

No further action is needed.


**step1_MD**

Contains the trajectory file obtained with the molecular dynamics (MD) simulation. No extra action is required here.

**step2_TDDFT_Excitations**

This folder contains two compressed files:

all_logfiles.tar.bz2

res.tar.bz2

To prepare this step, decompress both files:

tar -xvjf all_logfiles.tar.bz2
tar -xvjf res.tar.bz2

Ensure the extracted folders are named:

all_logfiles
res

**step3_NACs**

Contains multiple compressed files with nonadiabatic couplings and energy gap data.
Also includes scripts for plotting the energy gap and NAC probability distributions for two cases: (1) between the ground state and the first excited state (related to ground-state recombination) and (2) between all excited states (related to relaxation from a higher-energy initial state).

To prepare this step, decompress all .tar.bz2 files:

for f in *.tar.bz2; do tar -xvjf "$f"; done

Create a new folder called results and move all deompressed contents into it:

mv */* results/ 2>/dev/null
rmdir */ 2>/dev/null

**step4_NAMD**

Contains two compressed files:

step4_FSSH_TiO2_2_100K.tar.bz2
step4_DISH_TiO2_2_100K.tar.bz2

Also includes:

`GS_populations.py`— script for plotting the ground-state population evolution for both FSSH and DISH results.
`average_excitation_energy_decay_*.py` - 3 scripts that read data from Step 3 and Step 4, plots the average excitation energy decay to S1 and 
fits the decay with three functions (Stretched exponential, Bi-exponential, Single exponential). These scripts plot the raw average energy (cyan) with ±1σ standard deviation across trajectories (gray), smoothed average (red), and fit (black dashed). Horizontal lines mark the initial energy (green) and the population-weighted S₁ energy (blue), toward which the 
system decays. They are prepared for invoking just the DISH results but can be easily modified to account for FSSH evolution as well.

To prepare this step, decompress each file separately:

tar -xvjf step4_FSSH_TiO2_2_100K.tar.bz2
tar -xvjf step4_DISH_TiO2_2_100K.tar.bz2

Rename the extracted folders:

mv step4_FSSH_TiO2_2_100K FSSH_results
mv step4_DISH_TiO2_2_100K DISH_results
