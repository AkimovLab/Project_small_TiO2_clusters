# Project\_relaxation\_recombination\_TiO2\_nanoclusters



This project methodology is composed by 4 sequential steps: (i) generation of the nuclear guiding trajectory through ab initio Molecular
Dynamics (AIMD); (ii) computation of the TD-DFT excitations, giving rise to the excitation energies, Molecular Orbital overlaps and MO time-overlaps; (iii) nonadiabatic couplings (NACs) computation; and (iv) running the nonadiabatic molecular dynamics (NA-MD). In this way, the four directories contain the inputs and scripts required to run each step + some key outputs and extra scripts for visualizing some of the results.



We recommend visiting https://github.com/compchem-cybertraining/Tutorials\_Libra (specially the /master/6\_dynamics/2\_nbra\_workflows branch) for getting detailed
information about how to successfully perform the different steps.



Along with the different inputs and scripts for all four sequential steps, we provide a folder (full\_example\_TiO2\_2\_100K) containing the complete results for one representative case, which serves as a test system. This folder also includes several visualization scripts to plot the NAC and energy gap probability distributions, the ground-state population, and the average excitation energy over time—covering all properties analyzed in this project. Additional details are available in the README file within this directory. We recommend following this example rather than running the visualization scripts within each of the 4-step folders (as they require generating the proper results first).

