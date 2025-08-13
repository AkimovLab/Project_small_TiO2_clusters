# Project_relaxation_recombination_TiO2_nanoclusters

**scripts and selected outputs from computing the nonadiabatic couplings (NACs) between pairs of TD-DFT states.**

**scripts**:

The inputs in this case are part of the outcome from step2, so the required files are invoked
in the `step3_NACs.py` file. It is not needed to move such files to this directory.
There is not a "input" folder in this case.


We add the script `PD_NACs_Egaps.py` that (after conducting all step3
calculations) generate the plots for (at both 100 K and 300 K) probability distribution (PD) of 
all NACs between ground state (S0) and each of the first three excited states (S1-S3); NACs between 
all adjacent (i, i+1) states excluding S0; energy gaps between ground and first excited state; and 
energy gaps between all adjacent (i, i+1) states excluding S0.


**key_outputs**:

The outputs from this step (the NACs) are heavy files so are not included here. They
are easy to obtain by running the step3_NACs.py file with Python or calling it with
the proper submit (`submit_template.sh`) file. We have included PNG files
with the NACs and energy gaps PDs for all systems obtained running the
`PD_NACs_Egaps.py` script.



