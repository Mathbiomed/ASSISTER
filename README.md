# ASSISTER (Adaptive Simplification of StochastIc SystEm with Reversible binding)

This is a matlab code for efficient and accurate simulations of stochastic biochemical systems containing rapid reversible binding reactions. Detailed step-by-step manual will be uploaded soon.

## Code Description
1. Gillespie_Reduction.m
> Main function for this package. This function performs accurate and efficient stochastic simulations of models containing rapid reversible 
bindings, with the adaptive choice between two exclusively valid non-elementary propensities, the stochastic total quasi-steady-state approximation (stQSSA), and the stochastic low-state quasi-steady-state approximation (slQSSA). The below two functions are the auxiliary functions for this main function.

2. QSSA_Threshold.m
> This function determines which of the stQSSA and the slQSSA is valid and calculates the number of needed states (L) for the slQSSA to ensure the smaller error than the tolerance (epsilon) for given Kd value.

2. LQSSA.m
> This function calculates the L-state slQSSA of A in the reversible binding reaction (A + B <=> C) for given AT, BT, Kd, and L, where AT = A + C, BT = B + C, and Kd is the normalized dissociation constant (product of the system volume and the dissociation constant) of the reversible binding.

4. test.m
> This file contains an example code that runs the main function, Gillespie_Reduction.

5. Codes_for_figs.zip
> The matlab codes used for drawing figures in the manuscript and the supporting information.