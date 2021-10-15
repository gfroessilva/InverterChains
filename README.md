# String Stability in Microgrids using Frequency Controlled Inverter Chains
### Code for generating results of:

G. F. Silva, A. Donaire, M. M. Seron, A. McFadyen and J. Ford, "String Stability in Microgrids Using Frequency Controlled Inverter Chains," in IEEE Control Systems Letters, vol. 6, pp. 1484-1489, 2022, doi: 10.1109/LCSYS.2021.3114143.

The main file ```LCSS2021_InverterChains.m```:
1. sets up the system's parameters and the simulation settings
2. calls an optimisation problem through ```Optimisation_Vcte.mlx```
3. runs the simulation using ```ode15s``` and the dynamics function ```microgrids.m```
4. plots results

If you have questions, email g.froessilva@qut.edu.au.
