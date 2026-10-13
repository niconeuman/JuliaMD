# JuliaMD

This set of programs constitute a molecular dynamics environment based on [Molly.jl](https://github.com/JuliaMolSim/Molly.jl). JuliaMD contains utilities for generating a topology and detecting atom types using bonding criteria. It can also load a forcefield from an XML file (only GAFF has been tested so far), and assign the forcefield parameters to the appropriate bonds, angles and dihedrals, based on the GAFF atom types previously detected. It can also detect certain aromatic rings, and it can assign charges using the Charge Equilibration method [1,2]





1. Rappe, Anthony K., and William A. Goddard III. "Charge equilibration for molecular dynamics simulations." The Journal of Physical Chemistry 95, no. 8 (1991): 3358-3363.
2. Oda, Akifumi, and Shuichi Hirono. "Geometry-dependent atomic charge calculations using charge equilibration method with empirical two-center Coulombic terms." Journal of Molecular Structure: THEOCHEM 634, no. 1-3 (2003): 159-170.
