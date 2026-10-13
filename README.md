# JuliaMD

This set of programs constitute a molecular dynamics environment based on [Molly.jl](https://github.com/JuliaMolSim/Molly.jl). JuliaMD contains utilities for generating a topology and detecting atom types using bonding criteria. It can also load a forcefield from an XML file (only GAFF has been tested so far), and assign the forcefield parameters to the appropriate bonds, angles and dihedrals, based on the GAFF atom types previously detected. It can also detect certain aromatic rings, and it can assign charges using the Charge Equilibration method [1,2]

#Installation

Requirements: [Molly.jl](https://github.com/JuliaMolSim/Molly.jl), [EzXML.jl](https://github.com/JuliaIO/EzXML.jl), [Chemfiles.jl](https://github.com/chemfiles/Chemfiles.jl)

Download the contents of the folder and start julia from this folder, or alternatively add the folder to the appropriate path. The file buildMoleculeTopology.jl contains the main program, where some example structures are already defined, as well as some simulation options using Molly. This program can be modified by the user. In the future I will standardize an input file format to hide the complexity of the code. There are other programs with more limited functionality, such as fitFF2traj.jl, which is not finished, but is meant to analyze orca trajectories with a graphic interface.



1. Rappe, Anthony K., and William A. Goddard III. "Charge equilibration for molecular dynamics simulations." The Journal of Physical Chemistry 95, no. 8 (1991): 3358-3363.
2. Oda, Akifumi, and Shuichi Hirono. "Geometry-dependent atomic charge calculations using charge equilibration method with empirical two-center Coulombic terms." Journal of Molecular Structure: THEOCHEM 634, no. 1-3 (2003): 159-170.
