# cell-migration-scaffold
MATLAB code to simulate T-cell migration in microprinted 3D scaffolds and extracellular matrix.

The file example_trajectories.m generates different examples of scaffolds, simulate the cell migration inside the scaffold with a bias in the direction params.d and plots the results. The parameters are listed in the file params_reference.

The main function to simulate cell migration is helpers/movement_scaffold.m. Note that it works for any scaffold saved in .stl, as shown in the last example of example_trajectories.m.
