# Molecular-Mechanics
This is the project for the VU/UvA master course 'Scientific modelling and programming'. The topic is molecular mechanics and the minimization of the energy function of a certain molecule.

To use the code, use the command:

gfortran DataModule.f90 Calculations.f90 MinimizationModule.f90 EnergyModule.f90 main.f90

To change the molecule that you want to investigate, change the file name in DataModule.f90 in line 47 to the name of another input file.
