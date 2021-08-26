# LIOUSHELL
Computing true polar wander on dynamic planets

To compile and run the code simply type ./gun.bat in terminal in the folder with source code. Gfortran is chosen as the default compiler and it is free (the GNU Fortran compiler). If you want the executable to have a specific name and run on background, type ./gun.bat myname. For compilation with commercial ifort, type ./gun.bat myname ifort. However, gfortran compilation seems to provide the same speed upon execution. By default, LIOUSHELL uses openMPI and runs on 3 threads, remove the -fopenmp compilation flag for a serial run. To get basic figures for your simulation, type ./postpro.bat when your simulation is completed.

Two files are important to run a simulation: param.in and inistate.dat

param.in specifies input parameters used by the code. For examples of the full Liouville equation solutions described in "True Polar Wander on Dynamic Planets: Approximative Methods vs. Full Solution", copy param.ES (Fig. 1), param.EF (Fig. 2), or param.VP (Fig. 6) from subfolder param_files as param.in into the source code directory. See the default param.in in the source directory for a brief description of each of the typical input parameteres. 

Typically, each simulation is loaded from the hydrostatic shape of the studied body, stored in inistate.dat. Inistate.dat is generated by LIOUSHELL when savestate=1 (and loadstate=0, see the commented param.in). To spin a given model with a constant angular speed (and thus to obtain the hydrostatic state), run a simulation with loadstate=0, couple=0, and tstart= the expected relaxation time. For typical series of TPW simulations, the same precomputed hydrostatic shape loaded via loadstate=1 for all runs from the series. In subfolder inistates, the hydrostatic shapes for the Earth and Venus models from "True Polar Wander on Dynamic Planets: Approximative Methods vs. Full Solution" are included. Copy these as inistate.dat into the source directory in order to use it. For loadstate=1 simulations, the corresponding urad_ini.tail file should be placed on the same level as the source code directory - this file is needed if the energy components are to be plotted relative to the hydrostatic state only (see postpro.bat). Urad_ini.tail is the last line of urad_ini.dat that is generated when loadstate=0, i.e. that is generated during the initial relaxation process.

