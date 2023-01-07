# IGNM 
is a package for solving an ignition model that combines solid-phase (pyrolysis) and gas-phase combustion in a one-dimensional domain, suitable for reproducing cone calorimeter (and similar) experiments. The model features detailed chemical kinetics and relies on [Cantera](https://cantera.org) for evaluating the thermodynamic properties.

The function `gasPhaseSolver.solve` takes as input, the time-series of temperature `Tin`, mass flux `Jin`, mass fraction `Yin` at the inlet and solves the gas-phase combustion. 