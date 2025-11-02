# curved_shell_dynamics
Matlab code from my project work on non-linear vibrations in shell structures.

# Symbolic Derivations
The scripts analytically derive equations of motion for thin, doubly-curved shell elements with simply supported boundary conditions using different shell formulations as well as numbers of degrees of freedom, i.e. modes, with and without inertial coupling between modes. They also compute the Jacobians for each set of EOMs.

# Numerical Simulation
The numerical simulations have been performed in Julia for performance reasons. The results were exported to csv-files (not included here).

# Data Analysis
Scripts for data analysis are provided that parse the numerical results and create backbone-plots, bifurcation diagrams and orbit diagrams.
