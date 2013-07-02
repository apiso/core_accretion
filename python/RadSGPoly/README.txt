This folder contains the following scripts:

atmseries_poly.py		Used to generate evolutionary profiles for a series of cores at a fixed distance
atmseries_script.py		Used with 'execfile' to execute atmseries_poly at several distances
cooling_poly.py		Used to calculate cooling terms + cooling time
plot_errmap.py		Andrew's script for 2D model
plot_helper.py		Used to calculate crossover time at fixed distance + critical core mass at various AU. Generally imported in python notebooks
plot_profs.py			Andrew's script for 2D model
plot_shoot2.py		Andrew's script for 2D model
profiles_poly.py		Used to generate evolutionary profiles for a given a and Mc
shoot2.py			Andrew's script for 2D model
shooting_poly.py		Used to generate one static profiles for a given a, Mc and initial atmosphere mass
TPeffects.py			Old script used to explore the separate effects of disk T and P. No longer used.

and the following ipython notebooks:

TPeffects.ipynb			Generates plots that show the separate effects of disk T and P on crossover time
cooling_fixed_disklife.ipynb	Generates t vs M plots at fixed a, Mcrit vs a, etc. for different mu and delad
cooling_terms_poly.ipynb	Generates the cooling terms plot for the paper
delad_checks.ipynb		Old notebook plotting profiles. No longer used.
opacity_effects.ipynb		Plots L, t vs M for opacity reduced by 10 and 100
plots_poster.ipynb		Instantaneous and evolutionary plots for various a, Mc, numerical + analytic
stats_vs_evolve.ipynb		Andrew's notebook for 2D model
