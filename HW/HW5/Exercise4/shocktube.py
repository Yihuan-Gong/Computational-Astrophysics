import ulula.setups.shocktube as setup_shocktube
import ulula.run as ulula_run
import ulula.simulation as ulula_sim
import ulula.plots as ulula_plt
import matplotlib.pyplot as plt
import numpy as np

# setup the Sod shocktube problem in x-direction
setup = setup_shocktube.SetupSodX()


# specify the hydro schemes
hs = ulula_sim.HydroScheme(reconstruction = 'linear', limiter = 'mc', riemann='hll', time_integration='hancock', cfl = 0.8)

# run the simulation
sim = ulula_run.run(setup, hydro_scheme=hs, tmax=0.2, nx=200)



# plot the 1D images

def line_dir():
    return 0 # xdir

def true_solu(sim, x_plot, q_plot):
    return setup.trueSolution(sim, x_plot, q_plot) # Analitcal solution

q_plot = ['DN', 'VX', 'PR']
ulula_plt.plot1d(sim, q_plot=q_plot, plot_type='line', 
                 idir_func=line_dir, true_solution_func=true_solu)
plt.savefig("shocktube_mc_linear_hancock.png")
