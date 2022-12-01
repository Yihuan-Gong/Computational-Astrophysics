import ulula.setups.sedov_taylor as setup_sedov
import ulula.run as ulula_run
import ulula.simulation as ulula_sim
import ulula.plots as ulula_plt
import matplotlib.pyplot as plt
import numpy as np

# setup the Sod shocktube problem in x-direction
setup = setup_sedov.SetupSedov()

# specify the hydro schemes
hs = ulula_sim.HydroScheme(reconstruction = 'const', limiter = 'minmod', riemann='hll', time_integration='euler', cfl = 0.8)

# run the simulation
sim = ulula_run.run(setup, hydro_scheme=hs, tmax=0.02, nx=200)

# plot the 2D images
q_plot = ['DN','PR']
ulula_plt.plot2d(sim, q_plot=q_plot)
plt.savefig("sedov_img.png")

# TODO: plot the 2D profiles for density, pressure, and total velocity
