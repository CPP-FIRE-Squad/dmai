from stellarutil.simulation import Simulation
from stellarutil.graph import graph, star_scatter_plot, histogram, stars_scatter_plot

m10r = Simulation(simulation_name='m10v_res030md', species=['star'])
m10r = Simulation('m10v_res030md', species=['star'])
# sim = Simulation(
#     simulation_directory = '../data/m10r_res250md',
#     snapshot_directory='output',
#     ahf_path='../data/m10r_res250md/snapshot_600.AHF_halos',
# )

# Print hubble constant
print(m10r.h)
# Print field
print(m10r.get_field('nstar'))

# # Mvir vs Mstar
graph(m10r.get_field('Mvir'), m10r.get_field('Mstar'), "Mvir vs Mstar", showLine=False, logx=True, logy=True)
# position of all stars and color code by scale factor
halo = m10r.get_halo()
ages = [star.a for star in halo.stars]
star_scatter_plot(halo.stars, ages)
# star mass histogram
masses = [star.m for star in halo.stars]
histogram(masses, bins = 10, title='Mass Distribution', x_label='Mass') 
# graph two galaxies on the same 3d plot next to each other
halo2 = m10r.get_halo(1)
m10r.center_on(halo2, 0)
stars_scatter_plot(halo.stars, halo2.stars)