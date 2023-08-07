import gizmo_analysis as gizmo

# particles = gizmo.io.Read.read_snapshots(
#         simulation_directory = '../data/m10r_res250md', 
#         snapshot_directory = 'output',
#         species=['star'], 
#         snapshot_value_kind='index', 
#         snapshot_values=600
#     )

# h = gizmo.io.Read.read_header(
#         simulation_directory = '../data/m10r_res250md', 
#         snapshot_directory = 'output', 
#         snapshot_value_kind = 'index', 
#         snapshot_value = 600
#     )

# h = h['hubble']
# print(h)



from stellarutil.simulation import Simulation
from stellarutil.graph import graph

# sim = Simulation(
#     simulation_directory = '../data/m10r_res250md',
#     snapshot_directory='output',
#     ahf_directory='../data/m10r_res250md/snapshot_600.AHF_halos'
# )

sim = Simulation('m10r_res250md')
# Print hubble constant
print(sim.h)
# Mvir vs Mstar
graph(sim.get_field('Mvir'), sim.get_field('Mstar'), "Mvir vs Mstar", showLine=False)
# x vs y of all stars and color code by scale factor
stars = sim.get_stars_in_halo()
x_positions = [star.x for star in stars]
y_positions = [star.y for star in stars]
graph(x_positions, y_positions, "X vs Y", showLine=False)

# 2d histogram of mass of stars
# graph two galaxies on the same 3d plot next to each other



