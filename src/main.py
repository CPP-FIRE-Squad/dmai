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
from stellarutil.graph import graph, star_scatter_plot, histogram, stars_scatter_plot

# sim = Simulation(
#     simulation_directory = '../data/m10r_res250md',
#     snapshot_directory='output',
#     ahf_path='../data/m10r_res250md/snapshot_600.AHF_halos'
# )

sim = Simulation('m10r_res250md')

# Print hubble constant
print(sim.h)

# Mvir vs Mstar
# graph(sim.get_field('Mvir'), sim.get_field('Mstar'), "Mvir vs Mstar", showLine=False)

# x vs y of all stars and color code by scale factor
stars = sim.get_stars_in_halo()
x_positions = [star.x for star in stars]
y_positions = [star.y for star in stars]
z_positions = [star.z for star in stars]
ages = [star.a for star in stars]
# star_scatter_plot(x_positions, y_positions, z_positions, ages, .2)

# star was histogram
data = [star.m for star in stars]
# histogram(data, bins = 10, title='Mass Distribution', x_label='Mass') 

# graph two galaxies on the same 3d plot next to each other
stars2 = sim.get_stars_in_halo(1)
xc1 = sim.get_field('Xc(6)')[1] / sim.h - sim.get_field('Xc(6)')[0] / sim.h
yc1 = sim.get_field('Yc(7)')[1] / sim.h - sim.get_field('Yc(7)')[0] / sim.h
zc1 = sim.get_field('Zc(8)')[1] / sim.h - sim.get_field('Zc(8)')[0] / sim.h

for star in stars2:
    star.x = star.x - xc1
    star.y = star.y - yc1
    star.z = star.z - zc1

stars_scatter_plot(stars, stars2)