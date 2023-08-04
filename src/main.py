# import gizmo_analysis as gizmo


# particles = gizmo.io.Read.read_snapshots(
#         simulation_directory = 'm10r_res250md', snapshot_directory = 'output',
#         species=['star'], snapshot_value_kind='index', snapshot_values=600)

# h = gizmo.io.Read.read_header(
#         simulation_directory = '..data2/m10r_res250md', snapshot_directory = 'output', 
#         snapshot_value_kind = 'index', snapshot_value = 600)

# h = h['hubble']
# print(h)



# from stellarutil.simulation import Simulation, get_field, get_field_name
# sim = Simulation()
# print(sim.h)
# # simulation_directory = '/Users/cameronross/projs/galaxy_research/CameronRoss/data'
# # sim = Simulation(
# #     simulation_directory = '../data/old',
# #     snapshot_directory = '../data/old',
# #     ahf_directory = "../data/snapshot_600.z0.000.AHF_halos",
# #     snapshot_value = 600,
# #     snapshot_values = [0,1,2,3,4,5,6,7]
# # )


from stellarutil.simulation import Simulation
from stellarutil.graph import graph
sim = Simulation()
print(sim.h)

# Mvir vs Mstar
graph(sim.get_field('Mvir'), sim.get_field('Mstar'), "Mvir vs Mstar")

# x vs y of all stars and color code by scale factor

# 2d histogram of mass of stars

# graph two galaxies on the same 3d plot next to each other



