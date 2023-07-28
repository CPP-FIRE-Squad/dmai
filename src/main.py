import os
import sys

# Append the parent directory of 'src' to sys.path to enable relative imports
# current_dir = os.path.dirname(os.path.abspath(__file__))
# print(current_dir)
# parent_dir = os.path.dirname(current_dir)
# print(parent_dir)
# sys.path.append(parent_dir)

# Or use the command
        # set PYTHONPATH=/Users/cameronross/projs/galaxy_research/CameronRoss/modules
        # export PYTHONPATH=/Users/cameronross/projs/galaxy_research/CameronRoss/modules


import gizmo_analysis as gizmo

# Original way to do this
# particles = gizmo.io.Read.read_snapshots(
#         simulation_directory = '../data/subdir', snapshot_directory = '.',
#         species=['star'], snapshot_value_kind='index', snapshot_values=600)

# h = gizmo.io.Read.read_header(
#         simulation_directory = '../data/subdir', snapshot_directory = '.', 
#         snapshot_value_kind = 'index', snapshot_value = 600)



from stellarutil.simulation import Simulation, get_field, get_field_name

# simulation_directory = '/Users/cameronross/projs/galaxy_research/CameronRoss/data'
# sim = Simulation(
#     simulation_directory = '../data/old',
#     snapshot_directory = '../data/old',
#     ahf_directory = "../data/snapshot_600.z0.000.AHF_halos",
#     snapshot_value = 600,
#     snapshot_values = [0,1,2,3,4,5,6,7]
# )

sim = Simulation()
print(sim.h)


# # Print all particles you specified to keep track of
# print(sim.particles.keys())
# # Print properties that can be accessed for stars 
# print(sim.particles['star'].keys()) 
# # # Print the x pos of every star in the simulation
# print(sim.particles['star']['position'][:,0])
# # Print the n_star(64) column in the AHF file 
# print(sim.get_field('nstar')[0])
# Get a list of stars in the dark matter halo at index 0
# stars = sim.get_stars_in_halo(0, 10)
# print(f"Number of stars is {len(stars)}")
# print(stars[0].x)

# sum = 0
# for star in stars:
#     sum += star.m

# print(sum)
# print(sim.get_field('mstar')[0])