from stellarutil.simulation import Simulation

# This line does all the nasty set up for you.
sim = Simulation(
    simulation_directory = '../data',
    snapshot_directory = '../data',
    ahf_path = "../data/snapshot_600.z0.000.AHF_"
)

# Print hubble constant
print(sim.h) 
# Print all particles you specified in the line with "Take Note"
print(sim.particles.keys())
# Print properties that can be accessed for stars 
print(sim.particles['star'].keys()) 
# Print the x pos of every star in the simulation
print(sim.particles['star']['position'][:,0])
# Print the number of stars in each halo
print(sim.get_field("n_star"))