from stellarutil.simulation import Simulation, get_field, get_field_name

sim = Simulation()

# Print hubble constant
print(sim.h)
# Print all particles you specified to keep track of
print(sim.particles.keys())
# Print properties that can be accessed for stars 
print(sim.particles['star'].keys()) 
# Print the x pos of every star in the simulation
print(sim.particles['star']['position'][:,0])
# Print the n_star(64) column in the AHF file 
print(sim.get_field('nstar'))



# Get a list of stars in the dark matter halo at index 0
stars = sim.get_stars_in_halo(0)
print(f"Number of stars is {len(stars)}")
print(stars[0].x)

