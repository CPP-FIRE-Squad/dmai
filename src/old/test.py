import stellarutil.simulation as sim

simulation = sim.Simulation()
# Print hubble constant
print(simulation.h)
# Print the x pos of every star in the simulation
print(simulation.particles['star']['position'][:,0])
