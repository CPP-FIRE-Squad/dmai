from stellarutil.simulation import Simulation, get_field, get_field_name

simulation = Simulation()

# Now you can do: 
    # simulation.h
    # simulation.particles
    # simulation.ahf_data

# Print hubble constant
print(simulation.h)
# Print the x pos of every star in the simulation
print(simulation.particles['star']['position'][:,0])


# This function matches an input string with a column name
name = get_field_name(simulation.ahf_data, 'substruct')
print(f'The name is: {name}')

# Easy way to get access to AHF column
column = get_field(simulation.ahf_data, 'substruct')
print(column)
