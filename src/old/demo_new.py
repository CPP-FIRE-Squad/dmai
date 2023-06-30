from stellarutil.simulation import Simulation, get_field, get_field_name
from stellarutil.graph import graph
import random

simulation = Simulation()

# Now you can do: 
    # simulation.h
    # simulation.particles
    # simulation.ahf_data

# Print hubble constant
print(simulation.h)
# Print the x pos of every star in the simulation
print(simulation.particles['star']['position'][:,0])

# Look at how the axis names are formated
# Will modify so units can be easily specified
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]
graph(x,y, "Force vs Mass")


# This function returns the coordinates, velocity, mass, and age of each star in a galaxy
# Since no index is given, it defaults to 0
x,y,z,a,m,v = simulation.get_star_info()
# Graph - Star Velocity vs Star Scale Factor
graph(
    x=a, 
    y=v, 
    title="Star Velocity vs Star Scale Factor", 
    showLine=False
)

# This function matches an input string with a column name
name = get_field_name(simulation.ahf_data, 'substruct')
print(f'The name is: {name}')

# Easy way to get access to AHF column
column = get_field(simulation.ahf_data, 'substruct')
print(column)

# Let's graph random stuff for fun
data = simulation.ahf_data
for i in range(0,3):
    ran_1 = random.randint(1, 83)
    ran_2 = random.randint(1, 83)
    graph(
        x = get_field(data, ran_1), 
        y = get_field(data, ran_2), 
        title = f"{get_field_name(data, ran_2)} vs {get_field_name(data, ran_1)}",
        showLine =  False,
        windowTitle=f"Random Figure {i+1}"
    )
