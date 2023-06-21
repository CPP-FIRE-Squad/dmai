from stellarutil import cross as cu
import random

particles = cu.getParticles() # Get particle from sim file
h = cu.getH() # Get hubble constant from header file
data = cu.get_ahf_data("../snapshots/snapshot_600.z0.000.AHF_halos") # get data from AHF


# Clear the screen and print list of libraries
cu.clear_console()
cu.list_libraries()

# Look at how the axis names are formated
# Will modify so units can be easily specified
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]
cu.graph(x,y, "Force vs Mass")


# This function returns the coordinates, velocity, mass, and age of each star in a galaxy
# Will be cleaned up later
x,y,z,a,m,v = cu.getGalaxyStarInfo(data, particles, h)



# Graph - Star Velocity vs Star Scale Factor
cu.graph(x=a, y=v, title="Star Velocity vs Star Scale Factor", showLine=False)


# This function matches an input string with a column name
name = cu.getFieldName(data, 'substruct')
print(f'The name is: {name}')

# Easy way to get access to AHF column
column = cu.getField(data, 'substruct')
print(column)


# Let's graph random stuff for fun
for i in range(0,3):
    ran_1 = random.randint(1, 83)
    ran_2 = random.randint(1, 83)
    cu.graph(
        x = cu.getField(data, ran_1), 
        y = cu.getField(data, ran_2), 
        title = f"{cu.getFieldName(data, ran_2)} vs {cu.getFieldName(data, ran_1)}",
        showLine =  False,
        windowTitle=f"Random Figure {i+1}"
    )

cu.help()


