from stellarutil import cross as cu
import astropy.io.ascii as ascii
import gizmo_analysis as gizmo
import random

# Original way to do this
particles = gizmo.io.Read.read_snapshots(
        simulation_directory = '../data',
        snapshot_directory = '../data',
        species=['star'], 
        snapshot_value_kind='index',
        snapshot_values=600
    )
h = gizmo.io.Read.read_header(
        simulation_directory = '../data',
        snapshot_directory = '../data',
        snapshot_value_kind = 'index',
        snapshot_value = 600
    )
h = h['hubble']
data = ascii.read("../snapshots/snapshot_600.z0.000.AHF_halos")
data_filtered = data[(data.field('fMhires(38)') > 0.99)]
data = data_filtered





# The last way I was doing this
particles = cu.getParticles() # Get particle from sim file
h = cu.getH() # Get hubble constant from header file
data = cu.get_ahf_data("../snapshots/snapshot_600.z0.000.AHF_halos") # get data from AHF















# Print hubble constant
print(h)
# Print the x pos of every star in the simulation
print(particles['star']['position'][:,0])






# Look at how the axis names are formated
# Will modify so units can be easily specified
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]
cu.graph(x,y, "Force vs Mass")


def getGalaxyStarInfo(data, particles, h, index = 0):
    # Get the center of the halo
    xc = data.field('Xc(6)')[index] / h
    yc = data.field('Yc(7)')[index] / h
    zc = data.field('Zc(8)')[index] / h
    # Get the peculiar velocity of halo
    vxc = data.field('VXc(9)')[index] / h
    vyc = data.field('VYc(10)')[index] / h
    vzc = data.field('VZc(11)')[index] / h

    # Get the x,y,z positions of each star particle in the simulation, and normalize it with the galaxy center
    x = particles['star']['position'][:,0] - xc
    y = particles['star']['position'][:,1] - yc
    z = particles['star']['position'][:,2] - zc
    # Get the scalefactor (age) of each star in the simulation
    a = particles['star']['form.scalefactor']
    # Get the mass of each star in the simulation
    m = particles['star']['mass']

    # Get the x,y,z velocity of each star particle in the simulation, and normalize it with the velocity center
    vx = particles['star']['velocity'][:,0] - vxc
    vy = particles['star']['velocity'][:,1] - vyc
    vz = particles['star']['velocity'][:,2] - vzc

    # Get the stars in the galaxy
    from calculations import dist
    distances =  dist(x,y,z) # Get the distance of each particle from the center galaxy
    rgal = 0.15 * data.field('Rvir(12)')[index] / h # Get the radius of the galaxy that can actually hold stars
    x_gal = x[distances < rgal] # Filter out all stars that are too far away 
    y_gal = y[distances < rgal] # Filter out all stars that are too far away 
    z_gal = z[distances < rgal] # Filter out all stars that are too far away 
    a_gal = a[distances < rgal] # Filter out all stars that are too far away
    m_gal = m[distances < rgal] # Filter out all stars that are too far away

    vx_gal = vx[distances < rgal] # Filter out all stars that are too far away
    vy_gal = vy[distances < rgal] # Filter out all stars that are too far away
    vz_gal = vz[distances < rgal] # Filter out all stars that are too far away
    v_gal = dist(vx_gal, vy_gal, vz_gal)

    return (x_gal, y_gal, z_gal, a_gal, m_gal, v_gal)

# This function returns the coordinates, velocity, mass, and age of each star in a galaxy
# Since no index is given, it defaults to 0
x,y,z,a,m,v = getGalaxyStarInfo(data, particles, h)

# Graph - Star Velocity vs Star Scale Factor
cu.graph(
    x=a, 
    y=v, 
    title="Star Velocity vs Star Scale Factor", 
    showLine=False
)





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
