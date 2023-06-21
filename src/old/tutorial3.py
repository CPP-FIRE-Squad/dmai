from stellarutil import cross as cu
import random

# Get the data tables and constant we will be working with
particles = cu.getParticles()
h = cu.getH() 
data = cu.get_ahf_data("snapshot_600.z0.000.AHF_halos") # filter out anything not in that 1% - gotta ask what is happening here again

def getGalaxyStarInfo(index = 0):
    # Get the center of the given galaxy in the simulation
    xc = data.field('Xc(6)')[index] / h
    yc = data.field('Yc(7)')[index] / h
    zc = data.field('Zc(8)')[index] / h

    # Get the center of the given galaxy in the simulation
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
    distances = cu.dist(x,y,z) # Get the distance of each particle from the center galaxy
    rgal = 0.15 * data.field('Rvir(12)')[index] / h # Get the radius of the galaxy that can actually hold stars
    x_gal = x[distances < rgal] # Filter out all stars that are too far away 
    y_gal = y[distances < rgal] # Filter out all stars that are too far away 
    z_gal = z[distances < rgal] # Filter out all stars that are too far away 
    a_gal = a[distances < rgal] # Filter out all stars that are too far away
    m_gal = m[distances < rgal] # Filter out all stars that are too far away

    vx_gal = vx[distances < rgal] # Filter out all stars that are too far away
    vy_gal = vy[distances < rgal] # Filter out all stars that are too far away
    vz_gal = vz[distances < rgal] # Filter out all stars that are too far away
    v_gal = cu.dist(vx_gal, vy_gal, vz_gal)

    return (x_gal, y_gal, z_gal, a_gal, m_gal, v_gal)



cu.clear_console()

# Graph 1, Graph 2 - Two graphs of all stars in a galaxy, having color represent age
star_counts = cu.getField(data, 'n_star')
print(star_counts)
x,y,z,a,m,v = getGalaxyStarInfo()
x2,y2,z2,a2,m2,v2 = getGalaxyStarInfo(2)
cu.star_scatter_plot(x,y,z,a, 1)
cu.star_scatter_plot(x2,y2,z2,a2,1)


# Graph 3 - pie chart show all galaxies that have at least 3 stars
labels = ['Halo {} - {}'.format(i, count) for i, count in enumerate(star_counts) if count > 3]
cu.pie_chart(star_counts, labels) 


# Graph 4 - to demo function
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]
cu.graph(x,y, "Force vs Acceleration")


# Graph 5 - M_star vs M_vir (mass of all stars vs mass of halo)
M_vir = cu.getField(data, "M_vir")
M_star = cu.getField(data, "M_star")
cu.graph(x = M_vir, y = M_star, title = "M_star vs M_vir", showLine=False)

# Graph 6 - R_vir vs R_max (virial radius vs position of rotational curve maximum)
R_max = cu.getField(data, "R_max")
R_vir = cu.getField(data, "R_vir")
cu.graph(x = R_max, y = R_vir, title = "R_vir vs R_max", showLine=False)

# Graph 7 - V_max vs M_vir (max of rotational curve vs mass of halo)
M_vir = cu.getField(data, "M_vir")
V_max = cu.getField(data, "V_max")
cu.graph(x = M_vir, y = V_max, title = "V_max vs M_vir", showLine=False)


# Graph 8-11 - Random one cuz why not
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



x,y,z,a,m,v = getGalaxyStarInfo()
# Graph 12 - Star Mass vs Star Scale Factor
cu.graph(x=a, y=m, title="Star Mass vs Star Scale Factor", showLine=False)

# Graph 13 - Star Velocity vs Star Scale Factor
cu.graph(x=a, y=v, title="Star Velocity vs Star Scale Factor", showLine=False)

