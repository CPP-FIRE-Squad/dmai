from stellarutil import cross as cu
import numpy

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

# subtract x_c, restricvt @ 40%, get the mean of x,y,z, subtract the mean_x to every element in x and so on, then restrict @ 15% again
def getGalaxyStarInfo2(index = 0):
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

    rgalStar1 = 0.40 * data.field('Rvir(12)')[index] / h # Get the radius of the galaxy that can actually hold stars
    x_gal = x[distances < rgalStar1] # Filter out all stars that are too far away 
    y_gal = y[distances < rgalStar1] # Filter out all stars that are too far away 
    z_gal = z[distances < rgalStar1] # Filter out all stars that are too far away
    xMean = numpy.mean(x)
    yMean = numpy.mean(y)
    zMean = numpy.mean(z)
    print(x)
    print(xMean)
    x = x - xMean
    print(x)
    y = y - yMean
    z = z - zMean
    
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

# TODO - plt.colorbar - a legend of colors; have any age above .5 be a set color; have different colors for different age ranges
# Graph 1, Graph 2 - Two graphs of all stars in a galaxy, having color represent age
star_counts = cu.getField(data, 'n_star')
print(star_counts)
x,y,z,a,m,v = getGalaxyStarInfo(2)
# x2,y2,z2,a2,m2,v2 = getGalaxyStarInfo(2)
x2,y2,z2,a2,m2,v2 = getGalaxyStarInfo2(2)
cu.star_scatter_plot(x,y,z,a, 1)
cu.star_scatter_plot(x2,y2,z2,a2,1)



# Make a histogram instead of a py chart
# Graph 3 - pie chart show all galaxies that have at least 3 stars
labels = ['Halo {} - {}'.format(i, count) for i, count in enumerate(star_counts) if count > 3]
# cu.pie_chart(star_counts, labels) 


# Graph 5 - M_star vs M_vir (mass of all stars vs mass of halo?)
# Log slope this and find the slope of the linear line that is produced
# start putting units on all graph axis
M_vir = cu.getField(data, "M_vir")
M_star = cu.getField(data, "M_star")
cu.graph(x = M_vir, y = M_star, title = "M_star vs M_vir", showLine=False)
M_vir = numpy.log10(M_vir)
M_star = numpy.log10(M_star)
cu.graph(x = M_vir, y = M_star, title = "M_star vs M_vir", showLine=False)



# Make a histogram
# import scipy.stats.bin from discord
    # use sum and mean
    # for dark matter: x pos, y pos, and mas - use sum
    # for stars, x, y, and a - use sum
    # look at plt.imshow
    # stars - x,y, vz - use standard deviatian

    # pay attention to scaling on why the axis are not centered a 0,0