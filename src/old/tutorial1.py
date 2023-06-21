from stellarutil import cross as cu

# Get the data tables and constant we will be working with
particles = cu.getParticles()
h = cu.getH() 
data = cu.get_ahf_data("snapshot_600.z0.000.AHF_halos") # filter out anything not in that 1% - gotta ask what is happening here again

# Get the center of the biggest galaxy in the simulation
xc = data.field('Xc(6)')[0] / h
yc = data.field('Yc(7)')[0] / h
zc = data.field('Zc(8)')[0] / h

# Get the x,y,z positions of each star particle in the simulation, and normalize it with the galaxy center
x = particles['star']['position'][:,0] - xc
y = particles['star']['position'][:,1] - yc
z = particles['star']['position'][:,2] - zc

# Get the scalefactor (age) of the star
a = particles['star']['form.scalefactor']


# Get the stars in the galaxy
distances = cu.dist(x,y,z) # Get the distance of each particle from the center galaxy
rgal = 0.15 * data.field('Rvir(12)')[0] / h # Get the radius of the galaxy that can actually hold stars
x_gal = x[distances < rgal] # Filter out all stars that are too far away 
y_gal = y[distances < rgal] # Filter out all stars that are too far away 
z_gal = z[distances < rgal] # Filter out all stars that are too far away 
a_gal = a[distances < rgal] # Filter out all stars that are too far away

# Make the scatter plot
# cu.scatter_plot(x_gal,y_gal,z_gal)

print(particles['star'].keys())
print(particles.keys())
print(particles['gas']['mass']) # big list of just 261.751 repeated??? Is this a special case
print("The list of star positions:")
print(particles['star']['position'])
print("The center of the galaxy: on the x axis:")
print(xc)
print("The position of the first star:")
print(particles['star']['position'][0])
print("The position of the second star:")
print(particles['star']['position'][1])
print("##########")
print("The x position of the first star :")
print(x[0])
print("The x position of the first star (raw):")
print(particles['star']['position'][:,0][0])

