from stellarutil import cross as cu

# Get the data tables and constant we will be working with
particles = cu.getParticles()
h = cu.getH() 
data = cu.get_ahf_data("snapshot_600.z0.000.AHF_halos") # filter out anything not in that 1% - gotta ask what is happening here again


# 1, 2, 4, 6-13, 17, 18, 44, 45, 64, 65
# Make a nice print out that displays all the data in this list
fields = [1, 2, 4] + list(range(6, 14)) + [17, 18, 44, 45, 64, 65]
for field in fields:
    list = cu.getField(data, field)[:5]
    print(f"{cu.getFieldName(data, field)}:")
    print(f"[{ ' '.join(map(str, list)) }]\n")

# see how close rmax is to 15% of the radius
rgal = 0.15 * data.field('Rvir(12)')[0] / h
rmax = cu.getField(data, 'Rmax')[0] / h
print(f"15% of the radius of the biggest galaxy is: {rgal}")
print(f"The max radius of the biggest galaxy is: {rmax}")
print(f"The % difference between the two values is: {abs((rgal - rmax) / ((rgal + rmax) / 2)) * 100}%")