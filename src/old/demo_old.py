import astropy.io.ascii as ascii
import gizmo_analysis as gizmo

# Original way to do this
particles = gizmo.io.Read.read_snapshots(
        simulation_directory = '../data',
        snapshot_directory = '../data',
        species=['star'], # Take note
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
data = ascii.read("../data/snapshot_600.z0.000.AHF_halos")
data_filtered = data[(data.field('fMhires(38)') > 0.99)]
data = data_filtered

# Print hubble constant
print(h) 
# Print all particles you specified in the line with "Take Note"
print(particles.keys())
# Print properties that can be accessed for stars 
print(particles['star'].keys()) 
# Print the x pos of every star in the simulation
print(particles['star']['position'][:,0])
# Print the number of stars in each halo
print(data.field("n_star(64)"))