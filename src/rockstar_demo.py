import halo_analysis as halo, utilities as ut, numpy as np


simulation_directory = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/cr_heating_fix/m10v_r30'

# read a halo catalog at single snapshot (z = 0)
# read halo catalogs at all available snapshots by supplying None or 'all' as the input snapshot value.
# read_catalogs() returns this as a list of dictionaries, with the list index being the snapshot index.
# beware - this can take a while to read...
hal = halo.io.IO.read_catalogs('index', 600, simulation_directory)

# hal is a dictionary of halo properties
for k in hal.keys():
    print(k)

print('--------------------------------------')
# 3-D position (particle_number x dimension_number array) [kpc comoving]
hal['position']
# 3-D velocity (particle_number x dimension_number array) [km/s physical]
hal['velocity']
# DM mass of halo [M_sun]
# by default, I run the halo finder using 200m (200 x the mean matter density) for the default overdensity/virial definition
hal['mass']
# but Rockstar also stores halo DM mass based on different virial definitions [M_sun]
print('{}\n{}\n{}'.format(hal['mass.200m'], hal['mass.vir'], hal['mass.200c']))
# DM mass that is bound to halo [M_sun]
hal['mass.bound']
# halo radius [kpc physical] again using 200m for the overdensity/virial definition 
hal['radius']
# NFW scale radius [kpc physical]
hal['scale.radius']
# maximum of the circular velocity profile [km/s physical]
hal['vel.circ.max']
# standard deviation of the velocity (velocity dispersion) [km/s physical]
hal['vel.std']
# the fraction of DM mass within R_200m that is low-resolution DM
# this is a derived quantity, so you need to call via the .prop() function (see below)
hal.prop('lowres.mass.frac')
# index of the primary (most massive) host halo in the catalog
hal['host.index']
# 3-D distance to the primary host halo [kpc physical]
hal['host.distance']
# total (scalar) distance to the primary host halo [kpc physical]
# this is a derived quantity, so you need to call via the .prop() function (see below)
hal.prop('host.distance.total')
# 3-D velocity wrt the primary host halo [kpc physical]
# radial and tangential velocity wrt the primary host halo [kpc physical]
print(hal['host.velocity'])
print(hal['host.velocity.rad'])
print(hal['host.velocity.tan'])
# use .prop() to compute derived quantities
# this can handle simple arithmetic conversions, such as the log of the mass, or the ratio of masses
# see halo.io.HaloDictionaryClass for all options for derived quantities
print(hal.prop('host.distance.total'))
print(hal.prop('host.velocity.total'))
print(hal.prop('log mass'))
print(hal.prop('mass.bound / mass'))