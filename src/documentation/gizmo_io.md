# Table of Contents
- ParticleDictionaryClass Class
    - [_assign_ids_to_indices](#_assign_ids_to_indices): Assign to self an array [and dictionary] to point from a particle's id [and id.child]
to its array index in this dictionary catalog.
    - [_get_abundances](#_get_abundances): Get element mass fraction[s] or metallicity[s], either stored in dictionary or via
post-processing stored age-tracer mass weights.
    - [get_pointers_from_ids](#get_pointers_from_ids): For each input id [and child id], get a pointer to the array index [and species name] of
the particle in this dictionary catalog.
If running from within dictionary of single particle species (such as part['star']),
return only the pointer index.
If running from within the meta-dictionary of multiple species (such as part),
return the pointer index and species name.
    - [init ParticleDictionaryClass](#init-particledictionaryclass)
    - [prop](#prop): Get property, either stored in self's dictionary or derive it from stored properties.
Can compute basic mathematical manipulations/combinations, for example:
    'log temperature', 'temperature / density', 'abs position'
- ReadClass Class
    - [_adjust_particle_properties](#_adjust_particle_properties): Adjust properties for each species, including unit conversions, separating dark species by
mass, sorting by id, and subsampling.
    - [_assign_hosts_coordinates_from_halos](#_assign_hosts_coordinates_from_halos): Utility function for assign_hosts_coordinates().
Read and assign host halo positions and velocities from the halo catalog at that snapshot.
    - [_assign_hosts_coordinates_from_particles](#_assign_hosts_coordinates_from_particles): Utility function for assign_hosts_coordinates().
Compute and assign host galaxy positions and velocities from the particles.
    - [_check_particle_properties](#_check_particle_properties): Checks sanity of particle properties, print warning if they are outside given limits.
    - [_get_check_value](#_get_check_value)
    - [_get_cosmology](#_get_cosmology): Return cosmological parameters via Cosmology dictionary class.
If input simulation_directory, (try to) read cosmological parameters from the MUSIC
initial condition config file within this simulation_directory.
If cannot find the MUSIC config file, use what cosmological parameters are input,
and assume the rest from the AGORA simulation.
    - [_read_particles](#_read_particles): Read from snapshot file[s] all particles of species type[s] in self._species_read list.
    - [assign_hosts_coordinates](#assign_hosts_coordinates): Assign position [kpc comoving] and velocity [km / s] of each host galaxy/halo.
Use species_name, if defined, else default to stars for baryonic simulation or
dark matter for dark matter-only simulation.
    - [assign_hosts_rotation](#assign_hosts_rotation): Compute and assign rotation tensor and ratios of principal axes
(defined via the moment of inertia tensor) for each host galaxy.
By default, use stars for baryonic simulations, or if no stars in catalog, use gas.
    - [assign_particle_orbits](#assign_particle_orbits): Assign derived orbital properties wrt each host to each particle species.
    - [get_snapshot_file_names_indices](#get_snapshot_file_names_indices): Get name of file or directory (with relative path) and index for all snapshots in directory.
If input valid snapshot_index, get its file name (if multiple files per snapshot, get name
of 0th one).
If input snapshot_index as None or 'all', get name of file/directory and index for each
snapshot file/directory.
    - [init ReadClass](#init-readclass): Set properties for snapshot files.
    - [read_header](#read_header): Read header from snapshot file.
    - [read_snapshots](#read_snapshots): Read given properties for given particle species from simulation snapshot file[s].
Can read single snapshot or multiple snapshots.
If single snapshot, return as dictionary class;
if multiple snapshots, return as list of dictionary classes.
    - [read_snapshots_simulations](#read_snapshots_simulations): Read snapshots at the same redshift from different simulations.
Return as list of dictionaries.
- WriteClass Class
    - [rewrite_snapshot](#rewrite_snapshot): Read single snapshot.
Rewrite, deleting given species and/or adjusting particle properties.
    - [rewrite_snapshot_to_text](#rewrite_snapshot_to_text): Re-write snapshot to text file, one file per particle species.
    - [write_exsitu_flag](#write_exsitu_flag): Read single snapshot, with star coordinates at formation.
Apply a total distance treshold to define ex-situ stars.
Write text file that contains binary flag of whether star particle formed ex-situ.

## _assign_ids_to_indices

### Description:
Assign to self an array [and dictionary] to point from a particle's id [and id.child]
to its array index in this dictionary catalog.


## _get_abundances

### Description:
Get element mass fraction[s] or metallicity[s], either stored in dictionary or via
post-processing stored age-tracer mass weights.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| property_name | str |     name of property to get | required |
| indices | array [optional] |     indices of particles to get property of | required |

## get_pointers_from_ids

### Description:
For each input id [and child id], get a pointer to the array index [and species name] of
the particle in this dictionary catalog.
If running from within dictionary of single particle species (such as part['star']),
return only the pointer index.
If running from within the meta-dictionary of multiple species (such as part),
return the pointer index and species name.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| ids | array |     ids of particles | required |
| child_ids | array |     child ids of particles | required |

## init ParticleDictionaryClass

## prop

### Description:
Get property, either stored in self's dictionary or derive it from stored properties.
Can compute basic mathematical manipulations/combinations, for example:
    'log temperature', 'temperature / density', 'abs position'


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| property_name | str |     name of property | required |
| indices | array |     indices of particles to get properties of | required |
| _dict_only | bool |     require property_name to be in self's dict - avoids endless recursion<br />    primarily for internal/recursive usage of this function | required |

## _adjust_particle_properties

### Description:
Adjust properties for each species, including unit conversions, separating dark species by
mass, sorting by id, and subsampling.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dictionary class |     catalog of particles at snapshot | required |
| header | dict |     header dictionary | required |
| particle_subsample_factor | int |     factor to periodically subsample particles, to save memory | required |
| separate_dark_lowres | bool |     whether to separate low-resolution dark matter into separate dicts according to mass | required |
| sort_dark_by_id | bool |     whether to sort dark-matter particles by id | required |

## _assign_hosts_coordinates_from_halos

### Description:
Utility function for assign_hosts_coordinates().
Read and assign host halo positions and velocities from the halo catalog at that snapshot.


## _assign_hosts_coordinates_from_particles

### Description:
Utility function for assign_hosts_coordinates().
Compute and assign host galaxy positions and velocities from the particles.


## _check_particle_properties

### Description:
Checks sanity of particle properties, print warning if they are outside given limits.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dictionary class |     catalog of particles | required |

## _get_check_value

## _get_cosmology

### Description:
Return cosmological parameters via Cosmology dictionary class.
If input simulation_directory, (try to) read cosmological parameters from the MUSIC
initial condition config file within this simulation_directory.
If cannot find the MUSIC config file, use what cosmological parameters are input,
and assume the rest from the AGORA simulation.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| simulation_directory | str |     base directory of simulation | required |

## _read_particles

### Description:
Read from snapshot file[s] all particles of species type[s] in self._species_read list.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| simulation_directory | str |     directory of simulation | required |
| snapshot_directory | str |     directory of snapshot files within simulation_directory | required |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift' | required |
| snapshot_value | int or float |     input snapshot number kind: 'index', 'redshift' | required |
| properties | str or list |     name[s] of particle properties to read. 'all' or None = read all properties in snapshot | required |
| elements | str or list |     name[s] of elemental abundances to read. 'all' or None = read all elements in snapshot | required |
| convert_float32 | bool |     whether to convert all floats to 32 bit to save memory | required |
| header | dict |     snapshot file header information (from previous read of file) | required |

## assign_hosts_coordinates

### Description:
Assign position [kpc comoving] and velocity [km / s] of each host galaxy/halo.
Use species_name, if defined, else default to stars for baryonic simulation or
dark matter for dark matter-only simulation.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dictionary class |     catalog of particles at snapshot | required |
| method | str |     method to use to get host coordinates.<br />    if a string, tells the code which method to use:<br />        'track' : reads host coordinates from track/host_coordinates.hdf5,<br />            compiled during particle tracking using only stars in each host at z = 0<br />        'halo' : reads host halo coordinates from halo/rockstar_dm/catalog_hdf5/<br />        'mass' or 'potential' or 'massfraction.metals' : assign coordinates during read in<br />            via iterative zoom-in, weighting each particle by that property<br />    if True (default), will try a few methods in the following order of preference:<br />        if a baryonic simulation, try 'track' then 'mass'<br />        if a DM-only simulations, try 'halo' then 'mass' | required |
| species_name | str |     which particle species to use to define center<br />    relevant only if method is 'mass' or 'potential' or 'massfraction.metals' | required |
| part_indicess | array or list of arrays |     list of indices of particles to use to define host center coordinates<br />    if supply a list of arrays, use each list element for a different host | required |
| assign_formation_coordinates | bool |     whether to assign to stars their coordindates wrt each host galaxy at formation<br />    (if reading hosts coordinates from file) | required |
| velocity_distance_max | float |     maximum distance to keep particles to compute velocity | required |
| host_number | int |     number of hosts to assign | required |
| exclusion_distance | float |     radius around previous hosts' center position[s] to exclude particles in<br />    finding center of next host [kpc comoving] | required |
| simulation_directory | str |     directory of simulation | required |
| snapshot_directory | str |     directory of snapshot files, within simulation_directory | required |
| track_directory | str |     directory of files for particle pointers, formation coordinates, and host coordinates | required |

## assign_hosts_rotation

### Description:
Compute and assign rotation tensor and ratios of principal axes
(defined via the moment of inertia tensor) for each host galaxy.
By default, use stars for baryonic simulations, or if no stars in catalog, use gas.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dictionary class |     catalog of particles at snapshot | required |
| species_name | string |     name of particle species to use to determine rotation | required |
| distance_max | float |     maximum distance to select particles [kpc physical] | required |
| mass_percent | float |     keep particles within the distance that encloses mass percent [0, 100] of all particles<br />    within distance_max | required |
| age_percent | float |     keep youngest age_percent of (star) particles within distance cut | required |

## assign_particle_orbits

### Description:
Assign derived orbital properties wrt each host to each particle species.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dictionary class |     catalog of particles at snapshot | required |
| species | str or list |     particle species to compute | required |
| host_positions | array or array of arrays |     position[s] of hosts | required |
| host_velocities | array or array of arrays |     velocity[s] of hosts | required |

## get_snapshot_file_names_indices

### Description:
Get name of file or directory (with relative path) and index for all snapshots in directory.
If input valid snapshot_index, get its file name (if multiple files per snapshot, get name
of 0th one).
If input snapshot_index as None or 'all', get name of file/directory and index for each
snapshot file/directory.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| directory | str |     directory to check for files | required |
| snapshot_index | int |     index of snapshot: if None or 'all', get all snapshots in directory | required |
| snapshot_block_index | int |     index of file block (if multiple files per snapshot)<br />    if None or 'all', return names of all file blocks for snapshot | required |

## init ReadClass

### Description:
Set properties for snapshot files.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| verbose | bool |     whether to print diagnostics | required |

## read_header

### Description:
Read header from snapshot file.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift' | required |
| snapshot_value | int or float |     input snapshot number kind: 'index', 'redshift' | required |
| simulation_directory | str |     directory of simulation | required |
| snapshot_directory | str |     directory of snapshot files within simulation_directory | required |
| simulation_name | str |     name to store for future identification | required |
| snapshot_block_index | int |     index of file block (if multiple files per snapshot) | required |
| verbose | bool |     whether to print number of particles in snapshot | required |

## read_snapshots

### Description:
Read given properties for given particle species from simulation snapshot file[s].
Can read single snapshot or multiple snapshots.
If single snapshot, return as dictionary class;
if multiple snapshots, return as list of dictionary classes.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| species | str or list |     name[s] of particle species:<br />        'all' = all species in file<br />        'dark' = dark matter at highest resolution<br />        'dark2' = dark matter at lower resolution<br />        'gas' = gas<br />        'star' = stars<br />        'blackhole' = black holes, if snapshot contains them | required |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift', 'scalefactor' | required |
| snapshot_values | int or float or list |     index[s] or redshift[s] or scale-factor[s] of snapshot[s] | required |
| simulation_directory | str |     directory of simulation | required |
| snapshot_directory | str |     directory of snapshot files within simulation_directory | required |
| track_directory | str |     directory of files for particle pointers, formation coordinates, and host coordinates | required |
| simulation_name | str |     name to store for future identification | required |
| properties | str or list |     name[s] of particle properties to read. 'all' or None = read all properties in snapshot | required |
| elements | str or list |     name[s] of elemental abundances to read. 'all' or None = read all elements in snapshot | required |
| particle_subsample_factor | int |     factor to periodically subsample particles, to save memory | required |
| separate_dark_lowres | bool |     whether to separate low-resolution dark matter into separate dicts according to mass | required |
| sort_dark_by_id | bool |     whether to sort dark-matter particles by id | required |
| convert_float32 | bool |     whether to convert all floats to 32 bit to save memory | required |
| host_number | int |     number of hosts to assign and compute coordinates relative to | required |
| assign_hosts | bool or str |     whether to assign coordinates of each host.<br />    if a string, tells the code which method to use:<br />        'track' : reads host coordinates from track/host_coordinates.hdf5,<br />            compiled during particle tracking using only stars in each host at z = 0<br />        'halo' : reads host halo coordinates from halo/rockstar_dm/catalog_hdf5/<br />        'mass' or 'potential' or 'massfraction.metals': assign coordinates during read in<br />            via iterative zoom-in, weighting each particle by that property<br />    if True (default), will try a few methods in the following order of preference:<br />        if a baryonic simulation (or input species_name='star'), try 'track' then 'mass'<br />        if a DM-only simulations (or input species_name='dark'), try 'halo' then 'mass' | required |
| assign_hosts_rotation | bool |     whether to assign principal axes rotation tensor of each host galaxy | required |
| assign_orbits | bool |     whether to assign orbital properties wrt each host galaxy/halo | required |
| assign_formation_coordinates | bool |     whether to assign to stars their coordindates wrt each host galaxy at formation | required |
| assign_pointers | bool |     whether to assign pointers for tracking particles from z = 0 to this snapshot | required |
| check_properties | bool |     whether to check sanity of particle properties after read in | required |

## read_snapshots_simulations

### Description:
Read snapshots at the same redshift from different simulations.
Return as list of dictionaries.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| species | str or list |     name[s] of particle species to read | required |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift', 'scalefactor' | required |
| snapshot_value | int or float |     input snapshot number kind: 'index', 'redshift', 'scalefactor' | required |
| simulation_directories | list or dict |     list of simulation directories, or dict of simulation_directories: simulation_names | required |
| snapshot_directory | str |     directory of snapshot files within simulation_directory | required |
| track_directory | str |     directory of files for particle pointers, formation coordinates, and host coordinates | required |
| properties | str or list |     name[s] of particle properties to read. 'all' or None = read all properties in snapshot | required |
| elements | str or list |     name[s] of elemental abundances to read. 'all' or None = read all elements in snapshot | required |
| assign_hosts | bool or str |     whether to assign host coordinates.<br />    if a string, tells the code which method to use:<br />        'track' : reads host coordinates from track/host_coordinates.hdf5,<br />            compiled during particle tracking using only stars in each host at z = 0<br />        'halo' : reads host halo coordinates from halo/rockstar_dm/catalog_hdf5/<br />        'mass' or 'potential' : assign coordinates during read in via iterative zoom-in,<br />            weighting each particle by that property<br />    if True (default), will try a few methods in the following order of preference:<br />        if a baryonic simulation (or input species_name='star'), try 'track' then 'mass'<br />        if a DM-only simulations (or input species_name='dark'), try 'halo' then 'mass' | required |
| assign_hosts_rotation | bool |     whether to assign principal axes rotation tensor of each host galaxy/halo | required |
| assign_orbits | bool |     whether to assign orbital properties wrt each host galaxy/halo | required |
| assign_formation_coordinates | bool |     whether to assign to stars their coordindates wrt each host galaxy at formation | required |
| assign_pointers | bool |     whether to assign pointers for tracking particles from z = 0 to this snapshot | required |
| check_properties | bool |     whether to check sanity of particle properties after read in | required |

## rewrite_snapshot

### Description:
Read single snapshot.
Rewrite, deleting given species and/or adjusting particle properties.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| species | str or list |     name[s] of particle species to delete:<br />    'gas' = gas<br />    'dark' = dark matter at highest resolution<br />    'dark2' = dark matter at lower resolution<br />    'star' = stars<br />    'blackhole' = black holes | required |
| action | str |     what to do to snapshot file: 'delete', 'velocity' | required |
| value_adjust | float |     value by which to adjust property (if not deleting) | required |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift' | required |
| snapshot_value | int or float |     input snapshot number kind: 'index', 'redshift' | required |
| simulation_directory | str |     directory of simulation | required |
| snapshot_directory | str |     directory of snapshot files within simulation_directory | required |

## rewrite_snapshot_to_text

### Description:
Re-write snapshot to text file, one file per particle species.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| part | dict class |     catalog of particles at snapshot | required |

## write_exsitu_flag

### Description:
Read single snapshot, with star coordinates at formation.
Apply a total distance treshold to define ex-situ stars.
Write text file that contains binary flag of whether star particle formed ex-situ.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| snapshot_value_kind | str |     input snapshot number kind: 'index', 'redshift' | required |
| snapshot_value | int or float |     input snapshot number kind: 'index', 'redshift' | required |
| simulation_directory | str |     directory of simulation | required |
| track_directory | str |     directory of files for particle pointers, formation coordinates, and host coordinates | required |
| exsitu_distance | float |     minimum distance to define ex-situ stars [kpc physical or comoving] | required |
| exsitu_distance_scaling | bool |     whether to scale exsitu_distance with scale-factor at formation | required |
