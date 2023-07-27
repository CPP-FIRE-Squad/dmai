# Table of Contents
 - [getGalaxyStarInfo](#getgalaxystarinfo)
 - [get_ahf_data](#get_ahf_data): Read the content from the AHF file.
 - [get_field](#get_field): Return a list of all items in the field of a dataset
 - [get_field_name](#get_field_name): Return the correct name of a field.
 - [get_hubble_constant](#get_hubble_constant)
 - [get_particles](#get_particles)
- Simulation Class
    - [get_field](#get_field): Get the values in the column of the specified field from the .AHF_halos file.
    - [get_stars_in_halo](#get_stars_in_halo): Get the list of stars inside an indicated dark matter halo.
    - [init Simulation](#init-simulation): Initialize a new Simulation object.
- Star Class
    - [get_velocity](#get_velocity): Get the velocity of the star by calculating the magnitude of the velocity vector.
    - [init Star](#init-star): Initialize a new Star object.

## getGalaxyStarInfo

## get_ahf_data

### Description:
Read the content from the AHF file.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| path | str  |         the path to the file | required |
| filter | boolean |         should the data be filtered to only include rows with fMhires(38) > 0.99. | True |

## get_field

### Description:
Return a list of all items in the field of a dataset


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| data | 2d array, int, string |         The table of data  | [] |
| field | int, string  |         The name of the field. | required |

## get_field_name

### Description:
Return the correct name of a field.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| data | 2d array |         The table of data  | required |
| name | int, string  |         The name of the field. | required |

## get_hubble_constant

## get_particles

## get_field

### Description:
Get the values in the column of the specified field from the .AHF_halos file.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| field | string |     The name of the field. | required |

## get_stars_in_halo

### Description:
Get the list of stars inside an indicated dark matter halo.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| index | int |     The index of the dark matter halo. | 0 |
| restrict | float |     The restriction percentage. | 15 |

## init Simulation

### Description:
Initialize a new Simulation object.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| simulation_directory | string |     The path to the .hdf5 file. | '../data' |
| snapshot_directory | string |     The path to the snapshot_times.txt.  | '../data' |
| ahf_directory | string |     The path to the .AHF_halos file.  | '../data' |
| species | list |     name[s] of particle species:<br />        'all' = all species in file<br />        'dark' = dark matter at highest resolution<br />        'dark2' = dark matter at lower resolution<br />        'gas' = gas<br />        'star' = stars<br />        'blackhole' = black holes, if snapshot contains them | required |
| snapshot_values | int or float or list |     index[s] or redshift[s] or scale-factor[s] of snapshot[s] | required |

## get_velocity

### Description:
Get the velocity of the star by calculating the magnitude of the velocity vector.


## init Star

### Description:
Initialize a new Star object.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| x | float |     The x position of the star. | required |
| y | float |     The y position of the star. | required |
| z | float |     The z position of the star. | required |
| m | float |     The mass of the star. | required |
| a | float |     The scale factor of the star. | required |
| vx | float |     The x velocity of the star. | required |
| vy | float |     The y velocity of the star. | required |
| vz | float |     The z velocity of the star. | required |
