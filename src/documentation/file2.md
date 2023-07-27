# Table of Contents
 - [dist](#dist): Get the magnitude of a vector with three components (x, y, z).
 - [filter_list](#filter_list): Filter a list, comparing each corresponding distance to rgal.
If distances[i] < rgal, then the item will be added to the returned list.

## dist

### Description:
Get the magnitude of a vector with three components (x, y, z).


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| x | float, int  |         The x component of the vector. | required |
| y | float, int  |         The y component of the vector. | required |
| z | float, int  |         The z component of the vector. | required |

## filter_list

### Description:
Filter a list, comparing each corresponding distance to rgal.
If distances[i] < rgal, then the item will be added to the returned list.


### Parameters:


| Name | Type | Description | Default |
| --- | --- | --- | --- |
| list | numpy array  |         The list to filter. | required |
| distances | numpy array  |         The comparison list. | required |
| rgal | float, int  |         The constant to compare. | required |
