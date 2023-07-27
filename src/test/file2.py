import numpy as np

def dist(x, y, z):
    '''
    Get the magnitude of a vector with three components (x, y, z).

    Parameters:
    ----------
        x : float | int 
            The x component of the vector.
        y : float | int 
            The y component of the vector.
        z : float | int 
            The z component of the vector.

    Returns
    -------
        The magnitude of the vector.
    '''

    result = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    return result



def filter_list(list, distances, rgal):
    '''
    Filter a list, comparing each corresponding distance to rgal.
    If distances[i] < rgal, then the item will be added to the returned list.

    Parameters:
    ----------
        list : numpy array 
            The list to filter.
        distances : numpy array 
            The comparison list.
        rgal : float | int 
            The constant to compare.

    Returns
    -------
        The magnitude of the vector.
    '''
    if isinstance(list, np.ndarray):
        return list[distances < rgal]
    else:
        raise TypeError("Unsupported input type. Only NumPy arrays are accepted.")

