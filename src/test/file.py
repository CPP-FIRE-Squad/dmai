import astropy.io.ascii as ascii
import gizmo_analysis as gizmo
from stellarutil.calculations import dist, filter_list

#region talk to gizmo stuff


# Hello!!


def get_ahf_data(path, filter = True):
    '''
    Read the content from the AHF file.
    
    Parameters:
    ----------
        path : str 
            the path to the file
        filter : boolean = True
            should the data be filtered to only include rows with fMhires(38) > 0.99.

    Returns
    -------
        Return a 2D list of the data specified in the AHF file.
    '''

    data = ascii.read(path)
    if not filter:
        return data
    else:
        data_filtered = data[(data.field('fMhires(38)') > 0.99)]
        return data_filtered

def get_field(data, field):
    """
    Return a list of all items in the field of a dataset

    Parameters:
    ----------
        data : 2d array | int | string = []
            The table of data 
        field : int | string 
            The name of the field.

    Returns
    -------
        The ist of all items in the field.
    """
    
    # Get the correct name of the field
    field_name = get_field_name(data, field)
    print(f"{field} changed to {field_name}")
    # Store all the field data in a list called column
    column = data.field(field_name) 
    # Return the column
    return column

def get_field_name(data, name):
    '''
    Return the correct name of a field.

    Parameters:
    ----------
        data : 2d array
            The table of data 
        name : int | string 
            The name of the field.

    Returns
    -------
        The ist of all items in the field.
    '''
    name = str(name).lower()  # Convert field to string if it's an integer
    name = name.replace('_','')
    # Loop through all the field names
    for item in data.dtype.names:
        string = item.lower().replace('_','')
        if name in string:
            return item
            
    return None

def get_hubble_constant(simulation_directory, snapshot_directory, snapshot_value):
    '''
    '''
    header = gizmo.io.Read.read_header(
        simulation_directory = simulation_directory,
        snapshot_directory = snapshot_directory,
        snapshot_value_kind = 'index',
        snapshot_value = snapshot_value
    )
    return header['hubble']

def get_particles(simulation_directory, snapshot_directory, species, snapshot_values):
    return gizmo.io.Read.read_snapshots(
        simulation_directory = simulation_directory,
        snapshot_directory = snapshot_directory,
        species=species, 
        snapshot_value_kind='index',
        snapshot_values=snapshot_values
    )
#endregion


def getGalaxyStarInfo(data, particles, h, index = 0):
    # Get the center of the halo
    xc = data.field('Xc(6)')[index] / h
    yc = data.field('Yc(7)')[index] / h
    zc = data.field('Zc(8)')[index] / h
    # Get the peculiar velocity of halo
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
    from calculations import dist
    distances =  dist(x,y,z) # Get the distance of each particle from the center galaxy
    rgal = 0.15 * data.field('Rvir(12)')[index] / h # Get the radius of the galaxy that can actually hold stars
    x_gal = x[distances < rgal] # Filter out all stars that are too far away 
    y_gal = y[distances < rgal] # Filter out all stars that are too far away 
    z_gal = z[distances < rgal] # Filter out all stars that are too far away 
    a_gal = a[distances < rgal] # Filter out all stars that are too far away
    m_gal = m[distances < rgal] # Filter out all stars that are too far away
    vx_gal = vx[distances < rgal] # Filter out all stars that are too far away
    vy_gal = vy[distances < rgal] # Filter out all stars that are too far away
    vz_gal = vz[distances < rgal] # Filter out all stars that are too far away
    v_gal = dist(vx_gal, vy_gal, vz_gal)

    return (x_gal, y_gal, z_gal, a_gal, m_gal, v_gal)


class Star:

    def __init__(self, x = 0, y = 0, z = 0, m = 0, a = 0, vx = 0, vy = 0, vz = 0):
        """
        Initialize a new Star object.

        Parameters:
        ----------
        x : float
            The x position of the star.
        y : float
            The y position of the star.
        z : float
            The z position of the star.
        m : float
            The mass of the star.
        a : float
            The scale factor of the star.
        vx : float
            The x velocity of the star.
        vy : float
            The y velocity of the star.
        vz : float
            The z velocity of the star.

        Attributes:
        -----------
        x : float
            The x position of the star.
        y : float
            The y position of the star.
        z : float
            The z position of the star.
        m : float
            The mass of the star.
        a : float
            The scale factor of the star.
        vx : float
            The x velocity of the star.
        vy : float
            The y velocity of the star.
        vz : float
            The z velocity of the star.
        velocity : float
            The velocity of the star.
        """
        self.x = x
        self.y = y
        self.z = z
        self.m = m
        self.a = a
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.velocity = self.get_velocity()
    
    def get_velocity(self):
        """
        Get the velocity of the star by calculating the magnitude of the velocity vector.
        
        Returns
        -------
        The velocity of the star.
        """
        return dist(self.vx, self.vy, self.vz)
    
    def __str__(self):
        """
        The toString method for converting the object to a string.
        
        Returns
        -------
        A stringified version of the object.
        """
        output = f"Star:\n  Position: ({self.x}, {self.y}, {self.z}) [kpc]\n  Mass: {self.m} [unit]\n  Scale Factor (a): {self.a} [unit]\n  Velocity: {self.velocity()} [kpc/s]"
        return output

class Simulation:

    def __init__(self, simulation_directory = '../data', snapshot_directory = '../data', ahf_directory = "../data/snapshot_600.z0.000.AHF_halos", species = ['star'], snapshot_values = 600, snapshot_value = 600):
        """
        Initialize a new Simulation object.

        Parameters:
        ----------
        simulation_directory : string = '../data'
            The path to the .hdf5 file.
        snapshot_directory : string = '../data'
            The path to the snapshot_times.txt. 
        ahf_directory : string = '../data'
            The path to the .AHF_halos file. 
        species : list
            name[s] of particle species:
                'all' = all species in file
                'dark' = dark matter at highest resolution
                'dark2' = dark matter at lower resolution
                'gas' = gas
                'star' = stars
                'blackhole' = black holes, if snapshot contains them
        snapshot_values : int or float or list
            index[s] or redshift[s] or scale-factor[s] of snapshot[s]

        Attributes:
        -----------
        h : float
            The hubble constant.
        particles : float
            The data for all the indicated particles in the simulation.
        ahf_data : float
            The data within the .AHF_halos file. 
        """

        self.h = get_hubble_constant(simulation_directory, snapshot_directory, snapshot_value)
        self.particles = get_particles(simulation_directory, snapshot_directory, species, snapshot_values)
        self.ahf_data = get_ahf_data(ahf_directory)


    def get_stars_in_halo(self, index = 0, percentage = 15):
        """
        Get the list of stars inside an indicated dark matter halo.

        Parameters:
        ----------
        index : int = 0
            The index of the dark matter halo.
        restrict : float = 15
            The restriction percentage.
        
        Returns
        -------
        The list of all stars in the indicated dark matter halo.
        """
        # Get the center of the indicated dark matter halo
        xc = self.ahf_data.field('Xc(6)')[index] / self.h
        yc = self.ahf_data.field('Yc(7)')[index] / self.h
        zc = self.ahf_data.field('Zc(8)')[index] / self.h
        # Get the peculiar velocity of the indicated dark matter halo
        vxc = self.ahf_data.field('VXc(9)')[index] / self.h
        vyc = self.ahf_data.field('VYc(10)')[index] / self.h
        vzc = self.ahf_data.field('VZc(11)')[index] / self.h
        # Get the x,y,z positions of each star particle in the simulation
        # And normalize it with the center of the indicated dark matter halo
        x = self.particles['star']['position'][:,0] - xc
        y = self.particles['star']['position'][:,1] - yc
        z = self.particles['star']['position'][:,2] - zc
        # Get the scalefactor (age) of each star in the simulation
        a = self.particles['star']['form.scalefactor']
        # Get the mass of each star in the simulation
        m = self.particles['star']['mass']
        # Get the x,y,z velocity of each star particle in the simulation
        # And normalize it with the peculiar velocity of the indicated dark matter halo
        vx = self.particles['star']['velocity'][:,0] - vxc
        vy = self.particles['star']['velocity'][:,1] - vyc
        vz = self.particles['star']['velocity'][:,2] - vzc
        # Check AHF file for the number of stars in the the indicated dark matter halo
        num_stars = self.ahf_data.field('n_star(64)')[index]
        print(f"This dark matter halo has {num_stars} star(s) according to the AHF file.")
        # Ensure the percentage has proper bounds
        if percentage < 0:
            percentage = 0
        elif percentage > 100:
            percentage = 100
        # Loop through and capture the stars 
        while percentage <= 100:
            # Get the distance of each star from the center of the indicated dark matter halo
            distances =  dist(x,y,z) 
            # Get the radius of the galaxy that can actually hold stars
            rgal = (percentage / 100.0) * self.ahf_data.field('Rvir(12)')[index] / self.h 
            # Filter out all stars that are too far away 
            print(f"Filtering at: {percentage}%")
            x_gal = filter_list(x, distances, rgal)
            y_gal = filter_list(y, distances, rgal)
            z_gal = filter_list(z, distances, rgal)
            a_gal = filter_list(a, distances, rgal)
            m_gal = filter_list(m, distances, rgal)
            vx_gal = filter_list(vx, distances, rgal)
            vy_gal = filter_list(vy, distances, rgal)
            vz_gal = filter_list(vz, distances, rgal)

            print(f"Found {len(x_gal)} star(s) at {percentage}%")

            # Check to see if all the stars have been captured
            if len(x_gal) == num_stars:
                print(f"Found all star(s) at {percentage}%")
                break

            # Increase percentage for next iteration
            percentage += 5

        # All the lists are the same length
        # Loop through and make a list of stars
        stars = []
        for i in range(len(x_gal)):
            star = Star(x_gal[i], y_gal[i], z_gal[i], m_gal[i], a_gal[i], vx_gal[i], vy_gal[i], vz_gal[i])
            stars.append(star)
        # Return the list of stars in the indicated dark matter halo
        return stars


    def get_field(self, field):
        """
        Get the values in the column of the specified field from the .AHF_halos file.

        Parameters:
        ----------
        field : string
            The name of the field.
        
        Returns
        -------
        The list of values in that field.
        """
        # Get the correct name of the field
        field_name = get_field_name(self.ahf_data, field)
        # Store all the field data in a list called column
        column = self.ahf_data.field(field_name) 
        # Return the column
        return column