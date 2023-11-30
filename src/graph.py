import re, numpy as np,  matplotlib.colors as mcolors, matplotlib.pyplot as plt


def histogram(data, bins, title  = 'Historgram', x_label = 'Value', y_label = 'Frequency'):
    # Create a histogram
    plt.hist(data, bins=bins)
    # Add labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    # Show the histogram
    plt.show()

def stars_scatter_plot(stars1, stars2):
    '''

    '''
    # Create figure
    ax = plt.axes(projection='3d')
    colors1 = ['blue' for star in stars1]
    colors2 = ['red' for star in stars1]
    ax.scatter([star.x for star in stars1], [star.y for star in stars1], [star.z for star in stars1], colors1)
    ax.scatter([star.x for star in stars2], [star.y for star in stars2], [star.z for star in stars2], colors2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Scatter Plot')
    ax.grid(False)
    # Show figure
    plt.show()

def star_scatter_plot(stars, parts = None):
    '''
    Generate a scatter plot by supplying a list of x, y, and z values.

    Parameters:
    ----------
        x : list
            The list of x values.
        y : list
            The list of y values.
        z : list
            The list of z values.
        parts : list
            list of values (like form.scalar factor) associated with each star
        whole : number
            the whole number to compare to each part

    '''

    if parts is not None:
        # Normalize the parts values for colormap
        norm = mcolors.Normalize(vmin=0, vmax=0.5)
        # Create a colormap for the gradient
        cmap = plt.cm.get_cmap('viridis')
        # Get colors from the colormap based on normalized parts values
        colors = [cmap(norm(part)) for part in parts]

    # Create figure
    ax = plt.axes(projection='3d')
    ax.scatter([star.x for star in stars], [star.y for star in stars], [star.z for star in stars], c=colors, marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Scatter Plot')
    ax.grid(False)
    # Create color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label('Parts Value')
    # Show figure
    plt.show()

def scatter_plot(x, y, z, parts = None):
    '''
    Generate a scatter plot by supplying a list of x, y, and z values.

    Parameters:
    ----------
        x : list
            The list of x values.
        y : list
            The list of y values.
        z : list
            The list of z values.
        parts : list
            list of values (like form.scalar factor) associated with each star
        whole : number
            the whole number to compare to each part

    '''

    if parts is not None:
        # Normalize the parts values for colormap
        norm = mcolors.Normalize(vmin=0, vmax=0.5)
        # Create a colormap for the gradient
        cmap = plt.cm.get_cmap('viridis')
        # Get colors from the colormap based on normalized parts values
        colors = [cmap(norm(part)) for part in parts]

    # Create figure
    ax = plt.axes(projection='3d')
    ax.scatter(x, y, z, c=colors, marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Scatter Plot')
    ax.grid(False)
    # Create color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label('Parts Value')
    # Show figure
    plt.show()

def pie_chart(values, labels):
    '''
    Generate a pie chart by supplying a list of labels and values.

    Parameters:
    ----------
        labels : list
            The list of labels.
        values : list
            The list of values.
    '''
    sizes = [count for count in values if count > 3]
    explode = [0.1 if size == max(sizes) else 0 for size in sizes]
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', explode=explode)
    ax.legend(loc="best")
    plt.show()

def graph(x, y, title = None, x_label = None, y_label = None, windowTitle = "Figure", showLine = True, logx = False, logy = False):
    '''
    Generate a 2d graph with given x and y values.

    Parameters:
    ----------
        labels : list
            The list of labels.
        values : list
            The list of values.
    '''
    # Auto generate x_label and y_label with the given title
    if title is not None:
        pattern = r'\s*(vs\.|versus|vs)\s*'
        match = re.search(pattern, title, re.IGNORECASE)
        if match:
            before = title[:match.start()].strip()
            after = title[match.end():].strip()
            x_label = x_label if x_label is not None else after
            y_label = y_label if y_label is not None else before

    # Log the lists
    if logx:
        x = np.log10(x)
    if logy:
        y = np.log10(y)

    # Create a line plot
    if(showLine):
        #find line of best fit
        a, b = np.polyfit(x, y, 1)    
        # 'b.-' - This parameter is a shorthand notation to specify the line color, marker style, and line style of the plot.
        plt.plot(x, y, 'b.-')
    else:
        plt.scatter(x,y)

    # Set the x and y axis labels
    plt.xlabel(x_label if x_label is not None else 'X')
    plt.ylabel(y_label if y_label is not None else 'Y')
    # Set the title of the graph
    plt.title(title if title is not None else 'Y vs X')
    fig = plt.gcf()
    fig.canvas.manager.set_window_title(windowTitle)
    # Display the graph
    plt.show()
    # Return the plot
    return plt



def plot_gas_dens(part,data,header,runname,bins=70,ind=0,thick=None,vmin=None,vmax = None):
    """

    """
    import scipy

    ################### LOAD PARTICLE INFO #########################
    h = header['hubble']

    ###########################################

    maincen = np.array([data.field('Xc(6)')[ind]/h,data.field('Yc(7)')[ind]/h,data.field('Zc(8)')[ind]/h])
    mainrvir = data.field('Rvir(12)')[ind]/h

    if thick is None:
        thick = mainrvir

    pos = part['gas']['position'][:]
    pos = pos - maincen
    mass = part['gas']['mass'][:]


    x = pos[:,0] 
    y = pos[:,1] 
    z = pos[:,2] 

    x = x[np.abs(z) < thick]
    y = y[np.abs(z) < thick]
    mass = mass[np.abs(z) < thick]

    mass = mass[((np.abs(x)<=thick)&(np.abs(y)<=thick))]
    xx = x[((np.abs(x)<=thick)&(np.abs(y)<=thick))]
    yy = y[((np.abs(x)<=thick)&(np.abs(y)<=thick))]

    ret = scipy.stats.binned_statistic_2d(xx, yy, mass, statistic='sum', bins=bins)

    # Calc area per bin
    dx = np.diff(ret.x_edge)
    dy = np.diff(ret.y_edge)
    area = dx[:,  None] * dy


    plt.clf()
    fig = plt.figure(1, figsize = (10,10))
    ax1 = fig.add_subplot(1,1,1)
    ax1.patch.set_facecolor('black')

    if vmin is None:
        plt.imshow(np.log10(ret.statistic/area), cmap='seismic', interpolation='gaussian')
    else:
        plt.imshow(np.log10(ret.statistic/area), cmap='seismic', interpolation='gaussian',vmin = vmin,vmax = vmax)

    cbar = plt.colorbar()
    cbar.ax.setylabel(r'$\rm log ,\Sigma ,(M{\odot}/pc^2)$', fontsize = 22)
    plt.savefig("FDensHist%s%i%i.png"%(runname,thick,bins))
 