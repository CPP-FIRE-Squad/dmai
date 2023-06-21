import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pickle

with open("../data/data20_pickle", 'rb') as output:
    data = pickle.load(output)
with open("../data/part20_pickle", 'rb') as output:
    part = pickle.load(output)
with open("../data/header9_pickle", 'rb') as output:
    header = pickle.load(output)


h = header["hubble"]
part = part
header = header

# Initial Conditions
me = np.where(data.nstar > 400)
i = me[0][0]  # Id
R = 0.15 * data[i].rvir / h  # Radius we are looking at

# Centered Position
x = part.x - data[i].xc / h
y = part.y - data[i].yc / h
z = part.z - data[i].zc / h
# Centered Velocity
vz = part.vz - data[i].vcz / h

# Calculate distance and filter based on Restricted Radius
r = np.sqrt(x * x + y * y + z * z)
x = x[r < R]
y = y[r < R]
vz = vz[r < R]

# Perform binned statistic
ret = scipy.stats.binned_statistic_2d(x, y, vz, statistic='std', bins=52)[0]

def graph(inter_technique):
    fig = plt.figure(1, figsize=(10, 10))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.patch.set_facecolor('black')
    vmin = 0
    vmax = 20
    plt.imshow(ret, cmap='hot', vmin=vmin, vmax=vmax, aspect='auto', interpolation=inter_technique)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\rm  \,Velocity Dispersion \,(km/s)$', fontsize=22)
    plt.xlabel(inter_technique, fontsize=22)
    plt.ylabel('y', fontsize=22)
    plt.style.use('default')
    plt.show()


graph(None)
# graph('nearest')
# graph('bilinear')
# graph('spline16')
# graph('spline36')
# graph('hanning')
# graph('hamming')
# graph('hermite')
# graph('quadric')
# graph('kaiser')
# graph('catrom')
graph('gaussian')
# graph('bessel')
# graph('mitchell')
# graph('sinc')
# graph('lanczos')
# graph('blackman')
