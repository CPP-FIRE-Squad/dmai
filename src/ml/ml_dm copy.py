from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from stellarutil.simulation import Simulation
from stellarutil.graph import graph
from stellarutil.calculations import filter_list, dist
import numpy as np

# Step 1: Get the halo at index 0
m10r = Simulation('m10r_res250md', species=['star', 'dark'])
halo = m10r.get_halo()

def get_dm_mass(star):
    # Get the x,y,z positions of each dm particle in the simulation
    # And normalize it with the center of the indicated dark matter halo
    x = m10r.particles['dark']['position'][:,0] - halo.xc
    y = m10r.particles['dark']['position'][:,1] - halo.yc
    z = m10r.particles['dark']['position'][:,2] - halo.zc
    # Get the mass of each star in the simulation
    m = m10r.particles['dark']['mass']
    # Get the distance of each dm from the center of the indicated dark matter halo
    distances =  dist(x,y,z) 
    # Filter out all dm that are farther than the star's r
    dm_masses = filter_list(m, distances, star.get_r())
    # dark matter mass is the mass of each particle whose r < r_star
    return np.sum(dm_masses)

get_dm_mass(halo.stars[0])

# Step 2: Prepare the data for training
X = []
y = []
for star in halo.stars:
    X.append([star.x, star.y, star.vz, ((star.x**2 + star.y**2)**0.5) ])
    y.append(get_dm_mass(star))

# Step 3: Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 4: Initialize and train the KNN regressor
knn_regressor = KNeighborsRegressor(n_neighbors=5)
knn_regressor.fit(X_train, y_train)

# Step 5: Predict on the test data
y_pred = knn_regressor.predict(X_test)

# Step 6 - graph
fig = graph(y_test, y_pred, "Actual vs Predicted DM mass", "Actual", "Predicted", showLine=False, logx=True, logy=True)
fig.plot([0, 1], [0, 1], label='y = x', linestyle='--', color='red')