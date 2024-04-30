import matplotlib.pyplot as plt, numpy as np, joblib, os
from sklearn.model_selection import train_test_split 
from stellarutil import Simulation, Star 
from concurrent.futures import ThreadPoolExecutor
from IPython.display import clear_output

def get_sim_data(SIMULATION, HALO_START=0, HALO_END=100):
    X = []
    Y = []
    # Get the halo at index 0 restricted at 100%
    sim = Simulation(simulation_name=SIMULATION, species=['star', 'dark'])
    clear_output()
    # Get the x,y,z positions of each dm particle in the simulation, normalize it with halo center
    dark_x = sim.particles['dark']['position'][:,0]
    dark_y = sim.particles['dark']['position'][:,1]
    dark_z = sim.particles['dark']['position'][:,2]
    # Get the mass of each dm particle in the simulation
    dark_m = sim.particles['dark']['mass']
    # For each halo, get the stars and dark matter particles (X and Y)
    for i in range(HALO_START, HALO_END):
        halo = sim.get_halo(i)
        # Get the x,y,z positions of each dm particle in the simulation, normalize it with halo center
        halo_dark_x = dark_x - halo.xc
        halo_dark_y = dark_y - halo.yc
        halo_dark_z = dark_z - halo.zc
        # Get the distance of each dm particle from the center of the indicated dark matter halo
        halo_dark_distances = np.sqrt(np.square(halo_dark_x) + np.square(halo_dark_y) + np.square(halo_dark_z))
        # Get X - the features used by the ML to predict Y
        for star in halo.stars:
            X.append([star.x, star.y, star.vz, star.a, star.get_3DR(), star.get_3DR()])
        # Find the dark matter masses of each star in the halo, using multiple threads
        with ThreadPoolExecutor(max_workers=12) as executor:
            def get_dm_mass(star: Star):
                dm_masses = dark_m[halo_dark_distances < star.get_3DR()] # Filter out all dm that are farther than the star's r
                return np.sum(dm_masses) # dark matter mass is the mass of each particle whose r < r_star

            def process(star: Star):
                dm_mass = get_dm_mass(star)
                return dm_mass
            
            results = executor.map(process, halo.stars)
        # Get Y - the result of the multithreaded process
        for result in enumerate(results, start=1):
            index, data = result
            Y.append(data)
        
        print(f'Finished processing halo {i+1}/{HALO_END}')
    
    return np.array(X), np.array(Y)

def get_saved_data():
    X_BV = joblib.load('../data/pickle/big_victor/X_BV_0_to_100.pkl')
    Y_BV = joblib.load('../data/pickle/big_victor/Y_BV_0_to_100.pkl')
    X_LV = joblib.load('../data/pickle/little_victor/X_LV_0_to_100.pkl')
    Y_LV = joblib.load('../data/pickle/little_victor/Y_LV_0_to_100.pkl')
    X_LR = joblib.load('../data/pickle/little_romeo/X_LR_0_to_100.pkl')
    Y_LR = joblib.load('../data/pickle/little_romeo/Y_LR_0_to_100.pkl')
    X = X_BV + X_LV + X_LR
    Y = Y_BV + Y_LV + Y_LR
    return np.array(X), np.array(Y)

def split_data(X, Y, test_size=0.2, random_state=42):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=random_state)
    X_TRAIN = np.array(X_train) 
    Y_TRAIN = np.array(Y_train)
    X_TEST = np.array(X_test)
    Y_TEST = np.array(Y_test)
    return X_TRAIN, X_TEST, Y_TRAIN, Y_TEST

def graph(x,y, title, r):
    x = np.array(x)
    y = np.array(y)
    # Create the scatter plot
    plt.scatter(x, y, label='Data Points', c=r, vmin=0.011, vmax=1.5)
    plt.colorbar()
    # Get the max and min value
    minVal = min(min(x), min(y))
    maxVal = max(max(x), max(y))
    # Plot y=x line
    plt.plot([minVal, maxVal], [minVal, maxVal], color='green', label='y = x')
    # Add labels and legend
    plt.xlabel('Actual Mass [M☉]')
    plt.ylabel('Predicted Mass [M☉]')
    plt.title(title)
    plt.legend()
    plt.loglog()
    # Show the plot
    plt.show()

def accuracy(Y_TEST, Y_PRED):
    # Find the percent difference
    sum = 0
    for i in range(len(Y_PRED)):
        a = Y_PRED[i]
        b = Y_TEST[i]
        p = round(100 * ((2 * abs(a-b)) / (a+b)), 3)
        sum = sum + p
    # Accuracy is 100 - average percent difference
    return round(100 - (sum / len(Y_PRED)), 3)

