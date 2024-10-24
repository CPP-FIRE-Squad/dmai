import matplotlib.pyplot as plt, numpy as np

# These two graph functions will be used to produce the plots for the 2D Star Positions and Predicted vs Actual Dark Matter Masses
def graph(x,y,c, title, Y, Y_pred):
    # Create a figure
    fig = plt.figure(figsize=(15, 5))
    plt.suptitle(title, fontsize=16)
    # 2D Scatter Plot
    plt.subplot(1, 3, 1) # rows, cols, index
    plt.title("2D Star Positions")
    plt.scatter(x, y, label='Data Points', c=c, vmin=0, vmax=100)
    cb = plt.colorbar()
    cb.set_label('Percent Difference [%]')
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')
    # Predicted vs Actual
    plt.subplot(1, 3, 2)
    plt.title("Predicted vs Actual Dark Matter Masses")
    plt.scatter(Y, Y_pred, label='Data Points', c=c, vmin=0, vmax=100)
    cb = plt.colorbar()
    cb.set_label('Percent Difference [%]')
    # Get the max and min value and plot y=x line
    minVal = min(min(Y), min(Y_pred))
    maxVal = max(max(Y), max(Y_pred))
    plt.plot([minVal, maxVal], [minVal, maxVal], color='green', label='y = x')
    # Add labels and legend
    plt.xlabel('Actual Mass [M☉]')
    plt.ylabel('Predicted Mass [M☉]')
    plt.legend()
    plt.loglog()
    # Bin
    plt.subplot(1, 3, 3).hist(c, bins=20)
    plt.title("Histogram of Percent Differences")
    plt.xlabel('Percent Difference [%]')
    plt.ylabel('Frequency')
    # Show the plots
    plt.tight_layout()
    txt=f"Tested with {len(Y)} stars. The average percent difference is {np.mean(c):.2f}%."
    plt.figtext(0.5, -0.02, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.show()
    fig.savefig(f"../images/{title}_plot.png", bbox_inches="tight")

def graph2(x1, y1, x2, y2, title):
    # Create a figure
    fig = plt.figure(figsize=(15, 5))
    plt.suptitle(title, fontsize=16)
    plt.subplot(1, 2, 1) # rows, cols, index
    plt.scatter(x1, y1, color='red', label='Predicted')
    plt.scatter(x2, y2, color='blue', label='Actual')
    plt.xlabel('2D Radius [kpc]')
    plt.ylabel('DM Mass [M☉]')
    plt.yscale('log')
    plt.legend()
    # Calculate the percentage differences
    percent_diff = np.zeros(len(y2))
    non_zero_indices = np.logical_and(y2 != 0, y1 != 0)
    percent_diff[non_zero_indices] = (np.abs(y2[non_zero_indices] - y1[non_zero_indices]) / ((y2[non_zero_indices] + y1[non_zero_indices]) / 2)) * 100
    # Bin the x1 values and calculate the average percentage difference for each bin
    bins = np.linspace(np.min(x1), np.max(x1), 20)
    bin_indices = np.digitize(x1, bins)
    bin_averages = [percent_diff[bin_indices == i].mean() for i in range(1, len(bins))]
   # Plot the histogram of the average percentage differences
    plt.subplot(1, 2, 2)
    plt.hist(bins[:-1], bins, weights=bin_averages, edgecolor='black')
    plt.title("Histogram of Average Percent Differences")
    plt.xlabel('2D Radius [kpc]')
    plt.ylabel('Average % Difference')
    plt.show()
    fig.savefig(f"../images/{title}_curve.png", bbox_inches="tight")


# This graph3 function is used to generate graphs without a title or caption for the paper
def graph3(x1, y1, x2, y2, title, Y, Y_pred, c):
    # Create a figure
    fig = plt.figure(figsize=(15, 5))
    # Graph 1 - Predicted vs Actual
    plt.subplot(1, 2, 1)
    plt.scatter(Y, Y_pred, label='Data Points', c=c, vmin=0, vmax=100)
    cb = plt.colorbar()
    cb.set_label('Percent Difference [%]')
    # Get the max and min value and plot y=x line
    minVal = min(min(Y), min(Y_pred))
    maxVal = max(max(Y), max(Y_pred))
    plt.plot([minVal, maxVal], [minVal, maxVal], color='green', label='y = x')
    # Add labels and legend
    plt.xlabel('Actual Mass [M☉]')
    plt.ylabel('Predicted Mass [M☉]')
    plt.legend()
    plt.loglog()
    # Save graph 1
    # fig.savefig(f"../images/{title}_p_vs_a.png", bbox_inches="tight")
    # fig = plt.figure(figsize=(15, 5))
    
    # Graph 2 - Mass vs Radius
    plt.subplot(1, 2, 2) # rows, cols, index
    plt.scatter(x1, y1, color='red', label='Predicted')
    plt.scatter(x2, y2, color='blue', label='Actual')
    plt.xlabel('2D Radius [kpc]')
    plt.ylabel('DM Mass [M☉]')
    plt.yscale('log')
    plt.legend()
    
    # Save graph 2
    plt.show()
    fig.savefig(f"../images/{title}_2dr.png", bbox_inches="tight")
    fig = plt.figure(figsize=(15, 5))

    # Graph 3 - Plot the histogram of the average percentage differences
    # Calculate the percentage differences
    percent_diff = np.zeros(len(y2))
    non_zero_indices = np.logical_and(y2 != 0, y1 != 0)
    percent_diff[non_zero_indices] = (np.abs(y2[non_zero_indices] - y1[non_zero_indices]) / ((y2[non_zero_indices] + y1[non_zero_indices]) / 2)) * 100
    # Bin the x1 values and calculate the average percentage difference for each bin
    bins = np.linspace(np.min(x1), np.max(x1), 20)
    bin_indices = np.digitize(x1, bins)
    bin_averages = [percent_diff[bin_indices == i].mean() for i in range(1, len(bins))]
    plt.subplot(1, 2, 2)
    plt.hist(bins[:-1], bins, weights=bin_averages, edgecolor='black')
    # plt.title("Histogram of Average Percent Differences")
    plt.xlabel('2D Radius [kpc]')
    plt.ylabel('Average % Difference')
    # Save graph 3
    plt.show()
    fig.savefig(f"../images/{title}_2dr_histogram.png", bbox_inches="tight")



def graph4(title, Y, Y_pred, c, Y2, Y_pred2, c2):
    # Create a figure
    fig = plt.figure(figsize=(15, 5))

    # Graph 1 - Predicted vs Actual 2d
    plt.subplot(1, 2, 1)
    plt.scatter(Y, Y_pred, label='Data Points', c=c, vmin=0, vmax=100)
    cb = plt.colorbar()
    cb.set_label('Percent Difference [%]')
    # Get the max and min value and plot y=x line
    minVal = min(min(Y), min(Y_pred))
    maxVal = max(max(Y), max(Y_pred))
    plt.plot([minVal, maxVal], [minVal, maxVal], color='green', label='y = x')
    # Add labels and legend
    plt.xlabel('Actual Mass [M☉]')
    plt.ylabel('Predicted Mass [M☉]')
    plt.legend()
    plt.loglog()

    # Graph 2 - Predicted vs Actual 3d
    plt.subplot(1, 2, 2)
    plt.scatter(Y2, Y_pred2, label='Data Points', c=c2, vmin=0, vmax=100)
    cb = plt.colorbar()
    cb.set_label('3D Percent Difference [%]')
    # Get the max and min value and plot y=x line
    minVal = min(min(Y2), min(Y_pred2))
    maxVal = max(max(Y2), max(Y_pred2))
    plt.plot([minVal, maxVal], [minVal, maxVal], color='green', label='y = x')
    # Add labels and legend
    plt.xlabel('Actual Mass [M☉]')
    plt.ylabel('Predicted Mass [M☉]')
    plt.legend()
    plt.loglog()
    
    # Save graph
    plt.show()
    fig.savefig(f"../images/{title}_2dr.png", bbox_inches="tight")
    fig = plt.figure(figsize=(15, 5))


