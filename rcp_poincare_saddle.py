
"""
Created on Sat Sep 21 14:39:18 2024

@author: ar.ray
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp


time_series_data = np.loadtxt('tseries_pre.txt')                                             

colors = plt.cm.viridis(np.linspace(0, 1, 100))  # Different colors for each time series
# Initialize the plot
plt.figure(figsize=(10, 8))
all_crossings = []

for j in range(100):
    print(j)
    # Extract the ith time series columns
    time_series_0 = time_series_data[j*1000000:(j+1)*1000000, :]
    
    time_series_1 = time_series_0[:, [0, 1, 2]]
    
    # Filter time series based on condition z > 0.5
    condition = time_series_1[:, 2] > 0.5
    filtered_time_series = time_series_1[condition]

    # Exclude the first and last 10%
    n_rows = filtered_time_series.shape[0]
    start_index = int(0.2 * n_rows)
    end_index = int(0.8 * n_rows)
    time_series = filtered_time_series[start_index:end_index]

    # Define Poincaré section parameters
    poincare_y = 0.3
    crossings_x = []
    crossings_z = []

    # Detect crossings
    for i in range(1, len(time_series)):
        y_prev, y_next = time_series.T[1, i-1], time_series.T[1, i]
        
        # Check for crossing y = 0.3
        if (y_prev - poincare_y) * (y_next - poincare_y) < 0:
            # Interpolate crossing points
            x_crossing = np.interp(poincare_y, [y_prev, y_next], [time_series.T[0, i-1], time_series.T[0, i]])
            z_crossing = np.interp(poincare_y, [y_prev, y_next], [time_series.T[2, i-1], time_series.T[2, i]])
            crossings_x.append(x_crossing)
            crossings_z.append(z_crossing)

    all_crossings.append([crossings_x, crossings_z])
    # Plotting the Poincaré section for the current time series
    plt.plot(crossings_x, crossings_z, 'o', markersize=2, alpha=0.7, color=colors[j])

# Save all crossings data (x and z) to a single text file
all_crossings = np.array(all_crossings, dtype=object)  # Convert to a numpy array of objects
np.savetxt('all_crossings_1.txt', all_crossings, fmt='%s')


# Finalize the plot
plt.title(f'Poincaré Sections of RCP Ecological System at y = {poincare_y} for 100 Time Series')
plt.xlabel('x')
plt.ylabel('z')
plt.legend()
plt.grid(True)
plt.show()