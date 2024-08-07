# Now that we have the correct data format and the graph as a reference, we can proceed with the analysis.
# We will load the data, then find the yield part of the graph and the yield force.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, argrelextrema, savgol_filter

# Load the data from the text file
data_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S3\550mMS1.0.5Hz.nD.12nN_force_curves\Run1'
data_path = os.path.join(data_folder, 'approach__force_sep.txt')
data = pd.read_csv(data_path, delimiter='\s+', header=None, usecols=[0, 1])

# Since the data goes from right to left, we need to invert the zsep values
zsep = data.iloc[:, 0].values[::-1]
force = data.iloc[:, 1].values[::-1]

# Now we will find the index of the first significant force drop which is the yield point
# We will identify it by finding where the force drops to a small value after the peak
# For this, we can look for the first index where force falls below a certain threshold after the maximum force

# Find the maximum force which is the yield force and its index
yield_force_idx = np.argmax(force)
yield_force = force[yield_force_idx]
yield_zsep_start = zsep[yield_force_idx]

# Define a small threshold to determine when the force drops significantly
# Here, the threshold is chosen to be a small fraction of the yield force
# Adjust this as necessary to suit the specific characteristics of your data
threshold = yield_force * 0.05

# Find the index where the force first drops below this threshold after the yield point
drop_indices = np.where(force[yield_force_idx:] < threshold)[0]

if drop_indices.size > 0:
    # Correct the index to account for the offset after yield_force_idx
    first_drop_idx = drop_indices[0] + yield_force_idx
    yield_zsep_end = zsep[first_drop_idx]
    yield_width_nm = yield_zsep_end - yield_zsep_start
else:
    # If no drop below the threshold is found, we can't determine the width
    first_drop_idx = None
    yield_zsep_end = None
    yield_width_nm = None

# Plot the graph with the identified yield force and the yield part
plt.figure(figsize=(14, 7))

i_start = list(i for i,x in enumerate(zsep) if x>yield_zsep_end)[0]-1
i_end = 0

#smoothed_force = savgol_filter(force, window_length=11, polyorder=3)
#plt.plot(zsep[i_end:i_start], smoothed_force[i_end:i_start], marker='', linestyle='-', color='pink')
#plt.plot(zsep[i_end:i_start], force[i_end:i_start], marker='o', linestyle='-', color='blue')

#grad1 = np.gradient(smoothed_force, zsep)
#grad2 = np.gradient(grad1, zsep)
#plt.plot(zsep[i_end:i_start], grad2[i_end:i_start], marker='o', linestyle='-', color='blue')

n_snap_start = 31
n_snap_end   = 28

plt.plot(zsep[i_end:i_start], force[i_end:i_start], marker='o', linestyle='-', color='blue')
plt.title('Force vs Z-separation')
plt.xlabel('Z-separation (nm)')
plt.ylabel('Force (nN)')
plt.axvline(x=yield_zsep_start, color='r', linestyle='--', label='Start of Yield (Max Force)')
if first_drop_idx is not None:
    plt.axvline(x=yield_zsep_end, color='g', linestyle='--', label='End of Yield (Drop below threshold)')

plt.axvline(x=zsep[n_snap_start], color='blue', linestyle=':', label=f'zsep[{n_snap_start}] - {zsep[n_snap_start]} {force[n_snap_start]}')
plt.axvline(x=zsep[n_snap_end], color='black', linestyle=':', label=f'zsep[{n_snap_end}] - {zsep[n_snap_end]} {force[n_snap_end]}')

plt.legend()
plt.gca().invert_xaxis()  # Invert x-axis to match the original graph's direction
plt.show()

# Output the results
print(f"Yield Force: {yield_force} nN")
print(f"Z-separation at Yield Force: {yield_zsep_start} nm")
print(f"Width of the Yield Part: {yield_width_nm} nm")