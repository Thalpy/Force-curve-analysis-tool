# Let's load the data from the text file and plot the first two columns
import pandas as pd
import matplotlib.pyplot as plt
import os
# If the data contains noise, finding the second dip directly through derivative changes might not be effective.
# Instead, we can smooth the data first to reduce the impact of noise.
# Then, we'll look for local minima in the smoothed data.
from scipy.signal import savgol_filter, find_peaks

# Load the data from the text file
file_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S3\550mMS1.0.5Hz.nD.12nN_force_curves\Run1'
# Let's load both data files and find the points closest to 0 on the x-axis for each, and print those indexes.
# Load the data from the text files
deriv_path = os.path.join(file_folder, r'approach_force_curves\deriv.txt')
approach_path = os.path.join(file_folder, 'approach__force_sep.txt')
spring_constant = 0.156  # nN/nm

# Load the data
deriv_data = pd.read_csv(deriv_path, delim_whitespace=True, header=None)
approach_data = pd.read_csv(approach_path, delim_whitespace=True, header=None)

# Function to read the data, smooth it, and find the index of the second dip
def find_second_dip_smoothed(file_path):
    data = pd.read_csv(deriv_path, delimiter='\s+', header=None)
    # Apply Savitzky-Golay filter to smooth the data
    smoothed = savgol_filter(data.iloc[:, 1], window_length=51, polyorder=2)  # Window length and polyorder might need adjustment

    # After smoothing, we find the local minima
    # We're looking for minima, so we invert the smoothed data
    inverted_smoothed = -smoothed
    peaks, _ = find_peaks(inverted_smoothed, prominence=1)  # Prominence might need adjustment
    # The second dip will be the second peak in the inverted data
    if len(peaks) > 1:
        second_dip_index = peaks[1]  # Second peak index
        second_dip_value = data.iloc[second_dip_index, :]
    else:
        second_dip_index = None
        second_dip_value = None
    
    return second_dip_index, second_dip_value

# Find the second dips for the deriv data file
deriv_second_dip_idx, deriv_second_dip_val = find_second_dip_smoothed(deriv_path)

# Output the results
print(f"Deriv file second dip index: {deriv_second_dip_idx}, value: {deriv_second_dip_val.values if deriv_second_dip_idx else 'Not found'}")

# Repeat the same process for the approach__force_sep data file
approach_second_dip_idx, approach_second_dip_val = find_second_dip_smoothed('/mnt/data/approach__force_sep.txt')
print(f"Approach file second dip index: {approach_second_dip_idx}, value: {approach_second_dip_val.values if approach_second_dip_idx else 'Not found'}")


# Function to read the data and find the index of the point closest to zero on the x-axis
def find_closest_to_zero(file_path):
    data = pd.read_csv(file_path, delimiter='\s+', header=None, usecols=[0, 1])
    # Find the index of the point where the x value (column 0) is closest to zero
    index_closest_to_zero = (data.iloc[:, 0].abs()).idxmin()
    return index_closest_to_zero, data.iloc[index_closest_to_zero, :]

# Find the indexes and values for both files
deriv_idx, deriv_value = find_closest_to_zero(deriv_path)
approach_idx, approach_value = find_closest_to_zero(approach_path)

# Print the results
print(f"Deriv file index closest to zero: {deriv_idx}, value: {deriv_value.values}")
print(f"Approach file index closest to zero: {approach_idx}, value: {approach_value.values}")


# Read the file content assuming it's a three-column format with the first two columns as the data of interest
try:
    deriv_data = pd.read_csv(deriv_path, delimiter='\s+', header=None, usecols=[0, 1])
    curve_data = pd.read_csv(approach_path, delimiter='\s+', header=None, usecols=[0, 1])
    deriv_data = deriv_data * spring_constant # Convert to force
    # Plot the first two columns
    plt.figure(figsize=(10, 5))
    plt.plot(deriv_data[0], deriv_data[1], marker='o', linestyle='-')
    #plt.plot(curve_data[0], curve_data[1], marker='o', linestyle='-')
    plt.title('Force vs deriv def')
    plt.xlabel('Force (nN)')
    plt.ylabel('Deriv def (nm)')
    plt.grid(True)
    plt.show()
except Exception as e:
    print(f"An error occurred: {e}")
