import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

# Directory containing your data files
data_folder = r'C:\Users\phaso\OneDrive - University of Edinburgh\NewFits\Forcemaps\10mM_force_curves\approach'  # Replace with your folder path
retract_mode = False  # Set to True if processing retract data
output_folder = os.path.join(data_folder, 'output')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Function to process each file
def process_file(file_path, retract=False):
    # Read data
    data = pd.read_csv(file_path, comment='#', delim_whitespace=True, names=['Zpiezo', 'Defl', 'Force'])
    
    # Convert 'Force' to numeric, handling non-numeric entries
    data['Force'] = pd.to_numeric(data['Force'], errors='coerce')
    
    # Drop rows with NaN values (if any non-numeric entries were found)
    data.dropna(subset=['Force'], inplace=True)
    
    # Determine the region to calculate the average force
    if retract:
        # Use the last 50 data points
        avg_force = data['Force'].iloc[-100:-50].mean()
        floor_region_start = data['Zpiezo'].iloc[-100]
        floor_region_end = data['Zpiezo'].iloc[-50]
    else:
        # Use the first 50 data points
        avg_force = data['Force'].iloc[50:100].mean()
        floor_region_start = data['Zpiezo'].iloc[50]
        floor_region_end = data['Zpiezo'].iloc[100]
    
    # Shift the entire force dataset
    data['Adjusted Force'] = data['Force'] - avg_force
    
    # Find the minimum force in the adjusted dataset
    min_force = data['Adjusted Force'].min()
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(data['Zpiezo'], data['Adjusted Force'], label='Adjusted Force')
    plt.axvspan(floor_region_start, floor_region_end, color='yellow', alpha=0.5, label='Flooring Region')
    plt.xlabel('Zpiezo')
    plt.ylabel('Adjusted Force')
    plt.title(f'Adjusted Force Data for {os.path.basename(file_path)}')
    plt.legend()
    plot_file_path = os.path.join(output_folder, os.path.basename(file_path).replace('.txt', '.png'))
    plt.savefig(plot_file_path)
    plt.close()

    return os.path.basename(file_path), min_force


# Turn interactive plotting off
plt.ioff()

# Process each file and collect results
results = []
for file_path in glob.glob(os.path.join(data_folder, '*.txt')):
    filename, min_force = process_file(file_path, retract=retract_mode)
    results.append((filename, min_force))  # Make sure to append a tuple

# Create a DataFrame from results
results_df = pd.DataFrame(results, columns=['Filename', 'Minimum Force'])

# Save results to CSV
results_csv_path = os.path.join(output_folder, 'min_force_results.csv')
results_df.to_csv(results_csv_path, index=False)

# Now extract the 'Minimum Force' column as a NumPy array
min_forces_array = results_df['Minimum Force'].to_numpy()
min_forces_array = min_forces_array[~np.isnan(min_forces_array)]  # Remove NaN values

# Plot histogram
fig, ax = plt.subplots(figsize=(10, 4))
ax.hist(min_forces_array, bins=50)
mean_force = min_forces_array.mean()
std_force = min_forces_array.std()
ax.set_xlabel('Peak Attractive Force (nN)')
ax.set_ylabel('Counts')
ax.set_title(f'Average Peak Attactive Force = {mean_force:.3f} +/- {std_force:.3f} (nN)')

# Save the histogram figure
histogram_path = os.path.join(output_folder, 'contact_force_histogram.png')
plt.tight_layout()
fig.savefig(histogram_path)
plt.close(fig)