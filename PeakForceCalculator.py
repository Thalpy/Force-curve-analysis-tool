import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

# Directory containing your data files
data_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S3\550mMS1.0.5Hz.nD.12nN_force_curves\approach'  # Replace with your folder path
retract_mode = False  # Set to True if processing retract data
Zoom_window = True
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

    # Zoomed in plot
    # Use the window of data from 10nN on the y axis to a range of 1um on the x axis
    if Zoom_window:
        plt.figure(figsize=(10, 6))
        plt.plot(data['Zpiezo'], data['Adjusted Force'], label='Adjusted Force')
        #plt.axvspan(floor_region_start, floor_region_end, color='yellow', alpha=0.5, label='Flooring Region')
        plt.xlabel('Zpiezo')
        plt.ylabel('Adjusted Force')
        plt.title(f'Adjusted Force Data for {os.path.basename(file_path)}')
        plt.legend()
        plt.ylim(0, 10)
        #Find the x value where y passes 10nN
        x = data['Zpiezo'].loc[data['Adjusted Force'] > 10].iloc[0]
        #Set the x axis to a range of 0.1um from that value backwards
        plt.xlim(x-0.1, x)
        plot_file_path = os.path.join(output_folder, os.path.basename(file_path).replace('.txt', '_zoom.png'))
        plt.savefig(plot_file_path)
        plt.close()

    return os.path.basename(file_path), min_force


# Turn interactive plotting off
plt.ioff()

# Process each file and collect results
results = []
for file_path in glob.glob(os.path.join(data_folder, '*.txt')):
    result = process_file(file_path, retract=retract_mode)
    results.append(result)

# Save results to CSV
results_df = pd.DataFrame(results, columns=['Filename', 'Minimum Force'])
results_csv_path = os.path.join(output_folder, 'min_force_results.csv')
results_df.to_csv(results_csv_path, index=False)

