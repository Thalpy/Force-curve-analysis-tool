import os
import numpy as np

def process_folder(folder_path):
    min_forces = []
    
    # List all files in the directory
    file_list = [file for file in os.listdir(folder_path) if file.endswith('_retract.txt')]
    
    # Process each file
    for file_name in file_list:
        file_path = os.path.join(folder_path, file_name)
        try:
            # Load the data from the file
            data = np.loadtxt(file_path)

            # Normalise the left side of the data (assuming the force is in the thrid column)
            data[:, 2] -= np.mean(data[:50, 1])  # Adjust the 50 here if more points are needed

            # Find the minimum force in the dataset
            min_force = np.min(data[:, 2])  # Assuming force is in the third column

            # Append the minimum force to the list
            min_forces.append(min_force)

        except Exception as e:
            print(f"Error processing file {file_name}: {e}")

    # Print out the list of minimum forces
    print(min_forces)
    return min_forces

# Replace 'your_folder_path' with the path to the directory containing your files
your_folder_path = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\Forcemaps\Output_force_curves\retract'
min_forces_list = process_folder(your_folder_path)
