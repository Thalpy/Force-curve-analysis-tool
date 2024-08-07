import numpy as np

# Provided data
data = {
    "Shelf Force (nN)": [4.040380517, 2.1794461, 3.352544548, 3.55216903, 3.750639964, 
                         1.548608258, 2.764898754, 3.555564534, 1.849905243, 3.204899701, 3.078170327],
    "Shelf Range (nm)": [1.196935142, 1.533807763, 1.371978184, 1.566253863, 1.35238389, 
                         1.935991555, 1.404623789, 1.079712438, 1.407795322, 1.346135175, 1.216184347]
}

# Calculate the average and standard deviation for each category
average_force = np.mean(data["Shelf Force (nN)"])
std_dev_force = np.std(data["Shelf Force (nN)"], ddof=1)

average_range = np.mean(data["Shelf Range (nm)"])
std_dev_range = np.std(data["Shelf Range (nm)"], ddof=1)

print("Average Shelf Force (nN):", average_force)
print("Standard Deviation of Shelf Force (nN):", std_dev_force)
print("Average Shelf Range (nm):", average_range)
print("Standard Deviation of Shelf Range (nm):", std_dev_range)