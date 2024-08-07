import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Data extracted from the image
data = {
    "Concentration (mM)": [1.6, 5, 5, 10, 10, 25, 50, 50, 230, 550, 550],
    "Shelf Force (nN)": [4.040380517, 2.1794461, 3.352544548, 3.55216903, 3.750639964, 
                         1.548608258, 2.764898754, 3.555564534, 1.849905243, 3.204899701, 3.078170327],
    "Shelf Range (nm)": [1.196935142, 1.533807763, 1.371978184, 1.566253863, 1.35238389, 
                         1.935991555, 1.404623789, 1.079712438, 1.407795322, 1.346135175, 1.216184347]
}

# Create DataFrame
df = pd.DataFrame(data)

# Group by concentration and calculate mean and standard deviation
grouped = df.groupby('Concentration (mM)').agg(['mean', 'std'])

# Plotting the data with standard deviation error bars
fig, ax1 = plt.subplots()
ax1.set_xscale('log')  # Set the x-axis to a logarithmic scale

# Plot Shelf Force with error bars
color_force = 'tab:red'
ax1.set_xlabel('Concentration (mM)')
ax1.set_ylabel('Shelf Force (nN)', color=color_force)
ax1.errorbar(grouped.index, grouped['Shelf Force (nN)']['mean'], yerr=grouped['Shelf Force (nN)']['std'], 
             fmt='o', color=color_force, ecolor=color_force, elinewidth=2, capsize=5)
ax1.tick_params(axis='y', labelcolor=color_force)

# Create a twin Axes sharing the same x-axis for Shelf Range
ax2 = ax1.twinx()

# Plot Shelf Range with error bars
color_range = 'tab:blue'
ax2.set_ylabel('Shelf Range (nm)', color=color_range)
ax2.errorbar(grouped.index, grouped['Shelf Range (nm)']['mean'], yerr=grouped['Shelf Range (nm)']['std'], 
             fmt='o', color=color_range, ecolor=color_range, elinewidth=2, capsize=5)
ax2.tick_params(axis='y', labelcolor=color_range)

fig.tight_layout()
plt.show()
