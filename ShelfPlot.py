# Extracting the provided data from the image into a pandas DataFrame for plotting.
import pandas as pd
import matplotlib.pyplot as plt

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

# Plotting the data
# Since we have two y-axes, we'll create a twin axis
fig, ax1 = plt.subplots()
# Set the x axis to a logarithmic scale
ax1.set_xscale('log')

# Plot Shelf Force
color = 'tab:red'
ax1.set_xlabel('Concentration (mM)')
ax1.set_ylabel('Shelf Force (nN)', color=color)
ax1.scatter(df['Concentration (mM)'], df['Shelf Force (nN)'], color=color)
ax1.tick_params(axis='y', labelcolor=color)

# Create a twin Axes sharing the same x-axis
ax2 = ax1.twinx()

# Plot Shelf Range
color = 'tab:blue'
ax2.set_ylabel('Shelf Range (nm)', color=color)
ax2.scatter(df['Concentration (mM)'], df['Shelf Range (nm)'], color=color)
ax2.tick_params(axis='y', labelcolor=color)

# Show the plot
fig.tight_layout()
plt.show()
