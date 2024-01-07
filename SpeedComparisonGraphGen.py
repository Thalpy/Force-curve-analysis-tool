import matplotlib.pyplot as plt
import numpy as np

ionic_strengths_0_5Hz = [0.6, 1.6, 5, 10, 25, 50, 230, 550]
approach_forces_0_5Hz = [4.518, 3.568, 3.044, 3.335, 2.928666667, 3.0615, 1.467, -0.156]
std_dev_0_5Hz = [1.196644336, 0.571262199, 0.580405031, 0.746983266, 1.432576583, 0.657742351, 0.362660309, 0.47126072]

ionic_strengths_0_1Hz = [10]
approach_forces_0_1Hz = [4.064666667]
std_dev_0_1Hz = [0.532789514]

ionic_strengths_2Hz = [0.6, 1.6, 5, 10, 25, 230, 550]
approach_forces_2Hz = [7.557, 4.7655, 3.916, 3.878333333, 3.108333333, 2.4245, -0.0265]
std_dev_2Hz = [1.157, 0.982554324, 0.936915507, 0.986200284, 1.299710865, 0.930326018, 0.038483763]

# Data extracted from the provided text
approach_forces_0_5Hz_re = [-0.152666667, -0.1095, -0.328333333, -0.327666667, -2.377666667, -0.1885, -1.0375, -1.289333333]
std_dev_0_5Hz_re = [0.325552863, 0.147893543, 0.52553116,0.633612395, 1.388164976, 0.058843011, 0.378613919, 0.24962572]

approach_forces_0_1Hz_re = [-0.764]
std_dev_0_1Hz_re = [0.635746543]

approach_forces_2Hz_re = [-0.004, -0.119, -0.109333333, -0.267, -1.29, -0.7295, -1.14]
std_dev_2Hz_re = [0.016, 0.135329228, 0.145614789, 0.273550117, 0.94379341, 0.120702113, 0.179534119]

# Redefining the plotting function to color the error bars the same as the data points
def plot_log_linear(ionic_strengths, approach_forces, std_dev, title, color):
    plt.errorbar(ionic_strengths, approach_forces, yerr=std_dev, fmt='o', ecolor=color, capsize=5, color=color, label=title)
    plt.xscale('log')
    plt.xlabel('LiCl Concentration (mM)')
    plt.ylabel('Force (nN)')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.axhline(0, color='black', linewidth=0.5)

# Plotting the data for all frequencies
plt.figure(figsize=(12, 8))

# Plot for 0.5Hz with purple color as per user's request
plot_log_linear(ionic_strengths_0_5Hz, approach_forces_0_5Hz, std_dev_0_5Hz, '0.5Hz', 'purple')

# Plot for 0.1Hz, using a different color
plot_log_linear(ionic_strengths_0_1Hz, approach_forces_0_1Hz, std_dev_0_1Hz, '0.1Hz', 'green')

# Plot for 2Hz, using a different color and handling NaN values in the data
# Convert NaNs to a number that won't be plotted (out of range)
approach_forces_2Hz = [x if ~np.isnan(x) else -1 for x in approach_forces_2Hz]
plot_log_linear(ionic_strengths_2Hz, approach_forces_2Hz, std_dev_2Hz, '2Hz', 'red')

# Adding a legend to distinguish the different frequencies
plt.legend()

# Title for the overall plot
plt.title('Averaged Approach Force vs. LiCl Concentration at Different Tip Frequencies')

# Show the plot
plt.show()
plt.close()

# Plotting the data for all frequencies
plt.figure(figsize=(12, 8))

# Plot for 0.5Hz with purple color as per user's request
plot_log_linear(ionic_strengths_0_5Hz, approach_forces_0_5Hz_re, std_dev_0_5Hz_re, '0.5Hz', 'purple')

# Plot for 0.1Hz, using a different color
plot_log_linear(ionic_strengths_0_1Hz, approach_forces_0_1Hz_re, std_dev_0_1Hz_re, '0.1Hz', 'green')

# Plot for 2Hz, using a different color and handling NaN values in the data
# Convert NaNs to a number that won't be plotted (out of range)
approach_forces_2Hz = [x if ~np.isnan(x) else -1 for x in approach_forces_2Hz_re]
plot_log_linear(ionic_strengths_2Hz, approach_forces_2Hz_re, std_dev_2Hz_re, '2Hz', 'red')

# Adding a legend to distinguish the different frequencies
plt.legend()
# Title for the overall plot
plt.title('Averaged Retrace Attractive Force vs. LiCl Concentration at Different Tip Frequencies')
plt.show()