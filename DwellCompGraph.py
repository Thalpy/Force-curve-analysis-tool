import matplotlib.pyplot as plt
import numpy as np

approach = False

# standard data
ionic_strengths = np.array([0.6, 1.6, 5, 10, 25, 50, 230, 550])
approach_forces = np.array([4.518, 3.568, 3.044, 3.335, 2.928666667, 3.0615, 1.467, -0.156])
approach_std_dev = np.array([1.196644336, 0.571262199, 0.10209174, 0.12709399, 1.529998003, 0.300820902, 0.08472412, 0.004511646])


retract_forces = [-0.152666667, -0.1095, -0.328333333, -0.327666667, 
                        -2.377666667, -0.1885, -1.0375, -1.289333333]
retract_std_dev = [0.325552863, 0.147893543, 0.018316668, 0.003715339, 
                         0.372851406, 0.002227386, 0.078372887, 0.008847494]

retrace_ionic_strengths_5s = [0.6, 1.6, 5, 10, 25, 230, 550]
retrace_forces_5s = [-1.603, -0.211, -0.614, -0.335, -2.029, -1.892, -2.743]
std_dev_retrace_5s = [0.764, 0.135, 0.654, 0.208, 1.304, 0.336, 0.235]

approach_ionic_strengths_5s = [25, 230]
approach_forces_5s = [1.552, -0.058]
std_dev_approach_5s = [0.742, 0.096]


# Redefining the plotting function to color the error bars the same as the data points
def plot_log_linear(ionic_strengths, approach_forces, std_dev, title, color):
    plt.errorbar(ionic_strengths, approach_forces, yerr=std_dev, fmt='o', ecolor=color, capsize=5, color=color, label=title)
    plt.xscale('log')
    plt.xlabel('LiCl Concentration (mM)')
    plt.ylabel('Force (nN)')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.axhline(0, color='black', linewidth=0.5)

# Plotting the data for all
plt.figure(figsize=(12, 8))

if approach:
    # Plot for standard with purple color for standard
    plot_log_linear(ionic_strengths, approach_forces, approach_std_dev, 'No dwell', 'purple')
    # Plot for 5s, using a different color
    plot_log_linear(approach_ionic_strengths_5s, approach_forces_5s, std_dev_approach_5s, '5s dwell', 'green')
    # Title for the overall plot
    plt.title('Averaged Approach Force vs. LiCl Concentration at Different Dwelling Times')

else:
    # Plot for standard with purple color for standard
    plot_log_linear(ionic_strengths, retract_forces, retract_std_dev, 'No dwell', 'purple')
    plot_log_linear(retrace_ionic_strengths_5s, retrace_forces_5s, std_dev_retrace_5s, '5s dwell', 'green')
    # Title for the overall plot
    plt.title('Averaged Retrace Attractive Force vs. LiCl Concentration at Different Dwelling Times')

# Adding a legend to distinguish the different frequencies
plt.legend()

# Show the plot
plt.show()
