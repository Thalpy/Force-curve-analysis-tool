import numpy as np
import matplotlib.pyplot as plt

# Data for pH 5
sites_pH5 = [1, 2, 3, 4]
concentration_pH5 = ['pH 5', 'pH 5', 'pH 5', 'pH 5']
Fc_pH5 = [3.017, 4.444, 4.857, 5.92]
stdev_Fc_pH5 = [0.528, 0.799, 0.461, 0.696]
Fa_pH5 = [-0.593, -0.237, -2.698, -0.229]
stdev_Fa_pH5 = [0.47, 0.157, 0.392, 0.202]

# Data for 0.6 mM LiCl concentration
sites_LiCl = [1, 2, 3]
concentration_LiCl = [0.6, 0.6, 0.6]
approach_force_LiCl = [5.245, 3.644, 4.665]
stdev_approach_LiCl = [1.468, 0.332, 1.425]
retrace_force_LiCl = [-0.298, -0.037, -0.123]
stdev_retrace_LiCl = [0.509, 0.083, 0.228]

# Plotting the data for pH 5 and 0.6 mM LiCl concentration with error bars for approach and retrace forces
# We will place all the data points for each pH value on top of each other in the plot

# Preparing the plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 7))

# Approach forces
axes[0].errorbar(np.full_like(Fc_pH5, 5), Fc_pH5, yerr=stdev_Fc_pH5, fmt='o', label='pH 5', color='blue', ecolor='lightblue', capsize=5, capthick=2)
axes[0].errorbar(np.full_like(approach_force_LiCl, 7), approach_force_LiCl, yerr=stdev_approach_LiCl, fmt='o', label='0.6 mM LiCl', color='green', ecolor='lightgreen', capsize=5, capthick=2)

# Retrace forces
axes[1].errorbar(np.full_like(Fa_pH5, 5), Fa_pH5, yerr=stdev_Fa_pH5, fmt='o', label='pH 5', color='red', ecolor='pink', capsize=5, capthick=2)
axes[1].errorbar(np.full_like(retrace_force_LiCl, 7), retrace_force_LiCl, yerr=stdev_retrace_LiCl, fmt='o', label='0.6 mM LiCl', color='orange', ecolor='lightcoral', capsize=5, capthick=2)

# Setting titles and labels for each subplot
axes[0].set_title('Approach Repulsive Force between two pH concentrations')
axes[0].set_xlabel('pH')
axes[0].set_ylabel('Force (nN)')

axes[1].set_title('Retrace Attractive Force between two pH concentrations')
axes[1].set_xlabel('pH')
axes[1].set_ylabel('Force (nN)')

# Adding horizontal line at Force = 0 for reference and legends
for ax in axes:
    ax.axhline(0, color='grey', linewidth=0.5)
    ax.set_xticks([5, 7])
    ax.set_xticklabels(['pH 5', 'pH 7'])
    ax.legend()

plt.tight_layout()
plt.show()
