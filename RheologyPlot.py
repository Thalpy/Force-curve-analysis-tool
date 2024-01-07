import matplotlib.pyplot as plt
import numpy as np

# Raw data
concentrations = np.array([0.6, 1.6, 5, 10, 25, 50, 100, 150, 220, 230, 330, 550])
w_ro = np.array([np.nan, 0.39563, 0.32917, 0.36503, 0.26316, np.nan, 0.14561, 0.12239, 0.07074, 0.05637, 0.04542, 0.00964])
w_ro_std = np.array([np.nan, 0.01992, 0.02203, 0.02, 0.01013, np.nan, 0.0086, 0.00876, 0.00803, 0.00371, 0.00518, 0.00752])
w_afm = np.array([0.4358, 0.34416, 0.29362, 0.32169, 0.28249, 0.29531, np.nan, np.nan, np.nan, 0.1415, np.nan, -0.01505])
w_afm_std = np.array([0.11543, 0.0551, 0.05598, 0.07205, 0.13818, 0.06344, np.nan, np.nan, np.nan, 0.03498, np.nan, 0.04546])

# Filter out NaN values for continuous lines
valid_w_ro = ~np.isnan(w_ro)
valid_w_afm = ~np.isnan(w_afm)

# Plotting the data
plt.figure(figsize=(10, 8))

# Plot W* = 2.11 Rσ* (mN/m) with error bars and lines between points
plt.errorbar(concentrations[valid_w_ro], w_ro[valid_w_ro], yerr=w_ro_std[valid_w_ro], fmt='s-', color='black', label="W*=2.11 Rσ* (mN/m)", capsize=5)

# Plot W*(AFM) = Fc/(π R) (mN/m) with error bars and lines between points
plt.errorbar(concentrations[valid_w_afm], w_afm[valid_w_afm], yerr=w_afm_std[valid_w_afm], fmt='o-', color='red', label="W*(AFM)=Fc/(π R) (mN/m)", capsize=5)

# Set the scale of the x-axis to logarithmic
plt.xscale('log')

# Adding labels and title
plt.xlabel('I (mM)')
plt.ylabel('Critical surface energy W* (mN/m)')
plt.title('Critical surface energy vs. LiCl Concentration')

# Adding a black line at W* = 0
plt.axhline(0, color='black', linewidth=1)

# Adding a grid
plt.grid(True, which="both", ls="--", color='gray')

# Set background to white
plt.gca().set_facecolor('white')

# Adding a legend
plt.legend()

# Show the plot
plt.show()
