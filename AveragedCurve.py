import matplotlib.pyplot as plt
import numpy as np

#toggle
approach = False

# Data
ionic_strengths = np.array([0.6, 1.6, 5, 10, 25, 50, 230, 550])
approach_force = np.array([4.518, 3.568, 3.044, 3.335, 2.928666667, 3.0615, 1.467, -0.156])
standard_deviations_approach = np.array([1.196644336, 0.571262199, 0.580405031, 0.746983266,
                                         1.432576583, 0.657742351, 0.362660309, 0.47126072])
retrace_force = np.array([-0.152666667, -0.1095, -0.328333333, -0.327666667,
                          -2.377666667, -0.1885, -1.0375, -1.289333333])
standard_deviations_retrace = np.array([0.325552863, 0.147893543, 0.52553116,
                                        0.633612395, 1.388164976, 0.058843011,
                                        0.378613919, 0.24962572])


# Create the plot
plt.figure(figsize=(10, 6))

polynomial = None
if approach:
    # Plot the data with error bars
    plt.errorbar(ionic_strengths, approach_force, yerr=standard_deviations, fmt='o', color='red', ecolor='red', capsize=5)

    # Fit and plot a quadratic curve (second-order polynomial) for a slight curve at the end
    coefficients = np.polyfit(np.log10(ionic_strengths), approach_force, 2)
    polynomial = np.poly1d(coefficients)
    # Creating 100 points from the minimum x to the maximum x to plot the polynomial smoothly
    x_poly = np.linspace(min(ionic_strengths), max(ionic_strengths), 100)
    y_poly = polynomial(np.log10(x_poly))
    plt.plot(x_poly, y_poly, color='red', linestyle='--', linewidth=1)
    plt.title('Averaged Approach Force vs. LiCl Concentration')
else:
    # Plot the data with error bars
    plt.errorbar(ionic_strengths, retrace_force, yerr=standard_deviations_retrace, fmt='o', color='blue', ecolor='blue', capsize=5)

    # Fit and plot a quadratic curve (second-order polynomial) for a slight curve at the end
    coefficients = np.polyfit(np.log10(ionic_strengths), retrace_force, 2)
    polynomial = np.poly1d(coefficients)

    # Creating 100 points from the minimum x to the maximum x to plot the polynomial smoothly
    x_poly = np.linspace(min(ionic_strengths), max(ionic_strengths), 100)
    y_poly = polynomial(np.log10(x_poly))
    plt.plot(x_poly, y_poly, color='blue', linestyle='--', linewidth=1)
    plt.title('Averaged Retract Attractive Force vs. LiCl Concentration')

# Add a horizontal line at y=0 with a thinner line
plt.axhline(0, color='grey', linewidth=1)

# Set the x-axis to a logarithmic scale
plt.xscale('log')

# Add labels and a title to the plot
plt.xlabel('Ionic strength (mM)')
plt.ylabel('Force (nN)')

# Add a grid for easier readability
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Show the plot
plt.show()
