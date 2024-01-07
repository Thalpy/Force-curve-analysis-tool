import matplotlib.pyplot as plt
import numpy as np

# Toggle for showing approach or retrace graph
show_approach = True  # Set to False to show retrace graph

# Updated data
ionic_strengths = np.array([0.6, 1.6, 5, 10, 25, 50, 230, 550])
approach_force = np.array([4.518, 3.568, 3.044, 3.335, 2.928666667, 3.0615, 1.467, -0.156])
standard_deviations_approach = np.array([1.196644336, 0.571262199, 0.580405031, 0.746983266,
                                         1.432576583, 0.657742351, 0.362660309, 0.47126072])
retrace_force = np.array([-0.152666667, -0.1095, -0.328333333, -0.327666667,
                          -2.377666667, -0.1885, -1.0375, -1.289333333])
standard_deviations_retrace = np.array([0.325552863, 0.147893543, 0.52553116,
                                        0.633612395, 1.388164976, 0.058843011,
                                        0.378613919, 0.24962572])

# Plotting
plt.figure(figsize=(10, 6))
if show_approach:
    # Calculate weights as the inverse of the standard deviations for approach data
    weights_approach = 1 / (standard_deviations_approach + 1e-10)
    # Fit a polynomial to the approach force data using the weights
    coefficients_approach = np.polyfit(np.log(ionic_strengths), approach_force, deg=3, w=weights_approach)
    polynomial_approach = np.poly1d(coefficients_approach)
    # Creating points for the polynomial
    y_poly_approach = polynomial_approach(np.log(ionic_strengths))
    # Plot the approach data with error bars
    plt.errorbar(ionic_strengths, approach_force, yerr=standard_deviations_approach, fmt='o', color='red', ecolor='red', capsize=5, label='Approach')
    # Plot the polynomial trendline
    plt.plot(ionic_strengths, y_poly_approach, 'r--', linewidth=1, label='Approach Fit')
    plt.title('Averaged Approach Force vs. LiCl Concentration')
else:
    # Calculate weights as the inverse of the standard deviations for retrace data
    weights_retrace = 1 / (standard_deviations_retrace + 1e-10)
    # Fit a polynomial to the retrace force data using the weights
    coefficients_retrace = np.polyfit(np.log(ionic_strengths), retrace_force, deg=3, w=weights_retrace)
    polynomial_retrace = np.poly1d(coefficients_retrace)
    # Creating points for the polynomial
    y_poly_retrace = polynomial_retrace(np.log(ionic_strengths))
    # Plot the retrace data with error bars
    plt.errorbar(ionic_strengths, retrace_force, yerr=standard_deviations_retrace, fmt='o', color='blue', ecolor='blue', capsize=5, label='Retrace')
    # Plot the polynomial trendline
    plt.plot(ionic_strengths, y_poly_retrace, 'b--', linewidth=1, label='Retrace Fit')
    plt.title('Averaged Retrace Retension Force vs. LiCl Concentration')

# Add a horizontal line at y=0 with a thinner line
plt.axhline(0, color='black', linewidth=1)

# Set the x-axis to a logarithmic scale
plt.xscale('log')

# Add labels and a title to the plot
plt.xlabel('Ionic strength (mM)')
plt.ylabel('Force (nN)')


# Add a grid for easier readability
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Add a legend to the plot
plt.legend()

# Show the plot
plt.show()
