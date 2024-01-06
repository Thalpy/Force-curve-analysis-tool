# Importing required libraries for plotting
import matplotlib.pyplot as plt

# Given data
ionic_strengths = [0.6, 1.6, 5, 10, 25, 50, 230, 550]
ionic_strengths_2hz = [0.6, 1.6, 5, 10, 25, 230, 550]
std_dev_0_5Hz = [0.325552863, 0.147893543, 0.018316668, 0.003715339, 0.372851406, 0.002227386, 0.078372887, 0.008847494]
std_dev_2Hz = [0.016, 0.135329228, 0.145614789, 0.273550117, 0.94379341, 0.120702113, 0.179534119]


# Plotting the standard deviations for 0.5Hz and 2Hz on the same graph
plt.figure(figsize=(10, 5))
# Plotting standard deviations for 0.5Hz with log scale on x-axis
plt.scatter(ionic_strengths, std_dev_0_5Hz, color='blue', label='0.5Hz Std Dev')
# Plotting standard deviations for 2Hz with log scale on x-axis
plt.scatter(ionic_strengths_2hz, std_dev_2Hz, color='red', label='2Hz Std Dev')

# Setting the x-axis to logarithmic scale
plt.xscale('log')
# Handling missing data (np.nan) by not plotting it
plt.legend()
plt.xlabel('Ionic strength (mM)')
plt.ylabel('Standard deviation')
plt.title('Retract Standard Deviation Comparison at Different Tip Speeds Across Changing Ionic Strengths')
plt.grid(True)
plt.show()
