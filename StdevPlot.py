# Importing required libraries for plotting
import matplotlib.pyplot as plt

# Given data
ionic_strengths = [0.6, 1.6, 5, 10, 25, 50, 230, 550]
ionic_strengths_2hz = [0.6, 1.6, 5, 10, 25, 230, 550]
std_dev_0_5Hz = [1.196644336, 0.571262199, 0.580405031, 0.746983266, 1.432576583, 0.657742351, 0.362660309, 0.47126072]
std_dev_2Hz = [1.157, 0.982554324, 0.936915507, 0.986200284, 1.299710865, 0.930326018, 0.038483763]

# Given data

std_dev_0_5Hz_re = [0.325552863, 0.147893543, 0.52553116,0.633612395, 1.388164976, 0.058843011, 0.378613919, 0.24962572]
std_dev_2Hz_re = [0.016, 0.135329228, 0.145614789, 0.273550117, 0.94379341, 0.120702113, 0.179534119]


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
plt.title('Approach Standard Deviation Comparison at Different Tip Speeds Across Changing Ionic Strengths')
plt.grid(True)
plt.show()
plt.close()

#Retrace
# Plotting the standard deviations for 0.5Hz and 2Hz on the same graph
plt.figure(figsize=(10, 5))
# Plotting standard deviations for 0.5Hz with log scale on x-axis
plt.scatter(ionic_strengths, std_dev_0_5Hz_re, color='blue', label='0.5Hz Std Dev')
# Plotting standard deviations for 2Hz with log scale on x-axis
plt.scatter(ionic_strengths_2hz, std_dev_2Hz_re, color='red', label='2Hz Std Dev')

# Setting the x-axis to logarithmic scale
plt.xscale('log')
# Handling missing data (np.nan) by not plotting it
plt.legend()
plt.xlabel('Ionic strength (mM)')
plt.ylabel('Standard deviation')
plt.title('Retrace Standard Deviation Comparison at Different Tip Speeds Across Changing Ionic Strengths')
plt.grid(True)
plt.show()