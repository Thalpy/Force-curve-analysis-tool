import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Load the data
file_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S1\T6.2.550mM.Si.S1deflection_force_curves\Run3'
data = None
try:
    data_path = os.path.join(file_folder, 'approach__force_sep.txt')
    data = pd.read_csv(data_path, delim_whitespace=True, header=None)
except:
    data_path = os.path.join(file_folder, 'approach_force_sep.txt')
    data = pd.read_csv(data_path, delim_whitespace=True, header=None)


x = data.iloc[:, 0]
y = data.iloc[:, 1]

fig, ax = plt.subplots()
line, = ax.plot(x, y, 'o-', picker=5)  # The line through the data points with a picker tolerance
selected_points, = ax.plot([], [], 'ro', markersize=10)  # Points to be highlighted
#limit axis
plt.xlim(-1, 10)
# Set the title and axes labels
ax.set_title('Force vs Z-separation')
ax.set_xlabel('Z-separation (nm)')
ax.set_ylabel('Force (nN)')

# Variables to hold the click count and points
click_count = 0
points = []

# Function to be called when the mouse is clicked
def onclick(event):
    global click_count, points
    ix, iy = event.xdata, event.ydata
    # Find the nearest data point
    distance = np.sqrt((x - ix)**2 + (y - iy)**2)
    nearest_index = distance.argmin()
    nearest_x, nearest_y = x.iloc[nearest_index], y.iloc[nearest_index]

    # Update the points and click count
    click_count += 1
    points.append((nearest_x, nearest_y))
    selected_points.set_data([p[0] for p in points], [p[1] for p in points])
    fig.canvas.draw()

    if click_count == 2:
        # Calculate the average of the y values and the range of x values between the two points
        x_values = x[(x >= points[0][0]) & (x <= points[1][0])]
        y_values = y[x_values.index]
        avg_y = y_values.mean()
        range_x = x_values.max() - x_values.min()
        print(f'Average Force: {avg_y} nN, Drop range: {range_x} nm')
        click_count = 0  # Reset the click count
        points = []  # Reset the saved points
        selected_points.set_data([], [])  # Clear the selected points

# Connect the click event to the onclick function
fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()

