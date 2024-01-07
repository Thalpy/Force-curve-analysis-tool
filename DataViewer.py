import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Load the data
file_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S3\550mMS1.0.5Hz.nD.12nN_force_curves\Run1'
data_path = os.path.join(file_folder, 'approach__force_sep.txt')
data = pd.read_csv(data_path, delim_whitespace=True, header=None)

x = data.iloc[:, 0]
y = data.iloc[:, 1]


fig, ax = plt.subplots()
ax.plot(x, y, 'o-')
# limit the x and y axis
ax.set_xlim([-1, 10])


# Function to be called when the mouse is clicked
def onclick(event):
    ix, iy = event.xdata, event.ydata
    # Find the nearest data point
    distance = np.sqrt((x - ix)**2 + (y - iy)**2)
    nearest_index = distance.argmin()
    nearest_x, nearest_y = x.iloc[nearest_index], y.iloc[nearest_index]
    print(f'Nearest x = {nearest_x}, y = {nearest_y}')


# Connect the click event to the onclick function
fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
