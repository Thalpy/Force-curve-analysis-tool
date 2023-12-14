import pandas as pd
import tkinter as tk
from tkinter import filedialog

def bin_data(file_path, binned_file_path, bin_size=10):
    # Load the data
    data = pd.read_csv(file_path, sep="\t", header=0)  # Adjust separator if needed

    # Binning the data
    binned_data = data.groupby(data.index // bin_size).mean()

    # Saving the binned data
    binned_data.to_csv(binned_file_path, sep="\t", index=False)

def browse_file():
    # Set up the root Tkinter window
    root = tk.Tk()
    root.withdraw()

    # Open the file dialog
    file_path = filedialog.askopenfilename()
    
    # Check if a file was selected
    if file_path:
        # File selected, proceed with binning
        binned_file_path = file_path.rsplit('.', 1)[0] + '_binned.txt'
        bin_data(file_path, binned_file_path)
        print(f"Binned file saved as: {binned_file_path}")
    else:
        print("No file selected. Exiting.")

# Run the file browser
browse_file()
