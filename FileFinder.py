import os
import tkinter as tk
from tkinter import filedialog, ttk
from tkinter.font import Font
import re
import pandas as pd
import pyperclip
import subprocess

def browse_folder():
    folder_path = filedialog.askdirectory()
    if folder_path:
        print(f"Selected folder: {folder_path}")  # Debugging
        search_files(folder_path)

def search_files(folder_path):
    file_info_list = []
    seen_files = set()  # Set to track (filename, size) tuples

    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".txt"):
                full_path = os.path.join(root, file)
                file_size = os.path.getsize(full_path)
                file_key = (file, file_size)

                if file_key in seen_files:
                    # Skip this file; it's a duplicate
                    continue
                else:
                    # Add to the set of seen files
                    seen_files.add(file_key)

                if file_size > 15000:  # Check if file size is greater than 15KB
                    is_valid, file_info = is_target_file(file)
                    if is_valid:
                        file_info['has_folder'] = check_folder_existence(root, file)
                        file_info['filepath'] = full_path
                        file_info_list.append(file_info)
                        print(f"Found file: {full_path} of size {file_size}") # Debugging

    display_files(file_info_list)


def check_folder_existence(folder_path, filename):
    # The base name for the corresponding folder seems to be the file's base name plus "_force_curves"
    base_name = os.path.splitext(filename)[0]  # Removes the file extension
    corresponding_folder = f"{base_name}_force_curves"
    return os.path.isdir(os.path.join(folder_path, corresponding_folder))


def is_target_file(filename):
    # Exclude files with undesired patterns
    undesired_patterns = [
        r'pdist-',  # Excludes files starting with 'pdist-'
        r'^LICENSE',  # Excludes files starting with 'LICENSE'
        r'random-double-data',  # Excludes files with 'random-double-data'
        r'iris',  # Excludes files with 'iris'
        r'2018-05-23',  # Excludes files with a date like '2018-05-23'
        # Add any other patterns you want to exclude
    ]
    if any(re.search(pattern, filename) for pattern in undesired_patterns):
        return False, None
    
    # Skip unwanted files (e.g., in "_force_curves" folders)
    if '_force_curves' in filename:
        return False, None
    
    # Exclude files with approach or retract in the filename
    if 'approach' in filename or 'retract' in filename:
        return False, None
    
    # Exclude files with deriv in the filename
    if 'deriv' in filename:
        return False, None    


    # Regex patterns for different parts of the filename
    tip_pattern = r'T(\d+\.\d+)'  # e.g., T6.4
    tip_age_pattern = r'\.(\d+)\.'  # e.g., .2. or .4.
    concentration_pattern = r'(\d+\.?\d*mM)'  # e.g., 10mM or 0.5mM
    site_pattern = r'S\d+'  # e.g., S1
    frequency_pattern = r'(\d+\.?\d+Hz)'  # e.g., 0.5Hz
    dwell_time_pattern = r'(\d+D|nD)'  # e.g., 5D or nD
    ph_pattern = r'pH(\d+)'  # e.g., pH5
    acid_pattern = r'(Hcl)'  # e.g., Hcl
    # Include a pattern to match the force data, like '12nN'
    force_pattern = r'(\d+nN)'  # e.g., 12nN

    # Default values
    defaults = {
        'tip': 'unknown',
        'tip_age': 'unknown',
        'concentration': 'unknown',
        'site': 'unknown',
        'frequency': '0.5Hz',
        'dwell_time': 'nD',
        'ph': 'unknown',
        'acid': 'unknown'
    }

    file_info = defaults.copy()

    if not re.search(site_pattern, filename) and \
       not re.search(frequency_pattern, filename) and \
       not re.search(concentration_pattern, filename):
        return False, None

    # Extracting information using regex
    for pattern, key in [(tip_pattern, 'tip'), (tip_age_pattern, 'tip_age'),
                         (concentration_pattern, 'concentration'), (site_pattern, 'site'),
                         (frequency_pattern, 'frequency'), (dwell_time_pattern, 'dwell_time'),
                         (ph_pattern, 'ph'), (acid_pattern, 'acid')]:
        match = re.search(pattern, filename)
        if match:
            file_info[key] = match.group(1 if key in ['tip_age', 'ph'] else 0)

    force_match = re.search(force_pattern, filename)
    if force_match:
        file_info['force'] = force_match.group(1)
    else:
        file_info['force'] = 'unknown'

    return True, file_info


def on_double_click(event):
    selected_item = treeview.selection()[0]  # Get selected item
    filepath = treeview.item(selected_item)['values'][-1]  # Assuming filepath is the last value
    # Open the folder containing the file in Windows Explorer and select the file
    os.startfile(os.path.normpath(filepath))


def display_files(file_info_list):
    bold_font = Font(family="Helvetica", size=10, weight="bold")
    # Clear existing treeview
    for i in treeview.get_children():
        treeview.delete(i)

    # Define the columns with 'filename' as the first column
    columns = ('filename', 'tip', 'tip_age', 'concentration', 'site', 'frequency', 'dwell_time', 'ph', 'acid', 'processed')
    treeview['columns'] = columns

    # Format the columns
    treeview.column("#0", width=0, stretch=tk.NO)  # Hides the default first column
    treeview.column('filename', anchor=tk.W, width=120)  # Wider column for filename
    for col in columns[1:]:  # Start from second item to avoid 'filename' repetition
        treeview.column(col, anchor=tk.CENTER, width=80)
        # Bind the sort function to the column header
        treeview.heading(col, text=col.capitalize(), anchor=tk.CENTER,
                         command=lambda _col=col: treeview_sort_column(treeview, _col, False))

    # Insert data into the treeview
    for file_info in file_info_list:
        filename_only = os.path.basename(file_info.get('filepath', 'unknown'))  # Extract just the filename
        # Extract the filepath so it's just the containing folder path
        filepath = os.path.dirname(file_info.get('filepath', 'unknown'))
        processed = 'Yes' if file_info['has_folder'] else 'No'
        row = (
            filename_only,
            file_info.get('tip', 'unknown'),
            file_info.get('tip_age', 'unknown'),
            file_info.get('concentration', 'unknown'),
            file_info.get('site', 'unknown'),
            file_info.get('frequency', 'unknown'),
            file_info.get('dwell_time', 'unknown'),
            file_info.get('ph', 'unknown'),
            file_info.get('acid', 'unknown'),
            processed,
            filepath
        )
        # Create a tag for bold font if processed
        tag = 'processed' if file_info['has_folder'] else None
        if tag:
            treeview.tag_configure('processed', font=bold_font)

        treeview.insert('', tk.END, values=row, tags=(tag,))

    # Scrollbar for the Treeview
    scrollbar = ttk.Scrollbar(root, orient=tk.VERTICAL, command=treeview.yview)
    treeview.configure(yscroll=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

def on_file_select(event):
    selected_item = treeview.selection()[0]  # Get selected item
    filepath = treeview.item(selected_item)['values'][-1]  # Assuming filepath is the last column
    pyperclip.copy(filepath)  # Copy filepath to clipboard

def treeview_sort_column(tv, col, reverse):
    # Grab all the entries of the Treeview
    l = [(tv.set(k, col), k) for k in tv.get_children('')]
    # Sort the list depending on if the values are numbers or strings
    try:
        l.sort(key=lambda t: float(t[0]), reverse=reverse)
    except ValueError:  # if values are not numbers, sort as strings
        l.sort(key=lambda t: t[0], reverse=reverse)

    # Rearrange the items in sorted positions
    for index, (val, k) in enumerate(l):
        tv.move(k, '', index)

    # Reverse sort next time
    tv.heading(col, command=lambda: treeview_sort_column(tv, col, not reverse))

def main():
    #Very naughty
    global treeview
    global root
    root = tk.Tk()
    root.title("File Search Tool")
    root.geometry("1024x800")  # Set the window size to 1024x800

    browse_button = tk.Button(root, text="Browse Folder", command=browse_folder)
    browse_button.pack()

    treeview = ttk.Treeview(root)
    treeview.pack(expand=True, fill='both')
    treeview.bind('<Double-1>', on_double_click)  # Bind double-click event

    root.mainloop()

if __name__ == "__main__":
    main()