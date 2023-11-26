# Importing necessary libraries
import tkinter as tk
from tkinter import filedialog, messagebox, Label, Entry, Button
import os
import re

# Function to extract variables from the file contents
def extract_variables_from_log(text):
    patterns = {
        'name': r"name\s*=\s*'([^']+)'",
        'k_c': r"k_c\s*=\s*([\d.]+)",
        'binsize': r"binsize\s*=\s*([\d.]+)",
        'dfit_win': r"dfit_win\s*=\s*(\d+)",
        'dfit_off': r"dfit_off\s*=\s*(\d+)",
        'cfit_min': r"cfit_min\s*=\s*(\d+)",
        'cfit_max': r"cfit_max\s*=\s*(\d+)",
        'fitbin': r"fitbin\s*=\s*(\d+)",
        'cthresh': r"cthresh\s*=\s*(\d+)",
        'ext': r"ext\s*=\s*'([^']+)'",
        'clear': r"clear\s*=\s*(True|False)",
        'out': r"out\s*=\s*(True|False)",
    }
    variables = {key: re.search(pattern, text).group(1) for key, pattern in patterns.items() if re.search(pattern, text)}
    return variables

# Function to find the appropriate Run folder
def find_run_folder(file_path, mode):
    base_folder = file_path[:-4]
    force_curve_folder = base_folder + '_force_curves'

    if not os.path.exists(force_curve_folder):
        return None

    if mode == "Approach":
        pattern = "Run"
    else:
        pattern = "ret_Run"

    run_folders = [f for f in os.listdir(force_curve_folder) if os.path.isdir(os.path.join(force_curve_folder, f)) and pattern in f]
    
    # Finding the folder with the highest number or "(this one)"
    selected_run_folder = None
    max_run_number = -1
    for folder in run_folders:
        if "(this one)" in folder:
            selected_run_folder = folder
            break
        # Adjusted regex to match the folder naming pattern for both modes
        regex_pattern = r'Run(\d+)' if mode == "Approach" else r'ret_Run(\d+)_Ret'
        run_number_match = re.search(regex_pattern, folder)
        if run_number_match:
            run_number = int(run_number_match.group(1))
            if run_number > max_run_number:
                max_run_number = run_number
                selected_run_folder = folder

    return os.path.join(force_curve_folder, selected_run_folder) if selected_run_folder else None

# GUI Application
class MFPAnalysisApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MFP Analysis Tool")
        self.geometry("400x800")

        self.mode = tk.StringVar(value="Approach")
        self.approach_radio = tk.Radiobutton(self, text="Approach", variable=self.mode, value="Approach")
        self.retract_radio = tk.Radiobutton(self, text="Retract", variable=self.mode, value="Retract")

        self.approach_radio.pack()
        self.retract_radio.pack()

        # UI elements
        self.file_label = Label(self, text="No file selected")
        self.browse_button = Button(self, text="Browse", command=self.browse_file)
        self.run_button = Button(self, text="Run Analysis", command=self.run_analysis, state='disabled')

        # Layout
        self.file_label.pack(pady=10)
        self.browse_button.pack(pady=10)
        self.run_button.pack(pady=10)

        # Variables
        self.selected_file = None
        self.variable_entries = {}  # To store the entry widgets for variables

        self.clear_button = Button(self, text="Clear Fields", command=self.clear_fields)
        self.clear_button.pack(pady=10)
        self.create_variable_fields()

    def create_variable_fields(self):
        self.default_values = {
            'name': '', 'k_c': '0.185', 'binsize': '1.0', 'dfit_win': '45', 'dfit_off': '65',
            'cfit_min': '20', 'cfit_max': '60', 'fitbin': '30', 'cthresh': '50', 'ext': '.jpg',
            'clear': 'False', 'out': 'True'
        }
        variables = ['name', 'k_c', 'binsize', 'dfit_win', 'dfit_off', 'cfit_min', 'cfit_max', 'fitbin', 'cthresh', 'ext', 'clear', 'out']
        for var in variables:
            label = Label(self, text=var)
            entry = Entry(self, name=var)
            label.pack()
            entry.pack()
            self.variable_entries[var] = entry

    def browse_file(self):
        file_path = filedialog.askopenfilename(title="Select a file", filetypes=[("Text Files", "*.txt")])
        if file_path:
            if not file_path.endswith('.txt'):
                messagebox.showerror("Error", "Please select a .txt file")
                return

            self.selected_file = file_path  # Store the selected file path regardless of finding approachparameter_log.txt
            self.file_label.config(text=os.path.basename(file_path))
            self.run_button['state'] = 'normal'

            run_folder = find_run_folder(file_path, self.mode.get())
            if run_folder:
                log_file_path = os.path.join(run_folder, 'approachparameter_log.txt')
                if os.path.exists(log_file_path):
                    try:
                        with open(log_file_path, 'r') as file:
                            variables = extract_variables_from_log(file.read())
                            for var, entry in self.variable_entries.items():
                                if not entry.get():  # Update only if the field is blank
                                    entry.insert(0, variables.get(var, ''))
                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to load file: {e}")
                else:
                    messagebox.showinfo("Info", "approachparameter_log.txt not found. Fields will be left blank.")
            else:
                messagebox.showinfo("Info", "No valid Run folder found. Fields will be left blank.")
        else:
            messagebox.showerror("Error", "No file selected")


    def clear_fields(self):
        for entry in self.variable_entries.values():
            entry.delete(0, tk.END)

    def load_variables(self, file_path):
        # Clear previous entries
        for entry in self.variable_entries.values():
            entry.destroy()
        self.variable_entries.clear()

        # Read file and extract variables
        with open(file_path, 'r') as file:
            variables = extract_variables_from_log(file.read())

        # Create entry fields for each variable
        for var, value in variables.items():
            label = Label(self, text=var)
            entry = Entry(self)
            entry.insert(0, value)
            label.pack()
            entry.pack()
            self.variable_entries[var] = entry

    def run_analysis(self):
        if not self.selected_file:
            messagebox.showerror("Error", "No file selected")
            return
        # Extract values from entries
        variable_values = {var: entry.get() for var, entry in self.variable_entries.items()}

        # Fill in default values if fields are empty
        for var, entry in self.variable_entries.items():
            if not entry.get():
                entry.insert(0, self.default_values.get(var, ''))

        # Find the run folder
        run_folder = find_run_folder(self.selected_file)
        if run_folder:
            if self.mode.get() == "Approach":
                script_name = "run_mfp_analysis.py"
            else:
                script_name = "run_mfp_analysis_ret.py"
            # Here we would call the modified 'run_mfp_analysis.py' script with variable_values and run_folder
            # For now, just showing a message
            messagebox.showinfo("Info", f"Run folder: {run_folder}\nVariables: {variable_values}\nScript: {script_name}")
        else:
            messagebox.showerror("Error", "Run folder not found")

            


# Running the application
if __name__ == "__main__":
    app = MFPAnalysisApp()
    app.mainloop()