# Importing necessary libraries
import tkinter as tk
from tkinter import filedialog, messagebox, Label, Entry, Button, Checkbutton
import os
import re
import subprocess
import winsound

"""
    binsize: This parameter defines the size of the bins used for averaging data points. A larger binsize will result in a smoother curve but may obscure finer details.

    dfit_win: This parameter is used in the remove_farfield_drift function. It specifies the window size for fitting a line to a portion of the data far from the contact point to remove any background drift in the deflection signal. It's essentially defining the range of data used for this linear fit.
    The approach phase window.

    dfit_off: Also used in the remove_farfield_drift function, this parameter defines the offset from the guessed contact point to the start of the window for drift fitting. It helps in avoiding the influence of contact-related changes while fitting the drift.

    cfit_min and cfit_max: These parameters define the minimum and maximum deflection values within which the script tries to fit a line to identify the contact point. This fitting is important for accurately determining where the tip makes contact with the surface.

    fitbin: This parameter is used in the binning process (bin_z_df function) and defines the number of data points per bin. It's similar to binsize but specifically for the procedure that aligns the contact region.

    cthresh: This is a threshold for deflection used as a 'first pass' guess to identify the contact point. It helps in determining the initial guess for where the tip makes contact with the surface.

    Moving the Orange Area Left and Right: The orange area likely represents a range of data points used for a specific calculation, possibly related to the background drift removal or finding the contact point. To move this area left or right, you would adjust the parameters that define the window of data points considered for this part of the analysis. This could be the dfit_win and dfit_off parameters in your script, which define the window and offset for fitting the background drift.
    dfit_win = size of window
    dfit_off = offset from guessed contact point (Higher means more to the left)


    Moving the Purple Area Up and Down: The purple area seems to represent data points after the contact point where the tip is interacting with the surface. Moving this area up or down would involve adjusting the baseline or the deflection offset. This could be done by changing the way the background drift is removed or by adjusting the deflection threshold (cthresh) for contact.
    cfit_min = minimum deflection value for fitting a line to identify the contact point
    cfit_max = maximum deflection value for fitting a line to identify the contact point

    Increasing the Dots in the Purple Area: The density of dots in the purple area corresponds to the number of data points in the region of the curve that has been plotted. To increase the number of dots, you would likely need to reduce the fitbin parameter, which controls the number of data points per bin when averaging to smooth out the curve. A smaller fitbin value would result in less averaging and more individual data points being plotted.
"""

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

    # Define the pattern based on the mode
    if mode == "Approach":
        pattern = r'^Run\d+$'  # Pattern for "Approach" mode (e.g., "Run14")
    else:
       pattern = r'^ret_Run_Ret\d+.*$'  # Pattern for "Retract" mode (e.g., "ret_Run_Ret5")

    # Filter run folders based on the defined pattern
    run_folders = [f for f in os.listdir(force_curve_folder) if os.path.isdir(os.path.join(force_curve_folder, f)) and re.match(pattern, f)]

    # Finding the folder with the highest number or "(this one)"
    selected_run_folder = None
    max_run_number = -1
    for folder in run_folders:
        if "(this one)" in folder:
            selected_run_folder = folder
            break
        run_number_match = re.search(r'\d+', folder)
        if run_number_match:
            run_number = int(run_number_match.group(0))
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

        self.extra_var = tk.BooleanVar()
        self.extra_checkbox = Checkbutton(self, text="Extra", variable=self.extra_var)
        self.extra_checkbox.pack()

        self.dwell_var = tk.BooleanVar()  # Variable to track the state of the dwell checkbox
        self.dwell_checkbox = Checkbutton(self, text="Dwell", variable=self.dwell_var)
        self.dwell_checkbox.pack()  # Add the dwell checkbox to the GUI



    def create_variable_fields(self):
        # Define default values for 'Approach' mode
        approach_defaults = {
            'name': '', 'k_c': '0.185', 'dfit_win': '45', 'dfit_off': '65',
            'cfit_min': '20', 'cfit_max': '60', 'fitbin': '30', 'cthresh': '50',
            'ext': '.jpg', 'clear': 'False', 'out': 'True', 'binsize': '1.0'
        }
        
        # Define default values for 'Retract' mode
        retract_defaults = {
            'name': '', 'k_c': '0.184', 'dfit_win': '40', 'dfit_off': '85',
            'cfit_min': '20', 'cfit_max': '50', 'fitbin': '25', 'cthresh': '50',
            'ext': '.jpg', 'clear': 'False', 'out': 'True', 'binsize': '1.0'
        }
        
        # Use the correct set of default values based on the current mode
        if self.mode.get() == "Approach":
            self.default_values = approach_defaults
        else:  # Retract mode
            self.default_values = retract_defaults
        
        # Variables to create entry fields for
        variables = ['name', 'k_c', 'dfit_win', 'dfit_off', 'cfit_min', 'cfit_max',
                    'fitbin', 'cthresh', 'ext', 'clear', 'out', 'binsize']
        
        for var in variables:
            label = Label(self, text=var)
            entry = Entry(self, name=var)
            #entry.insert(0, self.default_values.get(var, ''))
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
                log_file_name = 'retract_parameter_log.txt' if self.mode.get() == 'Retract' else 'approachparameter_log.txt'
                log_file_path = os.path.join(run_folder, log_file_name)
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
        extra = self.extra_var.get()
        dwell = self.dwell_var.get()  # Get the state of the dwell checkbox

        # Fill in default values if fields are empty
        for var, entry in self.variable_entries.items():
            if not entry.get():
                entry.insert(0, self.default_values.get(var, ''))

        # Update variable_values with the new (or default) values
        variable_values = {var: entry.get() for var, entry in self.variable_entries.items()}

        # Set 'approach' based on the mode
        variable_values['approach'] = self.mode.get() == "Approach"

        # Prepare arguments for subprocess
        script_name = "run_mfp_analysis.py" # if self.mode.get() == "Approach" else "run_mfp_analysis_ret.py"
        args = [
            self.selected_file,  # sys.argv[1]
            str(variable_values.get('k_c', '')),  # sys.argv[2]
            str(variable_values.get('fitbin', '')),  # sys.argv[3]
            str(variable_values.get('cfit_min', '')),  # sys.argv[4]
            str(variable_values.get('cfit_max', '')),  # sys.argv[5]
            str(variable_values.get('cthresh', '')),  # sys.argv[6]
            str(variable_values.get('dfit_win', '')),  # sys.argv[7]
            str(variable_values.get('dfit_off', '')),  # sys.argv[8]
            str(variable_values.get('binsize', '')),  # sys.argv[9]
            str(variable_values.get('ext', '')),  # sys.argv[10]
            str(variable_values.get('out', '')).lower(),  # sys.argv[11]
            str(variable_values.get('clear', '')).lower(),  # sys.argv[12]
            str(variable_values.get('approach', '')).lower(),  # sys.argv[13]
            str(extra).lower(),  # sys.argv[14]
            str(dwell).lower()  # sys.argv[15] - New argument for dwell
        ]


        # Prepare confirmation message
        confirmation_message = f"Do you want to run {script_name} with the following parameters?\n\n"
        confirmation_message += "\n".join([f"{var}: {variable_values[var]}" for var in ['name', 'k_c', 'fitbin', 'cfit_min', 'cfit_max', 'cthresh', 'dfit_win', 'dfit_off', 'binsize', 'ext', 'approach', 'clear', 'out']])

        # Ask user to confirm running the script
        if messagebox.askyesno("Confirm", confirmation_message):
            # Run the script if user confirms
            subprocess.run(["python", script_name] + args)
            winsound.Beep(100, 100)
            winsound.Beep(200, 100)
            winsound.Beep(300, 100)
            winsound.Beep(400, 100)
            winsound.Beep(500, 100)
            winsound.Beep(700, 500)
            messagebox.showinfo("Info", f"Script {script_name} has been run with the selected parameters.")

        else:
            winsound.Beep(500, 400)
            winsound.Beep(200, 800)
            messagebox.showinfo("Cancelled", "Script execution cancelled.")





            


# Running the application
if __name__ == "__main__":
    app = MFPAnalysisApp()
    app.mainloop()