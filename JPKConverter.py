import os
import matplotlib.pyplot as plt

def convert_data(data_lines):
    converted_data = []
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 2:
            zpiezo, defl = parts[0], parts[1]
            zpiezo_um = -float(zpiezo) * 1E6  # Convert Zpiezo to micrometers
            force = 0.155 * float(defl) * 1E9  # Calculate the force
            converted_data.append([zpiezo_um, defl, force])  # Keep original zpiezo and defl
    return converted_data

def save_plot(approach_data, retract_data, filename):
    # Unpack the zpiezo, defl, and force data for approach and retract
    approach_zpiezo, approach_defl, approach_force = zip(*approach_data)
    retract_zpiezo, retract_defl, retract_force = zip(*retract_data)

    plt.figure()
    plt.plot(approach_zpiezo, approach_force, label='Approach')
    plt.plot(retract_zpiezo, retract_force, label='Retract')
    plt.xlabel('Zpiezo (um)')
    plt.ylabel('Force (N)')
    plt.title('Zpiezo vs Force')
    plt.legend()
    # Set the number of x-axis ticks
    plt.locator_params(axis='x', nbins=10)  # Adjust the number of bins as needed

    # Automatically use tight layout to prevent label overlap
    plt.tight_layout()

    # Use scientific notation for the x-axis
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Save the figure and close the plot to release memory
    plt.savefig(filename)
    plt.close()

def process_file(file_path, export_folder):
    try:
        with open(file_path, 'r') as file:
            contents = file.readlines()

        # Splitting data into approach and retract segments
        approach_lines, retract_lines = [], []
        current_segment = None
        for line in contents:
            if line.startswith("# segment:"):
                current_segment = line.strip().split(": ")[1]
                continue
            if line.startswith('#') or not line.strip():
                continue
            if current_segment == 'extend':
                approach_lines.append(line.strip())
            elif current_segment == 'retract':
                retract_lines.append(line.strip())

        approach_data = convert_data(approach_lines)
        retract_data = convert_data(retract_lines)
        # Since the data is backwards, we reverse the lists
        #approach_data.reverse()
        #retract_data.reverse()

        # Extracting the last three numbers from the filename
        file_number = file_path.split('_')[-1].split('.')[0]

        # Saving the converted data and plots
        approach_file = os.path.join(export_folder, 'approach', f"{file_number}_approach.txt")
        retract_file = os.path.join(export_folder, 'retract', f"{file_number}_retract.txt")
        plot_filename = os.path.join(export_folder, 'Raw_force_curves', f"{file_number}_force_curve.png")

        header = "# Zpiezo\tDefl\tForce\n"
        with open(approach_file, 'w') as file:
            file.write(header)
            for zpiezo, defl, force in approach_data:
                file.write(f"{zpiezo}\t{defl}\t{force}\n")  # Zpiezo and Defl are written as-is
                
        with open(retract_file, 'w') as file:
            file.write(header)
            for zpiezo, defl, force in retract_data:
                file.write(f"{zpiezo}\t{defl}\t{force}\n")

        save_plot(approach_data, retract_data, plot_filename)
        print(f"Processed and saved: {file_number}")

    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

def process_folder(input_folder, export_folder):
    # Creating export subfolders if they don't exist
    os.makedirs(os.path.join(export_folder, 'approach'), exist_ok=True)
    os.makedirs(os.path.join(export_folder, 'retract'), exist_ok=True)
    os.makedirs(os.path.join(export_folder, 'Raw_force_curves'), exist_ok=True)

    for filename in os.listdir(input_folder):
        if filename.endswith('.txt'):
            file_path = os.path.join(input_folder, filename)
            process_file(file_path, export_folder)


# Example usage
input_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\Silica Files for John\NewAFM (JPK)\0.5mM\petri\Maps\map-data-2019.06.16-17.35.14.780_processed-2023.12.20-07.25.51\processed_curves-2023.12.20-07.25.51'  # Replace with your input folder path
export_folder = r'F:\OneDrive\OneDrive - University of Edinburgh\NewFits\Forcemaps\Output_force_curves'  # Replace with your export folder path
process_folder(input_folder, export_folder)
