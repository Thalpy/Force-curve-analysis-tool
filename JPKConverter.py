import os
import matplotlib.pyplot as plt

def convert_data(data_lines):
    converted_data = []
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 2:
            zpiezo, defl = parts[0], parts[1]
            zpiezo_um = float(zpiezo) * 1E6  # Convert Zpiezo to micrometers
            force = 0.155 * float(defl) * 1E9  # Calculate the force
            converted_data.append([zpiezo_um, force])
    return converted_data

def save_plot(data, filename):
    zpiezo, force = zip(*data)
    plt.figure()
    plt.plot(zpiezo, force)
    plt.xlabel('Zpiezo (um)')
    plt.ylabel('Force (N)')
    plt.title('Zpiezo vs Force')
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

        # Extracting the last three numbers from the filename
        file_number = file_path.split('_')[-1].split('.')[0]

        # Saving the converted data and plots
        approach_file = os.path.join(export_folder, 'approach', f"{file_number}_approach.txt")
        retract_file = os.path.join(export_folder, 'retract', f"{file_number}_retract.txt")
        plot_filename = os.path.join(export_folder, 'Raw_force_curves', f"{file_number}_force_curve.png")

        header = "# Zpiezo\tForce\n"
        with open(approach_file, 'w') as file:
            file.write(header)
            for zpiezo, force in approach_data:
                file.write(f"{zpiezo:.6e}\t{force:.12e}\n")
        with open(retract_file, 'w') as file:
            file.write(header)
            for zpiezo, force in retract_data:
                file.write(f"{zpiezo:.6e}\t{force:.12e}\n")

        save_plot(approach_data, plot_filename)
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
