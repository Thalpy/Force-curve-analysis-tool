import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

# Constants
epsilon_r = 62  # Relative permittivity for 50:50 water-glycerol mixture # From that other peper
epsilon_0 = 8.854e-12 # Vacuum permittivity in F/m.
radius = 3.3e-6  # Radius of the sphere in meters (6 um diameter)
Psi_0 = -0.040  # Surface potential in volts (25 mV as an approximation for neutral pH) # -0.04 taken from https://nanocomposix.com/pages/silica-physical-properties
temperature = 300  # Room temperature in Kelvin (25°C)
NA = 6.022e23  # Avogadro's number
e = 1.602e-19  # Elementary charge in Coulombs
kB = 1.38064852e-23  # Boltzmann constant in m^2 kg s^-2 K^-1
h = 2e-9  # Separation distance in meters
concentration = 550  # mM
approach_dict = {
    0.6: r"F:\OneDrive\OneDrive - University of Edinburgh\NewFits\0.6mM\S1\T6.4.0.5mMSi.S1.deflection0000++_force_curves\Run4 (thiS)\approach__force_sep.txt",
    550: r"F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S1\T6.2.550mM.Si.S1deflection_force_curves\Run3\approach__force_sep.txt"
}
retract_dict = {
    0.6: r"F:\OneDrive\OneDrive - University of Edinburgh\NewFits\0.6mM\S1\T6.4.0.5mMSi.S1.deflection0000++_force_curves\ret_Run_Ret1\retract__force_sep.txt",
    550: r"F:\OneDrive\OneDrive - University of Edinburgh\NewFits\550mM\S1\T6.2.550mM.Si.S1deflection_force_curves\ret_Run3_Ret\retract__force_sep.txt"
}
# Assuming a typical Hamaker constant for glass-water interaction
Hamaker_constant = 0.63e-20  # in Joules # From psi_0 values silica.pdf
valency = 1 # aka little z
#Roughness value if needed - leave blank if not needed
roughness_value = 2.6e-10 + 6.45e-10 #From experimental data

#Retrace curve parameters
sep_initial = 1e-10  # 0.1 nm initial separation (close to surface)
adhesion_force = 7.03e-10  # From experimental data
tip_distance = 10e-9  # Retract 10 nm

circle_geom = 2 * 3.14 * radius

# Function to calculate Debye length (1/kappa) based on LiCl concentrationa
def debye_length(concentration, valence=1, temperature=298.15):
    # Convert concentration from mM to mol/m^3
    concentration_mol = concentration # It's already in mol/m^3 
    # Ionic strength
    I = 0.5 * valence**2 * concentration_mol
    # Debye length (1/kappa)
    #return np.sqrt(temperature * kB * epsilon_r * epsilon_0 / (2 * NA * (e**2) * concentration))
    # debye length
    debye_length = np.sqrt(epsilon_r * epsilon_0 * kB * temperature / (2 * NA * e**2 * I))
    print("Debye Length: " + str(debye_length) + " concentration: " + str(concentration))
    return debye_length

# Function to calculate surface charge density at the silica and glass interface
# This is a placeholder function, you would need to replace it with your actual method
# which might include a surface complexation model or experimental data
def calculate_surface_charge_density(pH, pzc):
    # Assuming a simple linear relationship between pH and surface charge density
    # This is a simplification and may not reflect the actual complex behavior
    charge_density = (pH - pzc) * -0.1  # This is an arbitrary scaling factor
    return charge_density

# Function to calculate the Stern layer potential (psi_d)
# For simplicity, let's assume psi_d is a fraction of the surface potential psi_0
# This fraction would typically be determined by fitting to experimental data
def stern_potential(psi_0, fraction=0.1):
    return psi_0 * fraction

# Function to calculate the surface potential psi_0 using the Gouy-Chapman-Stern model
def gouy_chapman_stern(charge_density, kappa, psi_d_guess):
    # Gouy-Chapman equation for psi_0 as a function of surface charge density and Debye length
    # The actual solution of this equation can be complex and often requires numerical methods
    def equation(psi_0):
        return charge_density - (2 * z * e * I) * np.sinh(z * e * psi_0 / (2 * k_B * T)) * np.exp(-kappa * psi_d_guess)
    
    # Solve for psi_0 using a numerical solver
    psi_0_solution = fsolve(equation, psi_d_guess)
    return psi_0_solution[0]

def electrostatic_repulsion_force(sep, concentration, roughness_value = None): #sphere plane fig 14.10 IMaSD
    kappa = 1 / debye_length(concentration)
    Z = 64 * np.pi * epsilon_0 * epsilon_r * (kB * temperature / e)**2 * np.tanh(valency * e * Psi_0 / (4 * kB * temperature))**2
    if roughness_value is not None:
        sep += roughness_value
    force_per_area = kappa * radius * Z * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(force_per_area))
    return force_per_area


# Assuming van der Waals interaction does not change with LiCl concentration for simplicity

# Calculation for a range of LiCl concentrations
#li_cl_concentrations = np.linspace(0.6, 550, 100)  # mM
#li_cl_concentrations = [0.6, 1.6, 5, 10, 25, 50, 230, 550]
li_cl_concentrations = np.arange(6e-1, 5.5e2, 1e-1)  # mM

# Placeholder function for van der Waals energy
def van_der_Waals_force(sep, roughness_value = None): # Sphere - plane fig 13.1 IMaSD
    if roughness_value is not None:
        sep += roughness_value
    force_vdw = - ((Hamaker_constant * radius) / (6 * sep ** 2)) # IMaSF ch13 eq 13.11b
    print(f"van der Waals Force for h {sep}: " + str(force_vdw))
    return force_vdw
    #return -Hamaker_constant / (12 * np.pi * h**2)

def van_der_Waals_force_plane_plane(sep): # plane - plane fig 13.1 IMaSD
    force_vdw = -Hamaker_constant / (6 * np.pi * sep**3)
    print(f"van der Waals Force for h {sep}: " + str(force_vdw))
    return force_vdw

# Total interaction energy function including LiCl concentration
def total_interaction_energy(sep, concentration, roughness_value = None):
    return van_der_Waals_force(sep, roughness_value) + electrostatic_repulsion_force(sep, concentration, roughness_value)

"""
Curve estimation
"""

def calculate_retrace_force(sep_initial, adhesion_force, concentration, retract_distance, roughness_value = None):
    # Initialize arrays to store force and distance values
    forces = []
    distances = np.linspace(sep_initial, sep_initial + retract_distance, num=100)
    
    for sep in distances:
        vdW_force = van_der_Waals_force(sep, roughness_value)
        electro_force = electrostatic_repulsion_force(sep, concentration, roughness_value)
        total_force = vdW_force + electro_force - adhesion_force  # Subtract the adhesion force
        forces.append(-total_force)
    
    return distances, forces

# Function to generate the approach curve
def generate_approach_curve(concentration, start_sep, end_sep, steps, roughness_value = None):
    
    separations = np.linspace(start_sep, end_sep, steps)
    forces = []
    for sep in separations:
        force_vdw = van_der_Waals_force(sep, roughness_value)
        force_el = electrostatic_repulsion_force(sep, concentration, roughness_value)
        total_force = force_vdw + force_el
        forces.append(total_force)
    return separations, forces

distances, forces = calculate_retrace_force(sep_initial, adhesion_force, concentration, tip_distance, roughness_value)

# Plotting the retrace curve
plt.figure()
# plotting the experimental data
if concentration in retract_dict:
    #data = pd.read_csv(retract_dict[concentration], sep='\t', header=None)
    data = pd.read_csv(retract_dict[concentration], delim_whitespace=True, header=None)
    separation_data = data.iloc[:, 0] * 1e-10  # Assuming first column is separation
    force_data = data.iloc[:, 1] * 1e-7      # Assuming second column is force
    plt.scatter(separation_data, force_data, color='red', label='Experimental Data', marker='o')

plt.plot(distances, forces, color='blue', label='Retrace curve')
plt.xlabel('Separation distance (m)')
plt.ylabel('Force (N)')
plt.title('DLVO Retrace Curve')

    #plt.plot(data[0], data[1], label='Retrace curve (experimental)', color='red')
plt.legend()
plt.grid(True)
plt.show()
plt.close()

separations, forces = generate_approach_curve(concentration, tip_distance, sep_initial, 100, roughness_value)
#plt.figure(figsize=(8, 5))
plt.plot(separations, forces, label='Approach curve', color='red')
plt.xlabel('Separation distance (nm)')
plt.ylabel('DLVO Force (N)')
plt.title('DLVO Approach and retrace Curve')
plt.legend()
plt.grid(True)
plt.show()

# Plot a graph of distance vs force with 3 curves; van der Waals, electrostatic, and combined
# generate a range of distances from 1nm to 20 nm NOT LINSPACE, use python range not numpy
#not np.arange
x = np.arange(1e-9, 15e-9, 1e-10)
y_vdw = []
y_electro = []
y_total = []
for i in x:
    print("doing " + str(i) + "m now")
    # If the force decreases past -5e-10, stop plotting
    tot_engy = total_interaction_energy(i, concentration)
    # van der Waals
    y_vdw.append(van_der_Waals_force(i))
    # electrostatic
    y_electro.append(electrostatic_repulsion_force(i, concentration))
    # combined
    
    y_total.append(tot_engy)

#plot curves
plt.figure(figsize=(10, 6))
plt.plot(x, y_vdw, label="van der Waals")
plt.plot(x, y_electro, label="electrostatic")
plt.plot(x, y_total, label="combined")
plt.xlabel("Distance (m)")
plt.ylabel("Force (N)")
#title the graph and add a legend with the concentration
plt.title("Force vs. Distance for LiCl Concentration " + str(concentration) + " mM")
# debye as a straight line downwards

debye_len = debye_length(concentration)
plt.plot([debye_len, debye_len], [-5e-10, 5e-10], label="Debye Length")
#plt.plot(debye_len, 0, 'o', label="Debye Length")
plt.legend()
plt.show()


# Calculate electrostatic repulsive force for a fixed separation distance across different concentrations
fixed_sep = 1e-9  # 1 nm
electrostatic_forces = [electrostatic_repulsion_force(fixed_sep, conc) for conc in li_cl_concentrations]

# Plotting electrostatic force vs concentration
plt.figure(figsize=(10, 6))
plt.plot(li_cl_concentrations, electrostatic_forces, marker='o')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("LiCl Concentration (mM)")
plt.ylabel("Electrostatic Repulsive Force (N)")
plt.title("Electrostatic Repulsive Force vs. LiCl Concentration")
plt.grid(True)
plt.show()

forces = []
# Calculation of forces and Debye length for each LiCl concentration
for conc in li_cl_concentrations:
    debye_len = debye_length(conc)  # Calculate Debye length
    #force = derjaguin_approximation_force(1e-9, radius, conc)  # Calculate force
    # find the
    force = total_interaction_energy(h, conc, roughness_value)
    forces.append(force)
    print(f"LiCl concentration: {conc} mM - Debye length: {debye_len} m - Force at contact: {force} N")

"""experimental_data = [
 [0.6, 5.245e-09, 1.468e-09],
 [1.6, 2.59e-09, 6e-10],
 [5, 2.216e-09, 5.33e-10],
 [10, 1.79e-09, 1.132e-09],
 [25, 2.189e-09, 1.693e-09],
 [50, 3.342e-09, 5.95e-10],
 [230, 2.18e-09, 2.78e-10],
 [550, -2.78e-10, 1.48e-10],
 [0.6, 3.644e-09, 3.32e-10],
 [1.6, 4.546e-09, 5.41e-10],
 [5, 3.133e-09, 7.15e-10],
 [10, 3.832e-09, 4.12e-10],
 [25, 2.746e-09, 1.067e-09],
 [50, 2.781e-09, 7.15e-10],
 [230, 7.54e-10, 4.31e-10],
 [550, -6.6e-11, 6.6e-11],
 [0.6, 4.665e-09, 1.425e-09],
 [5, 3.783e-09, 4.64e-10],
 [10, 4.383e-09, 4.72e-10],
 [25, 3.851e-09, 1.467e-09],
 [550, -1.24e-10, 8e-10]
]"""
 # LiCl concentration (mM), force at contact (N), standard deviation 

# Given data arrays
ionic_strengths = np.array([0.6, 1.6, 5, 10, 25, 50, 230, 550])
approach_force = np.array([4.518, 3.568, 3.044, 3.335, 2.928666667, 3.0615, 1.467, -0.156])
standard_deviations_approach = np.array([1.196644336, 0.571262199, 0.580405031, 0.746983266,
                                         1.432576583, 0.657742351, 0.362660309, 0.47126072])

# Convert the forces and standard deviations to the scientific format and combine them into the desired format
experimental_data = [[i, f, sd] for i, f, sd in zip(ionic_strengths, approach_force * 1e-9, standard_deviations_approach * 1e-9)]

# plot a graph
# make it a log scale on the x axis
plt.figure(figsize=(10, 6))
plt.plot(li_cl_concentrations, forces, label="Calculated Forces")
plt.xscale("log")
plt.xlabel("LiCl Concentration (mM)")
plt.ylabel("Force at Contact (N)")
plt.title("Force at Contact vs. LiCl Concentration")
#add the experimental data
honk = True
for data in experimental_data:
    if honk:
        plt.errorbar(data[0], data[1], yerr=data[2], fmt='o', color="red", ecolor='red', capsize=5, label='Experimental Data')
        honk = False
    plt.errorbar(data[0], data[1], yerr=data[2], fmt='o', color="red", ecolor='red', capsize=5)
plt.grid(True)
# Adding a black line at W* = 0
plt.axhline(0, color='black', linewidth=1)
# Adding a grid
plt.grid(True, which="both", ls="--", color='gray')
plt.legend()
plt.show()
