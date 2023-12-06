import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
epsilon_r = 62  # Relative permittivity for 50:50 water-glycerol mixture
epsilon_0 = 8.854e-12 # Vacuum permittivity in F/m.
radius = 1e-6  # Radius of the sphere in meters (2 um diameter)
Psi_0 = 0.025  # Surface potential in volts (25 mV as an approximation for neutral pH)
temperature = 300  # Room temperature in Kelvin (25°C)
NA = 6.022e23  # Avogadro's number
e = 1.602e-19  # Elementary charge in Coulombs
kB = 1.38064852e-23  # Boltzmann constant in m^2 kg s^-2 K^-1
h = 1e-9  # Separation distance in meters
concentration = 550  # mM
# Assuming a typical Hamaker constant for glass-water interaction
Hamaker_constant = 7.2e-20  # in Joules

circle_geom = 2 * 3.14 * radius

# Function to calculate Debye length (1/kappa) based on LiCl concentrationa
def debye_length(concentration, valence=1, temperature=298.15):
    # Convert concentration from mM to mol/m^3
    concentration_mol = concentration * 1e-3 
    # Ionic strength
    I = 0.5 * valence**2 * concentration_mol
    # Debye length (1/kappa)
    #return np.sqrt(temperature * kB * epsilon_r * epsilon_0 / (2 * NA * (e**2) * concentration))
    return np.sqrt((epsilon_0 * epsilon_r * kB * temperature) / (2000 * NA * e**2 * I))  # IMaSF ch13 eq 13.10


# Function to calculate electrostatic repulsion energy with LiCl concentration
#def electrostatic_repulsion_energ(sep, concentration):
#    """Calculate the electrostatic repulsion energy per unit area with LiCl concentration."""
#    kappa = 1 / debye_length(concentration)
#    force_electro = ((2 * 3.14 * radius * epsilon_r * epsilon_0 * Psi_0**2)/ kappa) * np.exp(-kappa * sep)
#    #force_electro = 64 * (concentration * 1e-3) * kB * T *
#    print("Electrostatic Force: " + str(force_electro))
#    return force_electro

# Function to calculate electrostatic repulsion energy with LiCl concentration
def electrostatic_repulsion_energy_alex(sep, concentration):
    kappa = 1 / debye_length(concentration)
    repulsion = (epsilon_r * epsilon_0 * Psi_0**2) * circle_geom * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(repulsion))
    return repulsion



# Assuming van der Waals interaction does not change with LiCl concentration for simplicity

# Calculation for a range of LiCl concentrations
#li_cl_concentrations = np.linspace(0.6, 550, 100)  # mM
li_cl_concentrations = [0.6, 1.6, 5, 10, 25, 50, 230, 550]

# Placeholder function for van der Waals energy
def van_der_Waals_energy(sep):
    force_vdw = - (Hamaker_constant * radius) / (6 * sep ** 2) # IMaSF ch13 eq 13.11b
    print(f"van der Waals Force for h {sep}: " + str(force_vdw))
    return force_vdw
    #return -Hamaker_constant / (12 * np.pi * h**2)

# Total interaction energy function including LiCl concentration
def total_interaction_energy(sep, concentration):
    return -van_der_Waals_energy(sep) + electrostatic_repulsion_energy_alex(sep, concentration)


# Print the results
#print("LiCl concentration:" + str(concentration) + " mM | Debye length: " + str(debye_length(concentration)) + " m | Force at contact: " + str(total_interaction_energy(1e-9, concentration) * 1e9) + " nN")
#1.6 should be 6.088369852 nm

# plot a graph from  30nm away to 1nm away
#x = np.linspace(1e-9, 10e-9, 100)
#y = total_interaction_energy(x, concentration)
#plt.plot(x, y)
#plt.xlabel("Distance (m)")
#plt.ylabel("Force (N)")
#plt.show()

# Plot a graph of distance vs force with 3 curves; van der Waals, electrostatic, and combined
# generate a range of distances from 1nm to 20 nm NOT LINSPACE, use python range not numpy
#not np.arange
x = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9, 8e-9, 9e-9, 10e-9, 11e-9, 12e-9, 13e-9, 14e-9, 15e-9, 16e-9, 17e-9, 18e-9, 19e-9, 20e-9]
y_vdw = []
y_electro = []
y_total = []
for i in x:
    print("doing " + str(i) + "m now")
    # van der Waals
    y_vdw.append(van_der_Waals_energy(i))
    # electrostatic
    y_electro.append(electrostatic_repulsion_energy_alex(i, concentration))
    # combined
    y_total.append(total_interaction_energy(i, concentration))
#plot curves
plt.figure(figsize=(10, 6))
plt.plot(x, y_vdw, label="van der Waals")
plt.plot(x, y_electro, label="electrostatic")
plt.plot(x, y_total, label="combined")
plt.xlabel("Distance (m)")
plt.ylabel("Force (N)")
# debye as a dot
debye_len = debye_length(concentration)
plt.plot(debye_len, 0, 'o', label="Debye Length")
plt.legend()
plt.show()


forces = []
# Calculation of forces and Debye length for each LiCl concentration
for conc in li_cl_concentrations:
    debye_len = debye_length(conc)  # Calculate Debye length
    #force = derjaguin_approximation_force(1e-9, radius, conc)  # Calculate force
    force = total_interaction_energy(1e-9, conc)
    forces.append(force)
    print(f"LiCl concentration: {conc} mM - Debye length: {debye_len} m - Force at contact: {force} N")

experimental_data = [
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
]
 # LiCl concentration (mM), force at contact (N), standard deviation 

# plot a graph
# make it a log scale on the x axis
plt.figure(figsize=(10, 6))
plt.plot(li_cl_concentrations, forces)
plt.xscale("log")
plt.xlabel("LiCl Concentration (mM)")
plt.ylabel("Force at Contact (N)")
plt.title("Force at Contact vs. LiCl Concentration")
#add the experimental data
for data in experimental_data:
    plt.errorbar(data[0], data[1], yerr=data[2], fmt='o', ecolor='red', capsize=5)
plt.grid(True)
plt.show()
