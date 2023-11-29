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
h = 1e-8  # Separation distance in meters
concentration = 50  # mM
# Assuming a typical Hamaker constant for glass-water interaction
Hamaker_constant = 7.2e-20  # in Joules

# Function to calculate Debye length (1/kappa) based on LiCl concentrationa
def debye_length(concentration, valence=1, temperature=298.15):
    # Convert concentration from mM to mol/m^3
    concentration_mol = concentration * 1e-3 
    # Ionic strength
    I = 0.5 * valence**2 * concentration_mol
    # Debye length (1/kappa)
    #return np.sqrt(temperature * kB * epsilon_r * epsilon_0 / (2 * NA * (e**2) * concentration))
    return np.sqrt((epsilon_0 * epsilon_r * kB * temperature) / (2 * e**2 * I))


# Function to calculate electrostatic repulsion energy with LiCl concentration
def electrostatic_repulsion_energy(h, concentration):
    """Calculate the electrostatic repulsion energy per unit area with LiCl concentration."""
    kappa = 1 / debye_length(concentration)
    #return (epsilon_r * epsilon_0 * Psi_0**2 / kappa) * np.exp(-kappa * h)
    #return ((64 * 3.14 * epsilon_r * epsilon_0 * kB * temperature * Psi_0**2)/(kappa**2)) * np.exp(-kappa * h)
    force_electro = ((2 * 3.14 * radius * epsilon_r * epsilon_0 * Psi_0**2)/ kappa) * np.exp(-kappa * h)
    print("Electrostatic Force: " + str(force_electro))
    return force_electro

# Assuming van der Waals interaction does not change with LiCl concentration for simplicity

# Calculation for a range of LiCl concentrations
#li_cl_concentrations = np.linspace(0.6, 550, 100)  # mM
li_cl_concentrations = [0.6, 1.6, 5, 10, 25, 50, 230, 550]

# Placeholder function for van der Waals energy
def van_der_Waals_energy(h):
    force_vdw = - (Hamaker_constant * radius) / (6 * h ** 2) # IMaSF ch13 eq 13.11b
    print("van der Waals Force: " + str(force_vdw))
    return force_vdw
    #return -Hamaker_constant / (12 * np.pi * h**2)

# Total interaction energy function including LiCl concentration
def total_interaction_energy(h, concentration):
    return van_der_Waals_energy(h) + electrostatic_repulsion_energy(h, concentration)


# Print the results
#print("LiCl concentration:" + str(concentration) + " mM | Debye length: " + str(debye_length(concentration)) + " m | Force at contact: " + str(total_interaction_energy(1e-9, concentration) * 1e9) + " nN")
#1.6 should be 6.088369852 nm

# plot a graph from  30nm away to 1nm away
x = np.linspace(1e-9, 30e-9, 100)
y = total_interaction_energy(x, concentration)
plt.plot(x, y)
plt.xlabel("Distance (m)")
plt.ylabel("Force (N)")
#plt.show()


forces = []
# Calculation of forces and Debye length for each LiCl concentration
for conc in li_cl_concentrations:
    debye_len = debye_length(conc)  # Calculate Debye length
    #force = derjaguin_approximation_force(1e-9, radius, conc)  # Calculate force
    force = total_interaction_energy(1e-9, conc)
    forces.append(force)
    print(f"LiCl concentration: {conc} mM - Debye length: {debye_len} m - Force at contact: {force} N")
    
# plot a graph
plt.figure(figsize=(10, 6))
plt.plot(li_cl_concentrations, forces)
plt.xlabel("LiCl Concentration (mM)")
plt.ylabel("Force at Contact (N)")
plt.title("Force at Contact vs. LiCl Concentration")
plt.grid(True)
plt.show()
