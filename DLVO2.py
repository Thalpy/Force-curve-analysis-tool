import numpy as np
import matplotlib.pyplot as plt

# Constants
epsilon_0 = 8.854187817e-12  # Permittivity of free space (F/m)
k_B = 1.380649e-23  # Boltzmann constant (J/K)

# Parameters
radius_sphere = 1e-6  # Sphere radius (m)
distances = np.linspace(30e-9, 1e-9, 1000)  # Separation distances from 30 nm to 1 nm (m)
Hamaker_constant = 1e-19  # Hamaker constant (J)

# LiCl concentration (mol/m^3)
LiCl_concentration = 550  # Converted from 0.6 mM to mol/m^3

# Temperature (K)
temperature = 298

# Calculate Debye length (m)
Debye_length = np.sqrt(epsilon_0 * k_B * temperature / (2 * LiCl_concentration * 1.60219e-19))

# Initialize arrays to store forces
repulsion = np.zeros_like(distances)
attraction = np.zeros_like(distances)
DLVO_force = np.zeros_like(distances)

# Calculate electrostatic repulsion and van der Waals attraction for each separation distance
for i, distance in enumerate(distances):
    repulsion[i] = (1 / (4 * np.pi * epsilon_0)) * (1 / Debye_length) * (np.exp(-distance / Debye_length) / distance)
    attraction[i] = -Hamaker_constant / (6 * distance)
    DLVO_force[i] = attraction[i] + repulsion[i]

# Plot the DLVO force as a function of separation distance
plt.figure(figsize=(10, 6))
plt.plot(distances * 1e9, DLVO_force, label='DLVO Force')
plt.xlabel('Separation Distance (nm)')
plt.ylabel('Force (N)')
plt.title('DLVO Force vs. Separation Distance')
plt.grid(True)
plt.legend()
plt.show()

def electrostatic_repulsion_force_plane_plane(sep, concentration): #plane plane fig 14.10 IMaSD
    kappa = 1 / debye_length(concentration)
    Z = 64 * np.pi * epsilon_0 * epsilon_r * (kB * temperature / e)**2 * np.tanh(valency * e * Psi_0 / (4 * kB * temperature))**2
    force_per_area = (kappa / (2 * np.pi)) * Z * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(force_per_area))
    return force_per_area

def electrostatic_repulsion_force_wowow(sep, concentration): #e-11 # Alt fig 14.10 IMaSD sphere - plane
    kappa = 1 / debye_length(concentration)
    rho_inf = NA * concentration
    #Calculate the electrostatic force per unit area using the Derjaguin approximation for a sphere-plane interaction
    W_flat = ((64 * kB * temperature * rho_inf * np.tanh(1 * e * Psi_0 / (4 * kB * temperature))**2) / kappa ) * np.exp(-kappa * sep)
    total_force = 2 * np.pi * radius * W_flat
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

def electrostatic_repulsion_force_nope(sep, concentration): #e-1
    kappa = 1 / debye_length(concentration)
    rho_inf = NA * concentration
    #Calculate the electrostatic force per unit area using the Derjaguin approximation for a sphere-plane interaction
    total_force = np.pi * epsilon_r * radius  * Psi_0**2 * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

def electrostatic_repulsion_force_oh(sep, concentration): #e-6 ???
    kappa = 1 / debye_length(concentration)
    total_force = 2 * epsilon_0 * sep * kappa ** 2 * Psi_0 ** 2 * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

# Updated function to calculate electrostatic repulsion force between a sphere and a plane
def electrostatic_repulsion_force(sep, concentration): #e-11 wrong shape curve
    kappa = 1 / debye_length(concentration)
    # Force per unit area
    force_per_area = ((epsilon_r * epsilon_0 * Psi_0**2 * kappa) / 2) * np.exp(-kappa * sep)
    # Approximating the contact area as a circle with radius equal to the separation distance
    contact_area = 2 * np.pi * radius
    # Total force
    total_force = force_per_area * contact_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force
    # WHY ARE YOU FLAAAAAT

# Updated function to calculate electrostatic repulsion force between a sphere and a plane
def electrostatic_repulsion_force_e34(sep, concentration): # e-34
    kappa = 1 / debye_length(concentration)
    # Force per unit area
    force_per_area = 2 * (epsilon_r * epsilon_0 * Psi_0**2)**2 * np.exp(-kappa * sep) # IMaSD Ch14 14.56
    # Approximating the contact area as a circle with radius equal to the separation distance
    contact_area = 2 * np.pi * radius
    # Total force
    total_force = force_per_area * contact_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

def electrostatic_repulsion_force_e20(sep, concentration): # e-20 #, Alex's attempt
    kappa = 1 / debye_length(concentration)
    # Force per unit area
    force_per_area = epsilon_r * epsilon_0 * Psi_0**2 * np.exp(-kappa * sep)
    # Derjaguin approximation to find the force
    total_force = 2 * np.pi * radius * force_per_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

# Function to calculate electrostatic repulsion force between a sphere and a plane
def electrostatic_repulsion_force_weird(sep, concentration): # e-12 but backwards??
    kappa = 1 / debye_length(concentration)
    sigma = epsilon_0 * epsilon_r * kappa * Psi_0
    # Derjaguin approximation to find the force
    force = 2 * np.pi * radius * epsilon_0 * epsilon_r * sigma**2 * ((kappa * np.exp(-kappa * sep)) / (1 + np.exp(-kappa * sep)))
    print(f"Electrostatic Force for distance h {sep}: " + str(force))
    return force

def electrostatic_repulsion_force_agh(sep, concentration): # e-10 but also wrong (u shaped)
    kappa = 1 / debye_length(concentration)
    # Force per unit area
    force_per_area = (2 * np.pi * radius * (epsilon_r * epsilon_0 * kappa * Psi_0)**2 * np.exp(-kappa * sep)) / (kappa * epsilon_0 * epsilon_r)
    # Derjaguin approximation to find the force
    total_force = force_per_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force

def electrostatic_repulsion_force_bad(sep, concentration): # like e20 or something idk
    kappa = 1 / debye_length(concentration)
    # Force per unit area
    force_per_area = 4 * np.pi ** 2 * radius * epsilon_r * epsilon_0 * Psi_0**2 * np.exp(-kappa * sep)
    # Derjaguin approximation to find the force
    total_force = force_per_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force


# These two seem closest

def electrostatic_repulsion_force(sep, concentration): #Right numbers wrong shape so close :(
    kappa = 1 / debye_length(concentration)
    #Calculate the electrostatic force per unit area using the Derjaguin approximation for a sphere-plane interaction
    
    Z = 64 * np.pi * epsilon_0 * epsilon_r * (kB * temperature / e)**2 * np.tanh(1 * e * Psi_0 / (4 * kB * temperature))**2
    #force_per_area = kappa * radius * Z * np.exp(-kappa * sep)
    W_flat = (kappa / 2 * np.pi) * Z * np.exp(-kappa * sep)
    force_per_area = 2 * np.pi * radius * W_flat
    # Derjaguin approximation to find the force
    total_force = force_per_area
    print(f"Electrostatic Force for distance h {sep}: " + str(total_force))
    return total_force