import numpy as np
import matplotlib.pyplot as plt
import math

# Constants
epsilon_r = 62  # Relative permittivity for 50:50 water-glycerol mixture # From that other peper
epsilon_0 = 8.854e-12 # Vacuum permittivity in F/m.
radius = 3.3e-6  # Radius of the sphere in meters (2 um diameter)
Psi_0 = -0.040  # Surface potential in volts (25 mV as an approximation for neutral pH) # -0.04 taken from https://nanocomposix.com/pages/silica-physical-properties
temperature = 300  # Room temperature in Kelvin (25°C)
NA = 6.022e23  # Avogadro's number
e = 1.602e-19  # Elementary charge in Coulombs
kB = 1.38064852e-23  # Boltzmann constant in m^2 kg s^-2 K^-1
h = 2e-9  # Separation distance in meters
concentration = 1.6  # mM
# Assuming a typical Hamaker constant for glass-water interaction
Hamaker_constant = 0.63e-20  # in Joules # From psi_0 values silica.pdf
valency = 1 # aka little z

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


# Function to calculate electrostatic repulsion energy with LiCl concentration
#def electrostatic_repulsion_energ(sep, concentration):
#    """Calculate the electrostatic repulsion energy per unit area with LiCl concentration."""
#    kappa = 1 / debye_length(concentration)
#    force_electro = ((2 * 3.14 * radius * epsilon_r * epsilon_0 * Psi_0**2)/ kappa) * np.exp(-kappa * sep)
#    #force_electro = 64 * (concentration * 1e-3) * kB * T *
#    print("Electrostatic Force: " + str(force_electro))
#    return force_electro

# Function to calculate electrostatic repulsion energy with LiCl concentration
#def electrostatic_repulsion_force(sep, concentration):
#    kappa = 1 / debye_length(concentration)
#    repulsion = (epsilon_r * epsilon_0 * Psi_0**2) * circle_geom * np.exp(-kappa * sep)
#    print(f"Electrostatic Force for distance h {sep}: " + str(repulsion))
#    return repulsion


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

def electrostatic_repulsion_force(sep, concentration): #sphere plane fig 14.10 IMaSD
    kappa = 1 / debye_length(concentration)
    Z = 64 * np.pi * epsilon_0 * epsilon_r * (kB * temperature / e)**2 * np.tanh(valency * e * Psi_0 / (4 * kB * temperature))**2
    force_per_area = kappa * radius * Z * np.exp(-kappa * sep)
    print(f"Electrostatic Force for distance h {sep}: " + str(force_per_area))
    return force_per_area

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
# Assuming van der Waals interaction does not change with LiCl concentration for simplicity

# Calculation for a range of LiCl concentrations
#li_cl_concentrations = np.linspace(0.6, 550, 100)  # mM
#li_cl_concentrations = [0.6, 1.6, 5, 10, 25, 50, 230, 550]
li_cl_concentrations = np.arange(6e-1, 5.5e2, 1e-1)  # mM

# Placeholder function for van der Waals energy
def van_der_Waals_force(sep): # Sphere - plane fig 13.1 IMaSD
    force_vdw = - ((Hamaker_constant * radius) / (6 * sep ** 2)) # IMaSF ch13 eq 13.11b
    print(f"van der Waals Force for h {sep}: " + str(force_vdw))
    return force_vdw
    #return -Hamaker_constant / (12 * np.pi * h**2)

def van_der_Waals_force_plane_plane(sep): # plane - plane fig 13.1 IMaSD
    force_vdw = -Hamaker_constant / (6 * np.pi * sep**3)
    print(f"van der Waals Force for h {sep}: " + str(force_vdw))
    return force_vdw

# Total interaction energy function including LiCl concentration
def total_interaction_energy(sep, concentration):
    return van_der_Waals_force(sep) + electrostatic_repulsion_force(sep, concentration)


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
    force = total_interaction_energy(h, conc)
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
