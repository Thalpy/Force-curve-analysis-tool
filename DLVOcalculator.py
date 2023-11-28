import numpy as np

# Constants
epsilon_r = 50  # Relative permittivity for 50:50 water-glycerol mixture
epsilon_0 = 8.854e-12 # Vacuum permittivity in F/m.
radius = 1e-6  # Radius of the sphere in meters (2 um diameter)
Psi_0 = 0.025  # Surface potential in volts (25 mV as an approximation for neutral pH)
temperature = 298.15  # Room temperature in Kelvin (25°C)
NA = 6.022e23  # Avogadro's number
e = 1.602e-19  # Elementary charge in Coulombs
kB = 1.38064852e-23  # Boltzmann constant in m^2 kg s^-2 K^-1

# Assuming a typical Hamaker constant for glass-water interaction
Hamaker_constant = 2.5e-20  # in Joules

# Function to calculate Debye length (1/kappa) based on LiCl concentrationa
def debye_length(concentration, valence=1, temperature=298.15):
    # Convert concentration from mM to mol/m^3
    concentration_mol = concentration * 1e-3
    # Ionic strength
    I = 0.5 * valence**2 * concentration_mol
    # Debye length (1/kappa)
    return np.sqrt(epsilon_0 * epsilon_r * kB * temperature / (2 * NA * e**2 * I))


# Function to calculate electrostatic repulsion energy with LiCl concentration
def electrostatic_repulsion_energy(h, concentration):
    """Calculate the electrostatic repulsion energy per unit area with LiCl concentration."""
    kappa = 1 / debye_length(concentration)
    return (epsilon_r * epsilon_0 * Psi_0**2 / kappa) * np.exp(-kappa * h)

# Assuming van der Waals interaction does not change with LiCl concentration for simplicity
# Otherwise, you would need to provide an equation that accounts for this

# Calculation for a range of LiCl concentrations
li_cl_concentrations = np.linspace(0.6, 550, 100)  # mM

# Placeholder function for van der Waals energy
def van_der_Waals_energy(h):
    return -Hamaker_constant / (12 * np.pi * h**2)

# Total interaction energy function including LiCl concentration
def total_interaction_energy(h, concentration):
    return van_der_Waals_energy(h) + electrostatic_repulsion_energy(h, concentration)

# Derjaguin approximation force calculation function
def derjaguin_approximation_force(h, radius, concentration):
    W_h = total_interaction_energy(h, concentration)
    return 2 * np.pi * radius * W_h

# Calculation of forces for each LiCl concentration
forces_at_contact = [derjaguin_approximation_force(1e-8, radius, conc) for conc in li_cl_concentrations]

# Output the results
for conc, force in zip(li_cl_concentrations, forces_at_contact):
    print(f"LiCl concentration: {conc} mM - Force at contact: {force} N")
