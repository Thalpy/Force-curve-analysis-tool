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