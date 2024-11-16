# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import sys

# Load data from .npz files
Hybrid_data = np.load("Data/Hybrid_data.npz")
C_grid = Hybrid_data["C_grid"]
D_grid = Hybrid_data["D_grid"]
combined_grid = Hybrid_data["combined_grid"]
SSA_X = Hybrid_data["SSA_X"]
PDE_X = Hybrid_data["PDE_X"]
time_vector = Hybrid_data["time_vector"]

# Load simulation parameters from JSON file
parameters = json.load(open("data/parameters.json"))
h = parameters["h"]
deltax = parameters["deltax"]
bar_positions = SSA_X  # Adjust bar positions for the SSA histogram display

# Load additional data sets
SSA_data = np.load("Data/Pure_SSA_data.npz")
SSA_grid = SSA_data["SSA_grid"]
Pure_PDE = np.load("Data/PDE_data.npz")
pure_PDE_grid = Pure_PDE["PDE_grid"]

# Display shapes of loaded data arrays for verification
print(f"shape of C_grid: {C_grid.shape}")
print(f"shape of D_grid: {D_grid.shape}")
print(f"shape of combined_grid: {combined_grid.shape}")
print(f"shape of SSA_X: {SSA_X.shape}")
print(f"shape of PDE_X: {PDE_X.shape}")
print(f"shape of time_vector: {time_vector.shape}")
print(f"shape of SSA_grid: {SSA_grid.shape}")

print(f"The conmbined grid: {combined_grid[:5,0:5]}")
print(f"The PDE grid: {C_grid[:5,0:5]}")
# Initialize analytical solution array
analytic_sol = np.zeros_like(C_grid)

# Retrieve parameters for analytical solution
production_rate = parameters["production_rate"]
degradation_rate = parameters["degredation_rate"]
initial_SSA = parameters["initial_SSA"]
concentration_threshold = parameters["threshold_conc"]

# Initial concentration for the analytical solution
initial_conc = initial_SSA[0] / h

# Calculate the analytical solution at each time point
for i in range(analytic_sol.shape[1]):
    analytic_sol[:, i] = (production_rate / degradation_rate + 
                          (initial_conc - production_rate / degradation_rate) * 
                          np.exp(-degradation_rate * time_vector[i]))

# Calculate total mass for SSA data
SSA_total_mass = np.sum(D_grid, axis=0)
pure_SSA_Mass = np.sum(SSA_grid, axis=0)

# Function to calculate total mass for continuous data
def calculate_mass_continuous(data_grid):
    total_mass = np.zeros(data_grid.shape[1])
    # Calculate total mass over continuous functions using the left-hand rule
    for i in range(data_grid.shape[1]):
        current_sum = 0
        for j in range(data_grid.shape[0] - 1):
            current_sum += data_grid[j, i] * deltax

        total_mass[i] = current_sum
        if i<20:
            print(f"Total mass at time {i}: {total_mass[i]}")
    return total_mass

# Calculate total mass over time for each grid
analytic_total_mass = calculate_mass_continuous(analytic_sol)
PDE_total_mass = calculate_mass_continuous(C_grid)
combined_total_mass = calculate_mass_continuous(combined_grid)
pure_PDE_total_mass = calculate_mass_continuous(pure_PDE_grid)


print(f"The combineed solution at first timestep: {combined_total_mass[0]}")
print(f"The PDE total mass at first timestep: {pure_PDE_total_mass[0]}")
# Adjust time vector for relative error calculation
adjusted_time_vector = time_vector[1:]

# Calculate relative error between the combined and analytical solutions
hybrid_relative_error = (combined_grid[:, 1:] - analytic_sol[:, 1:]) / np.abs(analytic_sol[:, 1:])
average_relative_error = np.mean(hybrid_relative_error, axis=0)

# Define interval for animation if provided in command-line arguments
if len(sys.argv) == 1:
    interval_number = 1
elif len(sys.argv) == 2:
    interval_number = int(sys.argv[1])
else:
    print("Usage: python animate.py [interval_number]")
    sys.exit(1)

# Plotting and Animation

# Set up the initial figure and axis for animation
fig, ax = plt.subplots()

# Initial bar plot for SSA data with custom bar widths
bar_SSA = ax.bar(bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='SSA (Bar Chart)')

# Initial plot for PDE and other data
line_PDE, = ax.plot(PDE_X, C_grid[:, 0], 'g--',label='PDE')
line_combined, = ax.plot(PDE_X, combined_grid[:, 0], 'k--', label='Combined')
line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic Solution', color='red')
line_pure_SSA, = ax.plot(SSA_X + h / 2, SSA_grid[:, 0] / h,'m', label='Pure SSA')
line_pure_PDE, = ax.plot(PDE_X, pure_PDE_grid[:, 0], 'g', label='Pure PDE')

# Set titles, labels, and axis limits
ax.set_xlabel('Spatial Domain')
ax.set_ylabel('Species Concentration')
ax.set_title('SSA and PDE Data Animation')
ax.set_xlim(0, parameters["domain_length"])
ax.set_ylim(0, 220)
ax.grid(True)
ax.legend()

# Plot concentration threshold line
ax.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold Concentration')

# Update function for animation frames
def update(frame):
    # Update SSA bar heights
    for bar, height in zip(bar_SSA, D_grid[:, frame] / h):
        bar.set_height(height)
    
    # Update lines for each dataset
    line_PDE.set_ydata(C_grid[:, frame])
    line_combined.set_ydata(combined_grid[:, frame])
    line_analytic.set_ydata(analytic_sol[:, frame])
    line_pure_SSA.set_ydata(SSA_grid[:, frame] / h)
    line_pure_PDE.set_ydata(pure_PDE_grid[:, frame])
    
    return (*bar_SSA, line_PDE, line_combined, line_analytic, line_pure_SSA, line_pure_PDE)

# Create animation
step_size = 1
ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), step_size), interval=interval_number)

# Display the animated plot
plt.show()

print(f"The combined mass: {combined_total_mass[:10]}")
# Plot total mass over time for each model
plt.figure()
plt.plot(time_vector, SSA_total_mass, 'b--',label='SSA part')
plt.plot(time_vector, analytic_total_mass,'r', label='Analytic')
plt.plot(time_vector, PDE_total_mass, 'g--',label='PDE part')
plt.plot(time_vector, combined_total_mass,'k--',label='Combined')
plt.plot(time_vector,pure_SSA_Mass, 'm',label = 'Pure SSA' )
plt.plot(time_vector, pure_PDE_total_mass,'g', label='Pure PDE')
plt.plot(time_vector,np.ones_like(pure_PDE_total_mass)*concentration_threshold,'k--',linewidth = 0.5,label = 'Threshold')
plt.xlabel('Time')
plt.ylabel('Total Mass')
plt.title('Total Mass over Time')
plt.legend()
plt.grid(True)
plt.show()

# Plot average relative error over time
plt.figure()
plt.plot(adjusted_time_vector, average_relative_error, label='Relative Error')
plt.xlabel('Time')
plt.ylabel('Average Relative Error')
plt.title('Relative Error of Hybrid Model')
plt.legend()
plt.grid(True)
plt.show()
