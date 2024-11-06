import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import sys

Hybrid_data = np.load("Data/Hybrid_data.npz")
C_grid = Hybrid_data["C_grid"]
D_grid = Hybrid_data["D_grid"]
combined_grid = Hybrid_data["combined_grid"]
SSA_X = Hybrid_data["SSA_X"]
PDE_X = Hybrid_data["PDE_X"]
time_vector = Hybrid_data["time_vector"]

parameters = json.load(open("data/parameters.json"))

h = parameters["h"]
bar_positions = SSA_X   # Shift left to center on the interval

SSA_data = np.load("Data/Pure_SSA_data.npz")
SSA_grid = SSA_data["SSA_grid"]

Pure_PDE = np.load("Data/PDE_data.npz")
pure_PDE_grid = Pure_PDE["PDE_grid"]

print(f"shape of C_grid: {C_grid.shape}")
print(f"shape of D_grid: {D_grid.shape}")
print(f"shape of combined_grid: {combined_grid.shape}")
print(f"shape of SSA_X: {SSA_X.shape}")
print(f"shape of PDE_X: {PDE_X.shape}")
print(f"shape of time_vector: {time_vector.shape}")
print(f"shape of SSA_grid: {SSA_grid.shape}")

analytic_sol = np.zeros_like(C_grid)
production_rate = parameters["production_rate"]
degradation_rate = parameters["degredation_rate"]
initial_SSA = parameters["initial_SSA"]
concentration_threshold = parameters["threshold_conc"]

initial_conc = initial_SSA[0] / h

for i in range(analytic_sol.shape[1]):
    analytic_sol[:, i] = production_rate / degradation_rate + (initial_conc - production_rate / degradation_rate) * np.exp(-degradation_rate * time_vector[i])

fig, ax = plt.subplots()

# Initial bar plot for SSA data with custom bar widths
bar_SSA = ax.bar(bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='SSA (Bar Chart)')

# Initial plot for PDE data
line_PDE, = ax.plot(PDE_X[:], C_grid[:, 0], label='PDE', color='green')
line_combined, = ax.plot(PDE_X, combined_grid[:, 0], 'k--', label='Combined')
line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic Solution', color='red')
line_pure_SSA, = ax.plot(SSA_X + h / 2, SSA_grid[:, 0] / h, label='Pure SSA', color='orange')
line_pure_PDE, = ax.plot(PDE_X, pure_PDE_grid[:, 0], 'g--', label='Pure PDE')

# Set titles and labels
ax.set_xlabel('Spatial Domain')
ax.set_ylabel('Species Concentration')
ax.set_title('SSA and PDE Data Animation')

# Set the axis limits
ax.set_xlim(0, parameters["domain_length"])  # Set based on the spatial domain [0, L]
ax.set_ylim(0, 220)  # Adjust based on the total range of values

# Add grid and legend
ax.grid(True)
ax.legend()

# Add concentration threshold line
ax.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold Concentration')

# Update function for animation
def update(frame):
    # Update SSA bar heights
    for bar, height in zip(bar_SSA, D_grid[:, frame] / h):
        bar.set_height(height)
    
    # Update PDE line
    line_PDE.set_ydata(C_grid[:, frame])
    line_combined.set_ydata(combined_grid[:, frame])
    line_analytic.set_ydata(analytic_sol[:, frame])
    line_pure_SSA.set_ydata(SSA_grid[:, frame] / h)
    line_pure_PDE.set_ydata(pure_PDE_grid[:, frame])
    
    # Return all updated artists as a tuple
    return (*bar_SSA, line_PDE, line_combined, line_analytic, line_pure_SSA, line_pure_PDE)
    #return (*bar_SSA, line_combined, line_analytic, line_pure_SSA, line_pure_PDE)

# Define interval_number
if len(sys.argv) == 1:
    interval_number = 1
elif len(sys.argv) == 2:
    interval_number = int(sys.argv[1])
else:
    print("Usage: python animate.py [interval_number]")
    sys.exit(1)

# Create animation
step_size = 4
ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), step_size), interval=interval_number)

# Show the animated plot
plt.show()

plt.figure()

analytic_av = np.mean(analytic_sol, axis=0)
PDE_av = np.mean(C_grid, axis=0)
combined_av = np.mean(combined_grid, axis=0)
SSA_av = np.mean(D_grid / h, axis=0)
pure_PDE_av = np.mean(pure_PDE_grid, axis=0)
pure_SSA_av = np.mean(SSA_grid / h, axis=0)

plt.plot(time_vector, analytic_av, label='Analytic', color='red')
plt.plot(time_vector, PDE_av, 'g', label='PDE_av')
plt.plot(time_vector, SSA_av, 'b', label='SSA_av')
plt.plot(time_vector, pure_PDE_av, 'g--', label='Pure PDE')
plt.plot(time_vector, combined_av, 'k--', label='Combined_av')
plt.plot(time_vector, pure_SSA_av, label='Pure SSA', color='orange')

# Add concentration threshold line
plt.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold Concentration')

plt.legend()
plt.xlabel('Time')
plt.ylabel('Average concentration over')
plt.show()