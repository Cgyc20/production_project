import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import sys

C_grid = np.load("Data/C_grid.npy")
D_grid = np.load("Data/D_grid.npy")
combined_grid = np.load("Data/combined_grid.npy")
SSA_X = np.load("Data/SSA_X.npy")
PDE_X = np.load("Data/PDE_X.npy")
time_vector = np.load("Data/time_vector.npy")
parameters = json.load(open("Data/parameters.json"))
h = parameters["h"]
bar_positions = SSA_X   # Shift left to center on the interval

SSA_grid = np.load("Data/Pure_SSA_grid.npy")

print(f"shape of C_grid: {C_grid.shape}")
print(f"shape of D_grid: {D_grid.shape}")
print(f"shape of combined_grid: {combined_grid.shape}")
print(f"shape of SSA_X: {SSA_X.shape}")
print(f"shape of PDE_X: {PDE_X.shape}")
print(f"shape of time_vector: {time_vector.shape}")
print(f"shape of SSA_grid: {SSA_grid.shape}")

analytic_sol = np.zeros_like(C_grid)
production_rate = parameters["production_rate"]
degredation_rate = parameters["degredation_rate"]
print(parameters)
initial_SSA = parameters["initial_SSA"]
initial_conc = initial_SSA[0]/h


for i in range(analytic_sol.shape[1]):
    analytic_sol[:,i] = production_rate/degredation_rate+(initial_conc-production_rate/degredation_rate)*np.exp(-degredation_rate*time_vector[i])

fig, ax = plt.subplots()

# Initial bar plot for SSA data with custom bar widths
bar_SSA = ax.bar(bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='SSA (Bar Chart)')

# Initial plot for PDE data
line_PDE, = ax.plot(PDE_X[:], C_grid[:, 0], label='PDE', color='red')
line_combined, = ax.plot(PDE_X , combined_grid[:, 0] , label='Combined', color='green')
line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic Solution', color='purple')
line_pure_SSA, = ax.plot(SSA_X + h / 2, SSA_grid[:, 0], label='Pure SSA', color='orange')

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

# Update function for animation
def update(frame):
    # Update SSA bar heights
    for bar, height in zip(bar_SSA, D_grid[:, frame]/h):
        bar.set_height(height)
    
    # Update PDE line
    line_PDE.set_ydata(C_grid[:, frame])
    line_combined.set_ydata(combined_grid[:, frame])
    line_analytic.set_ydata(analytic_sol[:, frame])
    line_pure_SSA.set_ydata(SSA_grid[:, frame])
    
    # Return all updated artists as a tuple
    return (*bar_SSA, line_PDE, line_combined, line_analytic, line_pure_SSA)

# Define interval_number
if len(sys.argv) == 1:
    interval_number = 1
elif len(sys.argv) == 2:
    interval_number = int(sys.argv[1])
else:
    print("Usage: python animate.py [interval_number]")
    sys.exit(1)

# Create animation
ani = FuncAnimation(fig, update, frames=range(len(time_vector)), interval=interval_number, blit=True)

# Show the animated plot
plt.show()


plt.figure()

analytic_av = np.mean(analytic_sol,axis=0)
PDE_av = np.mean(C_grid,axis=0)
combined_av = np.mean(combined_grid,axis=0)
SSA_av = np.mean(D_grid/h,axis=0)
pure_SSA_av = np.mean(SSA_grid/h,axis=0)

plt.plot(time_vector,analytic_av,label = 'Analytic')
plt.plot(time_vector,PDE_av, label = 'PDE_av')
plt.plot(time_vector,combined_av, label = 'Combined_av')
plt.plot(time_vector,SSA_av, label = 'SSA')

plt.legend()
plt.xlabel('Time')
plt.ylabel('Average concentration over')
plt.show()