import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from production_project import PDE
from scipy.integrate import quad

# Parameters
domain_length = 1
PDE_points = 100
total_time = 100
timestep = 0.01
production_rate = 10
degradation_rate = 0.1
diffusion_rate = 10e-3
K_1_K_2 = production_rate / degradation_rate

# Initial conditions
PDE_domain_1 = np.linspace(0, domain_length, PDE_points)
PDE_initial = 10 * np.sin(5 * PDE_domain_1 * np.pi / domain_length) + 10

# Define the integrand for A_0
def integrand_A0(x):
    return 10 * np.sin(5 * np.pi * x) + 10

# Define the integrand for A_k (k > 0)
def integrand(x, k):
    return 20 * (np.sin(5 * x*np.pi)-K_1_K_2) * np.cos(k * np.pi * x)

# Compute A_0 separately
A_0, error_A0 = quad(integrand_A0, 0, 1)
A_0-= K_1_K_2
print(f"A_0 = {A_0}")
n_max = 1001
# Compute A_k for k = 1 to 100
A_k_values = [A_0]  # Start with A_0
for k in range(1, n_max):
    A_k, error_Ak = quad(integrand, 0, 1, args=(k,))
    A_k_values.append(A_k)  # Append A_k to the list

# Optionally, display some of the results
for k, A_k in enumerate(A_k_values[:10]):  # Display first 10 coefficients as an example
    print(f"A_{k} = {A_k}")

# Define the analytic solution function
def u(x, t):
    # Compute the sum for each n term up to the length of A_k_values
    series_sum = np.zeros_like(x, dtype=np.float64)  # Initialize sum as zeros for each x
    for n in range(len(A_k_values)):
        A_n = A_k_values[n]
        decay_factor = np.exp(-(degradation_rate + (n**2 * np.pi**2 * diffusion_rate) / domain_length**2) * t)
        series_sum += A_n * decay_factor * np.cos(n * np.pi * x / domain_length)
    
    # Add the steady-state component k1/k2
    return series_sum + K_1_K_2

# Compute u(x, t) for each t in time_vector and store the results
time_vector = np.arange(0, total_time, timestep)
results = np.array([u(PDE_domain_1, t) for t in time_vector])

# Create the PDE model and run the simulation
Model = PDE(domain_length, PDE_points, total_time, timestep, production_rate, degradation_rate, diffusion_rate, PDE_initial)
PDE_grid = Model.run_simulation()
PDE_domain = Model.PDE_X

# Set up the figure, axis, and plot element
fig, ax = plt.subplots()
line_PDE, = ax.plot(PDE_domain[:], PDE_grid[:, 0], label='PDE', color='green')
line_analytic, = ax.plot(PDE_domain_1, results[0, :], label='Analytic Solution', color='red')

# Set titles and labels
ax.set_xlabel('Spatial Domain')
ax.set_ylabel('Species Concentration')
ax.set_title('PDE and Analytic Solution Animation')

# Set the axis limits
ax.set_xlim(0, domain_length)  # Set based on the spatial domain [0, L]
ax.set_ylim(0, 220)  # Adjust based on the total range of values

# Add grid and legend
ax.grid(True)
ax.legend()

# Update function for animation
def update(frame):
    # Update PDE line
    line_PDE.set_ydata(PDE_grid[:, frame])
    # Update analytic solution line
    line_analytic.set_ydata(results[frame, :])
    return line_PDE, line_analytic

# Create animation
step_size = 1
ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), step_size), interval=20, blit=True)

# Show the animated plot
plt.show()