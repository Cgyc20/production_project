import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import seaborn as sns

def load_data():
    """Load data from .npz files and JSON parameters."""
    Hybrid_data = np.load("Data/Hybrid_data.npz")
    SSA_data = np.load("Data/Pure_SSA_data.npz")
    PDE_data = np.load("Data/PDE_data.npz")
    parameters = json.load(open("data/parameters.json"))
    return Hybrid_data, SSA_data, PDE_data, parameters

def initialize_analytic_solution(C_grid):
    """Initialize the analytical solution array."""
    return np.zeros_like(C_grid)

def determine_coefficients(n, a, b, L, initial_conc):
    """Calculate coefficients for the analytical solution."""
    coef_list = [initial_conc*(b-a)/L]
    for i in range(1, n):
        coefficient = initial_conc * ((2/(i*np.pi*L)) * (np.sin(i*np.pi*b/L) - np.sin(i*np.pi*a/L)))
        coef_list.append(coefficient)
    return coef_list

def determine_coefficients(n, a, b, L, initial_conc):
    """Calculate coefficients for the analytical solution using vectorized operations."""
    i = np.arange(1, n)
    coefficients = initial_conc * ((2 / (i * np.pi * L)) * (np.sin(i * np.pi * b / L) - np.sin(i * np.pi * a / L)))
    return np.concatenate(([(b - a) / L], coefficients))

def calculate_analytic_solution(analytic_sol, coefficients, a0, PDE_X, domain_length, diffusion_rate, time_vector):
    """Calculate the analytical solution over time."""
    for i in range(analytic_sol.shape[1]):
        t = time_vector[i]
        for j in range(analytic_sol.shape[0]):
            x = PDE_X[j]
            u_xt = a0
            for n, coef in enumerate(coefficients):
                u_xt += coef * np.cos(n * np.pi * x / domain_length) * np.exp(-diffusion_rate * (n * np.pi / domain_length) ** 2 * t)
            analytic_sol[j, i] = u_xt
    return analytic_sol

def calculate_mass_continuous(data_grid, deltax):
    """Calculate total mass for continuous data."""
    return np.sum(data_grid, axis=0) * deltax

def calculate_mass_discrete(data_grid):
    """Calculate total mass for discrete data."""
    return np.sum(data_grid, axis=0)

def plot_initial_setup(ax, bar_positions, D_grid, h, PDE_X, C_grid, combined_grid, analytic_sol, concentration_threshold, domain_length):
    """Setup the initial plot and return plot elements."""
    bar_SSA = ax.bar(bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='SSA (Bar Chart)', alpha=0.7)
    line_PDE, = ax.plot(PDE_X, C_grid[:, 0], 'g', label='PDE', linewidth=2)
    line_combined, = ax.plot(PDE_X, combined_grid[:, 0], 'k--', label='Combined', linewidth=2)
    line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic', color='red', linewidth=2)
    threshold_line = ax.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold', linewidth=1.5)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')
    return bar_SSA, line_PDE, line_combined, line_analytic, threshold_line, time_text

def update(frame, bar_SSA, D_grid, h, line_combined, combined_grid, line_PDE, C_grid, line_analytic, analytic_sol, time_text, time_vector):
    """Update function for animation."""
    for bar, height in zip(bar_SSA, D_grid[:, frame] / h):
        bar.set_height(height)
    line_combined.set_ydata(combined_grid[:, frame])
    line_PDE.set_ydata(C_grid[:, frame])
    line_analytic.set_ydata(analytic_sol[:, frame])
    time_text.set_text(f'Time: {time_vector[frame]:.2f}')
    return (*bar_SSA, line_combined, line_PDE, line_analytic, time_text)

def plot_total_mass(time_vector, combined_total_mass, Hybrid_PDE_total_mass, pure_PDE_total_mass, Hybrid_SSA_mass, SSA_total_mass, production_rate, degradation_rate, concentration_threshold):
    """Plot total mass over time."""
    plt.figure(figsize=(12, 6))
    plt.plot(time_vector, combined_total_mass, 'k--', label='Combined (Dashed)', linewidth=2)
    plt.plot(time_vector, Hybrid_PDE_total_mass, 'g--', label='Hybrid PDE', linewidth=2)
    plt.plot(time_vector, Hybrid_SSA_mass, 'b--', label='Hybrid SSA', linewidth=2)
    plt.plot(time_vector, pure_PDE_total_mass, 'g', label='Pure PDE', linewidth=2)
    plt.plot(time_vector, SSA_total_mass, 'b', label='Pure SSA', linewidth=2)
    plt.axhline(y=production_rate / degradation_rate, color='gray', linestyle='--', label='Steady State', linewidth=1.5)
    plt.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold', linewidth=1.5)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Total Mass', fontsize=12)
    plt.title('Total Mass over Time', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

def main():
    sns.set_theme(style="whitegrid")
    Hybrid_data, SSA_data, PDE_data, parameters = load_data()
    
    C_grid = Hybrid_data["PDE_grid"]
    D_grid = Hybrid_data["SSA_grid"]
    combined_grid = Hybrid_data["combined_grid"]
    SSA_X = Hybrid_data["SSA_X"]
    PDE_X = Hybrid_data["PDE_X"]
    time_vector = Hybrid_data["time_vector"]
    SSA_grid = SSA_data["SSA_grid"]
    PDE_grid = PDE_data["PDE_grid"]
    
    h = parameters["h"]
    deltax = parameters["deltax"]
    bar_positions = SSA_X
    production_rate = parameters["production_rate"]
    degradation_rate = parameters["degradation_rate"]
    concentration_threshold = parameters["threshold_conc"]
    domain_length = parameters["domain_length"]
    diffusion_rate = parameters["diffusion_rate"]
    
    analytic_sol = initialize_analytic_solution(C_grid)
    initial_conc = 40 / h
    b = 0.55
    a = 0.45
    a0 = initial_conc * (b - a) / domain_length
    coefficients = determine_coefficients(50, a, b, domain_length, initial_conc)
    analytic_sol = calculate_analytic_solution(analytic_sol, coefficients, a0, PDE_X, domain_length, diffusion_rate, time_vector)
    
    analytic_total_mass = calculate_mass_continuous(analytic_sol, deltax)
    Hybrid_PDE_total_mass = calculate_mass_continuous(C_grid, deltax)
    pure_PDE_total_mass = calculate_mass_continuous(PDE_grid, deltax)
    combined_total_mass = calculate_mass_continuous(combined_grid, deltax)
    SSA_total_mass = calculate_mass_discrete(SSA_grid)
    Hybrid_SSA_mass = calculate_mass_discrete(D_grid)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_SSA, line_PDE, line_combined, line_analytic, threshold_line, time_text = plot_initial_setup(
        ax, bar_positions, D_grid, h, PDE_X, C_grid, combined_grid, analytic_sol, concentration_threshold, domain_length
    )
    
    steady_state_concentration = production_rate / degradation_rate
    y_max = max(np.max(combined_grid) * 1.1, steady_state_concentration * 1.1, concentration_threshold * 1.1)
    ax.set_ylim(0, y_max)
    steady_state_line = ax.axhline(y=steady_state_concentration, color='gray', linestyle='--', label='Steady State', linewidth=1.5)
    
    ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), 1), interval=40, fargs=(bar_SSA, D_grid, h, line_combined, combined_grid, line_PDE, C_grid, line_analytic, analytic_sol, time_text, time_vector))
    
    fig.subplots_adjust(right=0.8)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
    plt.show()
    
    plot_total_mass(time_vector, combined_total_mass, Hybrid_PDE_total_mass, pure_PDE_total_mass, Hybrid_SSA_mass, SSA_total_mass, production_rate, degradation_rate, concentration_threshold)

if __name__ == "__main__":
    main()