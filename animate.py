import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import seaborn as sns

def main():
    # Set seaborn style
    sns.set_theme(style="whitegrid")

    # Load data from .npz files
    Hybrid_data = np.load("Data/Hybrid_data.npz")
    C_grid = Hybrid_data["PDE_grid"]
    D_grid = Hybrid_data["SSA_grid"]
    combined_grid = Hybrid_data["combined_grid"]
    SSA_X = Hybrid_data["SSA_X"]
    PDE_X = Hybrid_data["PDE_X"]
    time_vector = Hybrid_data["time_vector"]

    SSA_data = np.load("Data/Pure_SSA_data.npz")
    SSA_grid = SSA_data["SSA_grid"]

    PDE_data = np.load("Data/PDE_data.npz")
    PDE_grid = PDE_data["PDE_grid"]

    # Load simulation parameters from JSON file
    parameters = json.load(open("data/parameters.json"))
    h = parameters["h"]
    deltax = parameters["deltax"]
    bar_positions = SSA_X

    # Initialize analytical solution array
    analytic_sol = np.zeros_like(C_grid)

    # Retrieve parameters for analytical solution
    production_rate = parameters["production_rate"]
    degradation_rate = parameters["degradation_rate"]
    initial_SSA = parameters["initial_SSA"]
    concentration_threshold = parameters["threshold_conc"]
    domain_length = parameters["domain_length"]

    # Calculate analytical solution
    initial_conc = initial_SSA[0] / h
    for i in range(analytic_sol.shape[1]):
        analytic_sol[:, i] = (
            production_rate / degradation_rate
            + (initial_conc - production_rate / degradation_rate) * np.exp(-degradation_rate * time_vector[i])
        )

    # Function to calculate total mass for continuous data
    def calculate_mass_continuous(data_grid, deltax):
        return np.sum(data_grid, axis=0) * deltax
    
    def calculate_mass_discrete(data_grid):
        return np.sum(data_grid,axis=0)

    # Calculate total mass for all solutions
    analytic_total_mass = calculate_mass_continuous(analytic_sol, deltax)
    Hybrid_PDE_total_mass = calculate_mass_continuous(C_grid, deltax)
    pure_PDE_total_mass = calculate_mass_continuous(PDE_grid, deltax)
    combined_total_mass = calculate_mass_continuous(combined_grid, deltax)

    # Calculate total mass for pure SSA and pure PDE
    SSA_total_mass = calculate_mass_discrete(SSA_grid)
    Hybrid_SSA_mass = calculate_mass_discrete(D_grid)

    # Calculate relative error for combined solution
    relative_error_combined = (combined_total_mass - analytic_total_mass) / analytic_total_mass

    # Plotting and Animation
    fig, ax = plt.subplots(figsize=(10, 6))

    # Initial SSA bar plot
    bar_SSA = ax.bar(
        bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='SSA (Bar Chart)', alpha=0.7
    )

    # Continuous plots
    line_PDE, = ax.plot(PDE_X, C_grid[:, 0], 'g', label='PDE', linewidth=2)
    line_combined, = ax.plot(PDE_X, combined_grid[:, 0], 'k--', label='Combined', linewidth=2)
    line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic', color='red', linewidth=2)

    # Threshold line
    threshold_line = ax.axhline(y=concentration_threshold, color='purple', linestyle='--', label='Threshold', linewidth=1.5)

    # Axis labels and title
    ax.set_xlabel('Spatial Domain', fontsize=12)
    ax.set_ylabel('Species Concentration', fontsize=12)
    ax.set_title('Hybrid simulation', fontsize=14)
    ax.set_xlim(0, domain_length)
    ax.set_ylim(0, max(np.max(combined_grid) * 1.1, concentration_threshold * 1.1))
    ax.grid(True, linestyle='--', alpha=0.6)

    # Add a text annotation for the timestamp
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')

    # Steady-state concentration
    steady_state_concentration = production_rate / degradation_rate

    # Adjust y-axis limit to ensure steady state is included
    y_max = max(np.max(combined_grid) * 1.1, steady_state_concentration * 1.1, concentration_threshold * 1.1)
    ax.set_ylim(0, y_max)

    # Add a steady-state line
    steady_state_line = ax.axhline(
        y=steady_state_concentration,
        color='gray',
        linestyle='--',
        label='Steady State',
        linewidth=1.5,
    )

    # Update function for animation
    def update(frame):
        for bar, height in zip(bar_SSA, D_grid[:, frame] / h):
            bar.set_height(height)
        line_combined.set_ydata(combined_grid[:, frame])
        line_PDE.set_ydata(C_grid[:, frame])
        line_analytic.set_ydata(analytic_sol[:, frame])
        
        # Update the timestamp
        time_text.set_text(f'Time: {time_vector[frame]:.2f}')
        
        return (*bar_SSA, line_combined, line_PDE, line_analytic, time_text, threshold_line, steady_state_line)

    # Create animation
    ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), 1), interval=10)

    # Set legend position fixed
    # Adjust the figure layout to make space for the legend
    fig.subplots_adjust(right=0.8)

    # Set legend position outside the main pane
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)

    # Display the animation
    plt.show()

    # Additional plot: Total mass over time
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

    # Plot relative error for combined solution
    plt.figure(figsize=(12, 6))
    plt.plot(time_vector, relative_error_combined, 'k--', label='Relative Error (Combined)', linewidth=2)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Relative Error', fontsize=12)
    plt.title('Relative Error of Combined Solution over Time', fontsize=14)
    plt.legend(fontsize=10)
    plt.ylim(-0.02,0.04)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

if __name__ == "__main__":
    main()