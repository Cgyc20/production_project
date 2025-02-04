import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import seaborn as sns
import pandas as pd

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

    # Load SSA events and PDE update times separately
    SSA_events = np.load("Data/SSA_events_logs.npy", allow_pickle=True)
    PDE_update_times = np.load("Data/PDE_update_times.npy", allow_pickle=True)

    SSA_data = np.load("Data/Pure_SSA_data.npz")
    SSA_grid = SSA_data["SSA_grid"]

    PDE_data = np.load("Data/PDE_data.npz")
    PDE_grid = PDE_data["PDE_grid"]

    # Load simulation parameters from JSON file
    parameters = json.load(open("data/parameters.json"))
    h = parameters["h"]
    print(f"the h value in animation: {h}")
    deltax = parameters["deltax"]
    diffusion_rate = parameters["diffusion_rate"]
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
    relative_error_combined = np.abs((combined_total_mass - analytic_total_mass) / analytic_total_mass)
    relative_error_SSA = np.abs((SSA_total_mass - analytic_total_mass) / analytic_total_mass)

    # Plotting and Animation
    fig, ax = plt.subplots(figsize=(12, 8)) 

    # Initial SSA bar plot
    bar_SSA = ax.bar(
        bar_positions, D_grid[:, 0] / h, width=h, color='blue', align='edge', label='Hybrid SSA (Bar Chart)', alpha=0.7
    )

    # Initial pure SSA bar plot
    bar_pure_SSA = ax.bar(
        bar_positions, SSA_grid[:, 0] / h, width=h, color='cyan', align='edge', label='Pure SSA (Bar Chart)', alpha=0.5
    )

    # Continuous plots
    line_PDE, = ax.plot(PDE_X, C_grid[:, 0], 'g--', label='Hybrid PDE', linewidth=2)
    line_combined, = ax.plot(PDE_X, combined_grid[:, 0], 'k--', label='Combined', linewidth=2)
    # line_analytic, = ax.plot(PDE_X, analytic_sol[:, 0], label='Analytic', color='red', linewidth=2)
    line_pure_PDE, = ax.plot(PDE_X, PDE_grid[:, 0], 'g', label='Pure PDE', linewidth=2)

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
    ax.set_ylim(-20, y_max)

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
        for bar, height in zip(bar_pure_SSA, SSA_grid[:, frame] / h):
            bar.set_height(height)
        line_combined.set_ydata(combined_grid[:, frame])
        line_PDE.set_ydata(C_grid[:, frame])
        
        line_pure_PDE.set_ydata(PDE_grid[:, frame])
        
        # Update the timestamp
        time_text.set_text(f'Time: {time_vector[frame]:.2f}')
        
        return (*bar_SSA, *bar_pure_SSA, line_combined, line_PDE, line_pure_PDE, time_text, threshold_line, steady_state_line)

    # Create animation
    ani = FuncAnimation(fig, update, frames=range(0, len(time_vector), 1), interval=10)

    # Set legend position fixed
    # Adjust the figure layout to make space for the legend
    fig.subplots_adjust(right=0.8)

    # Set legend position outside the main pane
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)

    # Display the animation
    # plt.show()

    # Additional plot: Total mass over time
    plt.figure(figsize=(8, 6))

    plt.plot(time_vector, combined_total_mass, 'k--', label='Combined (Dashed)', linewidth=2)
    plt.plot(time_vector, Hybrid_PDE_total_mass, 'g--', label='Hybrid PDE', linewidth=2)
    plt.plot(time_vector, Hybrid_SSA_mass, 'b--', label='Hybrid SSA', linewidth=2)
    plt.plot(time_vector, pure_PDE_total_mass, 'g', label='Pure PDE', linewidth=2)
    plt.plot(time_vector, SSA_total_mass, 'b', label='Pure SSA', linewidth=2)
    plt.axhline(y=domain_length*(production_rate / degradation_rate), color='gray', linestyle='--', label='Steady State', linewidth=1.5)
    plt.axhline(y=domain_length*concentration_threshold, color='purple', linestyle='--', label='Threshold', linewidth=1.5)

    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Total Mass', fontsize=12)
    plt.title('Total Mass over Time', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    # plt.show()


    """Working out the average wave speed for each one"""

    #wavespeed_PDE = k_2/k_1 dM/dt
    
    def workout_wavespeed(total_mass_vector):
        """Works out the wave speed."""
        deriv = np.gradient(total_mass_vector, time_vector)
        wavespeed_vector = (degradation_rate / production_rate) * deriv
        return wavespeed_vector

    def moving_average(data, window_size):
        """Computes the moving average of the input data."""
        return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

    # Calculate wave speeds
    wavespeed_PDE = workout_wavespeed(pure_PDE_total_mass)
    wavespeed_SSA = workout_wavespeed(SSA_total_mass)
    wavespeed_hybrid = workout_wavespeed(combined_total_mass)


    # Apply moving average to smooth the wave speeds
    window_size = 10  # Adjust window size for desired smoothing
    wavespeed_PDE_smooth = moving_average(wavespeed_PDE, window_size)
    wavespeed_SSA_smooth = moving_average(wavespeed_SSA, window_size)
    wavespeed_hybrid_smooth = moving_average(wavespeed_hybrid, window_size)

    # Adjust time_vector for the moving average (because it's shorter after smoothing)
    time_vector_smooth = time_vector[:len(wavespeed_PDE_smooth)]

    max_wave_speed = np.ones_like(time_vector_smooth) * 2*np.sqrt(diffusion_rate*production_rate) #Theorical max

    # Create the side-by-side plots
    plt.figure(figsize=(16, 6))

    # Find the common y-axis limits
    min_y = min(wavespeed_PDE.min(), wavespeed_SSA.min(), wavespeed_hybrid.min(),
                wavespeed_PDE_smooth.min(), wavespeed_SSA_smooth.min(), wavespeed_hybrid_smooth.min())
    max_y = max(wavespeed_PDE.max(), wavespeed_SSA.max(), wavespeed_hybrid.max(),
                wavespeed_PDE_smooth.max(), wavespeed_SSA_smooth.max(), wavespeed_hybrid_smooth.max())


    #plot 0
    print(f"the shape of the D -grid is {D_grid.shape}")
    print(f"The shape of the timevector is {time_vector.shape}")

    plt.plot(time_vector,D_grid[0,:],'b--',label = 'Hybrid SSA')
    plt.plot(time_vector,SSA_grid[0,:],'g',label = 'Pure SSA')
    #Now plotting the per compartment threshold
    plt.plot(time_vector, np.ones_like(time_vector)*concentration_threshold*h, label = 'Compartment Threshold')
    plt.plot(time_vector,combined_grid[0,:]*h,'k--',label = 'Hybrid Combined')
    plt.legend()
    plt.show()
    plt.figure()

    # Plot 1: Original Wave Speeds
    plt.subplot(1, 2, 1)
    plt.plot(time_vector, wavespeed_PDE, 'g', label='PDE', linewidth=2)
    plt.plot(time_vector, wavespeed_SSA, 'b', label='SSA', linewidth=2)
    plt.plot(time_vector, wavespeed_hybrid, 'r', label='Hybrid', linewidth=2)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Wave Speed', fontsize=12)
    plt.title('Original Wave Speeds', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ylim(min_y, max_y)  # Set the same y-axis limits

    # Plot 2: Smoothed Wave Speeds
    
    plt.subplot(1, 2, 2)
    plt.plot(time_vector_smooth, wavespeed_PDE_smooth, 'g', label='PDE (Smoothed)', linewidth=2)
    plt.plot(time_vector_smooth, wavespeed_SSA_smooth, 'b', label='SSA (Smoothed)', linewidth=2)
    plt.plot(time_vector_smooth, wavespeed_hybrid_smooth, 'r', label='Hybrid (Smoothed)', linewidth=2)
    plt.plot(time_vector_smooth,max_wave_speed,'k--', label = 'Max theoretical wavespeed')
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Wave Speed', fontsize=12)
    plt.title('Smoothed Wave Speeds', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ylim(min_y, max_y)  # Set the same y-axis limits

    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()

    # print(f"The wavespeed PDE_smooth: {wavespeed_PDE_smooth}")
    # # Plotting
    # fig, axs = plt.subplots(2, 1, figsize=(15, 8), sharex=True)

    # # PDE Usage Heatmap
    # im1 = axs[0].imshow(PDE_usage, aspect='auto', cmap='Blues', extent=[time_vector[0], time_vector[-1], 0, compartments])
    # axs[0].set_title('PDE Usage Over Time')
    # axs[0].set_ylabel('Compartments')
    # plt.colorbar(im1, ax=axs[0], label='PDE Activation (1=Active)')

    # # SSA Usage Heatmap
    # im2 = axs[1].imshow(SSA_usage, aspect='auto', cmap='Reds', extent=[time_vector[0], time_vector[-1], 0, compartments])
    # axs[1].set_title('SSA Usage Over Time')
    # axs[1].set_xlabel('Time')
    # axs[1].set_ylabel('Compartments')
    # plt.colorbar(im2, ax=axs[1], label='SSA Activation (1=Active)')

    # plt.tight_layout()
    # plt.show()

    all_events = [event for run in SSA_events for event in run]
    
    df = pd.DataFrame(all_events, columns=['Time', 'Compartment', 'Event Type'])

    print(df.dtypes)
    print(df.head())
    df['Time'] = pd.to_numeric(df['Time'], errors='coerce')
    df['Compartment'] = pd.to_numeric(df['Compartment'], errors='coerce')

    # Ensure 'Event Type' is a string
    df['Event Type'] = df['Event Type'].astype(str)

    # Ensure Event Type is string
    df['Event Type'] = df['Event Type'].astype(str)
    dictionary_of_events = df['Event Type'].unique()
    print(dictionary_of_events)

    print(f"The PDE update times are {PDE_update_times}")

    def plot_heatmap_of_event(event_type):
        """Plots a heatmap of the specified event type with PDE update density."""
        
        df_filtered = df[df['Event Type'] == event_type]
        # Create    bins for time and compartments
        time_bins = np.linspace(df_filtered['Time'].min(), df_filtered['Time'].max(), 128)
        compartment_bins = np.arange(df_filtered['Compartment'].min(), df_filtered['Compartment'].max() + 1)
    
        # Create a pivot table for event counts
        heatmap_data = pd.crosstab(
            pd.cut(df_filtered['Time'], bins=time_bins),
            pd.cut(df_filtered['Compartment'], bins=compartment_bins)
        )
    
        # Plotting the heatmap
        fig, ax1 = plt.subplots(figsize=(14, 8))
        sns.heatmap(heatmap_data.T, cmap='YlGnBu', cbar_kws={'label': 'Event Count'}, ax=ax1, linewidths=0.1)
    
        # Plot PDE update times on a secondary axis
      
        # Legends and labels
        ax1.set_title(f'Event Density Heatmap for: {event_type}')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Compartment')
 
    
        plt.tight_layout()
        plt.show()

    def plot_PDE_dist(PDE_update_times):
        """Here we plot the distribution of the PDE_distribution"""
    
        # Ensure PDE_update_times is a 1D array or list
        if isinstance(PDE_update_times, np.ndarray):
            PDE_update_times = PDE_update_times.flatten()
        elif not isinstance(PDE_update_times, list):
            raise ValueError("PDE_update_times must be a list or a 1D numpy array")
    
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(PDE_update_times, bins=128, color='green', alpha=0.7)
        ax.set_title('PDE Update Time Distribution')
        ax.set_xlabel('Time')
        ax.set_ylabel('Frequency')
        plt.show()
    # Call the function for 'diffusion'
    
    plot_heatmap_of_event('J degredation')
    # plot_PDE_dist(PDE_update_times)

if __name__ == "__main__":
    main()