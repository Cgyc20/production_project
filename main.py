import numpy as np
from production_project import Hybrid, Stochastic, PDE
import subprocess
import sys

def main():
    # Open file and read parameters into a dictionary
    with open("parameter_input.dat", "r") as f:
        parameters_dict = {line.split()[0]: line.split()[1] for line in f}

    # Convert values to the appropriate types, with error handling
    try:
        domain_length = float(parameters_dict.get('domain_length', 0))
        compartment_number = int(parameters_dict.get('compartment_number', 0))
        if compartment_number == 0:
            raise ValueError("compartment_number cannot be zero.")
        
        PDE_multiple = int(parameters_dict.get('PDE_multiple', 0))
        total_time = int(parameters_dict.get('total_time', 0))
        timestep = float(parameters_dict.get('timestep', 0))
        particles_per_compartment_thresh = int(parameters_dict.get('particles_per_compartment_thresh', 0))
        gamma = float(parameters_dict.get('gamma', 0))
        production_rate = float(parameters_dict.get('production_rate', 0))
        degradation_rate = float(parameters_dict.get('degradation_rate', 0))
        repeats = int(parameters_dict.get('repeats', 0))
        # Ensure number_particles_per_cell is defined and not zero
        number_particles_per_cell = int(parameters_dict.get('number_particles_per_cell', 0))
        if number_particles_per_cell == 0:
            raise ValueError("number_particles_per_cell cannot be zero.")
        
        # Calculate diffusion_rate using domain_length
        diffusion_rate = (domain_length ** 2) * (10e-3)

        # Print the values to confirm
        print(f"domain_length: {domain_length}")
        print(f"compartment_number: {compartment_number}")
        print(f"PDE_multiple: {PDE_multiple}")
        print(f"total_time: {total_time}")
        print(f"timestep: {timestep}")
        print(f"particles_per_compartment_thresh: {particles_per_compartment_thresh}")
        print(f"gamma: {gamma}")
        print(f"production_rate: {production_rate}")
        print(f"degradation_rate: {degradation_rate}")
        print(f"diffusion_rate: {diffusion_rate}")
        print(f"number_particles_per_cell: {number_particles_per_cell}")
        print(f"The steady state is: {production_rate/degradation_rate}")

    except ValueError as e:
        print("Error:", e)
        return

    """Initialise the hybrid model"""
    # np.random.seed(0)
    SSA_initial = np.zeros((compartment_number), np.int64) # Initial conditions (within each cell)

    
    SSA_initial[compartment_number//2-1:compartment_number//2+1] = number_particles_per_cell


    # multiply_vector = np.arange(0, compartment_number)%2
    
    # SSA_initial = SSA_initial* multiply_vector


    Model = Hybrid(domain_length, compartment_number, PDE_multiple, total_time, timestep, particles_per_compartment_thresh, gamma, production_rate, degradation_rate, diffusion_rate, SSA_initial, use_c_functions=False) # Define the hybrid model

    Hybrid_SSA, Hybrid_PDE, Hybrid_combined = Model.run_simulation(number_of_repeats=repeats)
    Model.save_simulation_data(Hybrid_SSA, Hybrid_PDE, Hybrid_combined, datadirectory='data')

    SSA_model = Stochastic(domain_length, compartment_number, total_time, timestep, production_rate, degradation_rate, diffusion_rate, SSA_initial, use_c_functions=True)
    SSA_grid = SSA_model.run_simulation(number_of_repeats=repeats)
    SSA_model.save_simulation_data(SSA_grid, datadirectory='data') # ignore


    PDE_points = Model.PDE_M
    PDE_initial = np.ones_like(PDE_points) * number_particles_per_cell / Model.h
    print(PDE_initial)
    PDE_Model = PDE(domain_length, PDE_points, total_time, timestep, production_rate, degradation_rate, diffusion_rate, PDE_initial)
    PDE_grid = PDE_Model.run_simulation()
    PDE_Model.save_simulation_data(PDE_grid, datadirectory='data')

    print(f"PDE grid at timestep one: {PDE_grid[:,0]}")
    print(f"Stochastic grid at timestep 1: {SSA_grid[:,0]}")
    print(f"SSA in hybrid model at timestep 1: {Hybrid_SSA[:,0]}")

def plot():
    subprocess.run(['python', 'animate.py'])

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        if sys.argv[1] == 'plot':
            plot()
        elif sys.argv[1] == 'run':
            main()
        elif sys.argv[1] == 'run_and_plot':
            main()
            plot()
        else:
            print("Invalid argument. Use 'run', 'plot', or 'run_and_plot'.")
    else:
        print("No argument provided. Use 'run', 'plot', or 'run_and_plot'.")