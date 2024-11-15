import numpy as np
from production_project import Hybrid, Stochastic, PDE
import subprocess


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

except ValueError as e:
    print("Error:", e)



# domain_length = 1 #Length of the domain
# compartment_number = 8 #Number of compartments
# PDE_multiple = 8 #How many PDE's points per cell (ie we make it a multiple times more fine)
# total_time = 200 #The total time to run the system
# timestep = 0.05 #The time step
# particles_per_compartment_thresh = 10 #The max number of particles per compartment in which the conversion begins
# gamma = 0.2 #The rate of conversion
# production_rate = 2 #The rate of production across the entire sim (this is later changed to be per cell, multiplied by h)
# degredation_rate = 0.01 #The rate of degredation
# diffusion_rate = (domain_length**2)*(10e-3) #The rate of diffusion (Scale down by L^2) Look at courant number
# number_particles_per_cell = 1 #Number of particles initially per compartment

repeats = 100

SSA_initial= np.ones((compartment_number), np.int64) * number_particles_per_cell #Initial conditions (within each cell) 

Model = Hybrid(domain_length, compartment_number, PDE_multiple,total_time, timestep, particles_per_compartment_thresh, gamma, production_rate, degradation_rate, diffusion_rate, SSA_initial)

D_grid,C_grid,combined_grid = Model.run_simulation(number_of_repeats=repeats)
Model.save_simulation_data(D_grid,C_grid,combined_grid, datadirectory='data')

SSA_model =  Stochastic(domain_length, compartment_number, total_time, timestep, production_rate, degradation_rate, diffusion_rate, SSA_initial) 
SSA_grid = SSA_model.run_simulation(number_of_repeats=repeats)
SSA_model.save_simulation_data(SSA_grid, datadirectory='data') #ignore

PDE_points = Model.PDE_M
PDE_initial = np.ones_like(PDE_points) * number_particles_per_cell/Model.h
PDE_Model = PDE(domain_length, PDE_points, total_time, timestep, production_rate, degradation_rate, diffusion_rate, PDE_initial)
PDE_grid = PDE_Model.run_simulation()
PDE_Model.save_simulation_data(PDE_grid, datadirectory='data')
subprocess.run(['python','animate.py'])