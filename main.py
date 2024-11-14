import numpy as np
from production_project import Hybrid, Stochastic, PDE
import subprocess


domain_length = 1 #Length of the domain
compartment_number = 8 #Number of compartments
PDE_multiple = 16 #How many PDE's points per cell (ie we make it a multiple times more fine)
total_time = 200 #The total time to run the system
timestep = 0.05 #The time step

repeats = 100

# Make it total number per compartment

particles_per_compartment_thresh = 
threshold_conc = 
threshold_conc = 100 #The threshold for the SSA to switch to the continuous regime
particles_per_compartment_thresh = threshold_conc*(domain_length/compartment_number)
print(f"The threshold particles per compartment will be: {particles_per_compartment_thresh}")
gamma = 0.4 #The rate of conversion
production_rate = 2 #The rate of production across the entire sim (this is later changed to be per cell, multiplied by h)
degredation_rate = 0.01 #The rate of degredation
diffusion_rate = (domain_length**2)*(10e-3) #The rate of diffusion (Scale down by L^2) Look at courant number
number_particles_per_cell = 2 #Number of particles initially per compartment
SSA_initial= np.ones((compartment_number), np.int64) * number_particles_per_cell #Initial conditions (within each cell) 

Model = Hybrid(domain_length, compartment_number, PDE_multiple,total_time, timestep, threshold_conc, gamma, production_rate, degredation_rate, diffusion_rate, SSA_initial)

D_grid,C_grid,combined_grid = Model.run_simulation(number_of_repeats=repeats)
Model.save_simulation_data(D_grid,C_grid,combined_grid, datadirectory='data')

SSA_model =  Stochastic(domain_length, compartment_number, total_time, timestep, production_rate, degredation_rate, diffusion_rate, SSA_initial) 
SSA_grid = SSA_model.run_simulation(number_of_repeats=repeats)
SSA_model.save_simulation_data(SSA_grid, datadirectory='data') #ignore

PDE_points = Model.PDE_M
PDE_initial = np.ones_like(PDE_points) * number_particles_per_cell/Model.h
PDE_Model = PDE(domain_length, PDE_points, total_time, timestep, production_rate, degredation_rate, diffusion_rate, PDE_initial)
PDE_grid = PDE_Model.run_simulation()
PDE_Model.save_simulation_data(PDE_grid, datadirectory='data')
subprocess.run(['python','animate.py'])