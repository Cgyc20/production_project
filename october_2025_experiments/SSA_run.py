import numpy as np
from production_project import Stochastic, PDE

def main():
    # ---- PARAMETERS ---- #
    domain_length = 5.0
    h = 0.1
    PDE_multiple = 8
    total_time = 8
    timestep = 0.005
    production_rate = 10
    degradation_rate = 0.01
    diffusion_rate = 1e-2
    repeats = 500
    number_particles_per_cell = 100

    # ---- DERIVED VALUES ---- #
    compartment_number = int(domain_length / h)
    print(f"Compartment number: {compartment_number}")
    print(f"Steady state: {production_rate / degradation_rate}")

    # ---- INITIAL CONDITIONS ---- #
    SSA_initial = np.zeros(compartment_number, dtype=np.int64)
    SSA_initial[0:5] = number_particles_per_cell

    # ---- STOCHASTIC MODEL ---- #
    SSA_model = Stochastic(domain_length, compartment_number, total_time, timestep,
                           production_rate, degradation_rate, diffusion_rate, SSA_initial)
    SSA_grid = SSA_model.run_simulation(number_of_repeats=repeats)
    SSA_model.save_simulation_data(SSA_grid, datadirectory='data')
    print("SSA simulation complete.")

    # ---- PDE MODEL ---- #
    PDE_points = compartment_number * PDE_multiple  # approximate fine grid for PDE
    PDE_initial = np.zeros(PDE_points)
    PDE_initial[0:5*PDE_multiple] = number_particles_per_cell / h

    PDE_model = PDE(domain_length, PDE_points, total_time, timestep,
                    production_rate, degradation_rate, diffusion_rate, PDE_initial)
    PDE_grid = PDE_model.run_simulation()
    PDE_model.save_simulation_data(PDE_grid, datadirectory='data')
    print("PDE simulation complete.")

    # ---- PRINT INITIAL STATES ---- #
    print(f"PDE grid at t=0: {PDE_grid[:, 0]}")
    print(f"SSA at t=0: {SSA_grid[:, 0]}")

if __name__ == "__main__":
    main()