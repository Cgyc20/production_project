import numpy as np
from production_project import Hybrid, Stochastic, PDE

def main():
    """
    Self-contained demonstration of Hybrid, Stochastic, and PDE simulations
    using hardcoded parameters.
    """

    # --------------------------
    # Simulation Parameters
    # --------------------------
    domain_length = 5.0
    h = 0.1
    compartment_number = int(domain_length / h)

    PDE_multiple = 8
    total_time = 10
    timestep = 0.008
    particles_per_compartment_thresh = 50
    gamma = 1.0
    production_rate = 10.0
    degradation_rate = 0.01
    diffusion_rate = 1e-2
    repeats = 10
    number_particles_per_cell = 1 #This will be plugged into the initial conditions

    # --------------------------
    # Initial Conditions
    # --------------------------
    SSA_initial = np.zeros(compartment_number, dtype=np.int64)
    SSA_initial[0] = number_particles_per_cell  # start with one particle in the first compartment

    # --------------------------
    # Hybrid Model Simulation
    # --------------------------
    hybrid_model = Hybrid(
        domain_length, compartment_number, PDE_multiple, total_time, timestep,
        particles_per_compartment_thresh, gamma, production_rate, degradation_rate,
        diffusion_rate, SSA_initial, use_c_functions=True
    )

    Hybrid_SSA, Hybrid_PDE, Hybrid_combined = hybrid_model.run_simulation(number_of_repeats=repeats)
    hybrid_model.save_simulation_data(Hybrid_SSA, Hybrid_PDE, Hybrid_combined, datadirectory='data', filename='Hybrid_data')

    # --------------------------
    # Pure Stochastic Simulation
    # --------------------------
    stochastic_model = Stochastic(
        domain_length, compartment_number, total_time, timestep,
        production_rate, degradation_rate, diffusion_rate, SSA_initial
    )

    SSA_grid = stochastic_model.run_simulation(number_of_repeats=repeats)
    stochastic_model.save_simulation_data(SSA_grid, datadirectory='data')

    # --------------------------
    # PDE Simulation
    # --------------------------
    PDE_points = hybrid_model.PDE_M
    PDE_initial = np.zeros(PDE_points)
    PDE_initial[0:PDE_multiple] = number_particles_per_cell / hybrid_model.h

    pde_model = PDE(
        domain_length, PDE_points, total_time, timestep,
        production_rate, degradation_rate, diffusion_rate, PDE_initial
    )

    PDE_grid = pde_model.run_simulation()
    pde_model.save_simulation_data(PDE_grid, datadirectory='data')

    # --------------------------
    # Quick Summary Printout
    # --------------------------
    print("Simulation complete!")
    print(f"Steady-state estimate (production/degradation): {production_rate / degradation_rate}")
    print(f"Initial SSA: {SSA_initial}")
    print(f"Hybrid SSA final mean: {Hybrid_SSA[:, -1].mean():.2f}")
    print(f"PDE final mean: {PDE_grid[:, -1].mean():.2f}")

if __name__ == "__main__":
    main()
