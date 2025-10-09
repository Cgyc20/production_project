import numpy as np
from production_project import Hybrid

def main():
    # ---- PARAMETERS ---- #
    domain_length = 5.0
    h = 0.1
    PDE_multiple = 8
    total_time = 8
    timestep = 0.005
    particles_per_compartment_thresh = 50
    production_rate = 10
    degradation_rate = 0.01
    diffusion_rate = 1e-2
    number_particles_per_cell = 100

    repeats_per_batch = 100
    total_repeats = 500

    # ---- DERIVED VALUES ---- #
    compartment_number = int(domain_length / h)
    print(f"Compartment number: {compartment_number}")
    print(f"Steady state: {production_rate / degradation_rate}")

    # ---- INITIAL CONDITIONS ---- #
    SSA_initial = np.zeros(compartment_number, dtype=np.int64)
    SSA_initial[0:5] = number_particles_per_cell

    # ---- GAMMA LOOP ---- #
    for gamma in range(1, 11):
        print(f"\nRunning simulations for gamma = {gamma}")

        # Initialize Hybrid model once per gamma
        hybrid_model = Hybrid(
            domain_length, compartment_number, PDE_multiple,
            total_time, timestep, particles_per_compartment_thresh,
            gamma, production_rate, degradation_rate,
            diffusion_rate, SSA_initial, use_c_functions=True
        )

        # ---- BATCH LOOP ---- #
        # ---- BATCH LOOP ---- #
        for batch_start in range(0, total_repeats, repeats_per_batch):
            print(f"  Running batch {batch_start // repeats_per_batch + 1}")
            Hybrid_SSA, Hybrid_PDE, Hybrid_combined= hybrid_model.run_simulation(
                number_of_repeats=repeats_per_batch
            )

            # Pass only base filename, let function handle directory and extension
            filename = f'Hybrid_data_gamma_{gamma}_batch_{batch_start // repeats_per_batch + 1}'
            hybrid_model.save_simulation_data(Hybrid_SSA, Hybrid_PDE, Hybrid_combined,
                                            datadirectory='data', filename=filename)
            print(f"  Saved batch to {filename}")


if __name__ == "__main__":
    main()