import numpy as np
import time
import csv
import os
from production_project import Hybrid

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
    number_particles_per_cell = 100

    repeats = 50
    base_threshold = 50
    base_gamma = 1

    # ---- DERIVED VALUES ---- #
    compartment_number = int(domain_length / h)
    SSA_initial = np.zeros(compartment_number, dtype=np.int64)
    SSA_initial[0:5] = number_particles_per_cell

    print(f"Compartment number: {compartment_number}")
    print(f"Steady state: {production_rate / degradation_rate}\n")

    # ---- OUTPUT DIRECTORY ---- #
    save_dir = "/Users/charliecameron/CodingHub/PhD/RCM/1D/production_project/october_2025_experiments"
    os.makedirs(save_dir, exist_ok=True)

    # =====================================================
    # GAMMA SWEEP
    # =====================================================
    print("=== Gamma Sweep ===")
    gamma_results = [("gamma", "avg_time_per_sim")]

    for gamma in range(1, 11):  # gamma = 1..10
        hybrid_model = Hybrid(
            domain_length, compartment_number, PDE_multiple,
            total_time, timestep, base_threshold,
            gamma, production_rate, degradation_rate,
            diffusion_rate, SSA_initial, use_c_functions=True
        )

        start_time = time.time()
        hybrid_model.run_simulation(number_of_repeats=repeats)
        elapsed = time.time() - start_time
        avg_time = elapsed / repeats

        gamma_results.append((gamma, avg_time))
        print(f"Gamma = {gamma:<3} | Avg time per sim: {avg_time:.4f} s")

    gamma_path = os.path.join(save_dir, "gamma.csv")
    with open(gamma_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(gamma_results)
    print(f"\nGamma results saved to {gamma_path}")

    # =====================================================
    # THRESHOLD SWEEP
    # =====================================================
    print("\n=== Threshold Sweep ===")
    threshold_results = [("threshold", "avg_time_per_sim")]

    for threshold in range(10, 81, 10):  # threshold = 10,20,...,80
        hybrid_model = Hybrid(
            domain_length, compartment_number, PDE_multiple,
            total_time, timestep, threshold,
            base_gamma, production_rate, degradation_rate,
            diffusion_rate, SSA_initial, use_c_functions=True
        )

        start_time = time.time()
        hybrid_model.run_simulation(number_of_repeats=repeats)
        elapsed = time.time() - start_time
        avg_time = elapsed / repeats

        threshold_results.append((threshold, avg_time))
        print(f"Threshold = {threshold:<3} | Avg time per sim: {avg_time:.4f} s")

    threshold_path = os.path.join(save_dir, "threshold.csv")
    with open(threshold_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(threshold_results)
    print(f"\nThreshold results saved to {threshold_path}")


if __name__ == "__main__":
    main()
