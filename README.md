
# Hybrid Stochastic-PDE Simulation

## Overview

This project provides a **self-contained Python framework** for simulating reaction-diffusion systems using **Hybrid, Stochastic, and PDE models**. It is designed for researchers and practitioners who want to explore hybrid stochastic-PDE dynamics in spatially distributed systems.

The framework implements the **Spatial Regime Conversion Method (SRCM)** for the Fisher-KPP adaptive switching problem (Test Case 4), as described in:

*Cameron, Yates, and Smith, "Regime Conversion Method."*

- **Hybrid Model:** Combines discrete stochastic simulations (SSA) with continuous PDEs for efficient and accurate modeling.
- **Stochastic Model:** Pure SSA simulation of particles in compartments.
- **PDE Model:** Pure partial differential equation simulation of the system.

The code is fully self-contained with hardcoded parameters for demonstration purposes.

## Requirements

- Python 3.8+
- Base packages:
  - `numpy`
  - `tqdm`
  - `os`
  - `json`
  - `ctypes`
- `production_project` package (contains `Hybrid`, `Stochastic`, `PDE` classes)
- Optional: C compiler for hybrid model acceleration

> **Note:** The hybrid model uses `ctypes` to call a compiled C library (`.so`) for speed. On macOS, this `.so` is provided. For Windows, you may need to recompile the library or simply disable C acceleration by setting `use_c_functions=False`.

---


## Usage

Run the main simulation script:

```bash
python main.py
```

### How it works

`main.py` is a self-contained demonstration. Here’s what it does:

1. **Set Simulation Parameters**

```python
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
number_particles_per_cell = 1
```

2. **Set Initial Conditions**

```python
SSA_initial = np.zeros(compartment_number, dtype=np.int64)
SSA_initial[0] = number_particles_per_cell
```

3. **Run Hybrid Simulation**

```python
hybrid_model = Hybrid(
    domain_length, compartment_number, PDE_multiple, total_time, timestep,
    particles_per_compartment_thresh, gamma, production_rate, degradation_rate,
    diffusion_rate, SSA_initial, use_c_functions=True
)

Hybrid_SSA, Hybrid_PDE, Hybrid_combined = hybrid_model.run_simulation(number_of_repeats=repeats)
hybrid_model.save_simulation_data(Hybrid_SSA, Hybrid_PDE, Hybrid_combined, datadirectory='data', filename='Hybrid_data')
```

4. **Run Pure Stochastic Simulation**

```python
stochastic_model = Stochastic(
    domain_length, compartment_number, total_time, timestep,
    production_rate, degradation_rate, diffusion_rate, SSA_initial
)

SSA_grid = stochastic_model.run_simulation(number_of_repeats=repeats)
stochastic_model.save_simulation_data(SSA_grid, datadirectory='data', filename = 'pure_SSA_data')
```

5. **Run Pure PDE Simulation**

```python
PDE_points = hybrid_model.PDE_M
PDE_initial = np.zeros(PDE_points)
PDE_initial[0:PDE_multiple] = number_particles_per_cell / hybrid_model.h

pde_model = PDE(
    domain_length, PDE_points, total_time, timestep,
    production_rate, degradation_rate, diffusion_rate, PDE_initial
)

PDE_grid = pde_model.run_simulation()
pde_model.save_simulation_data(PDE_grid, datadirectory='data', filename = 'PDE_data')
```

6. **Completion Message**

```python
print("Simulation complete!")
```

---

## Visualisation

After running `main.py`, you can animate the results:

```bash
python animate.py
```

This will create visualisations of the hybrid simulation outputs. Make sure `main.py` has been run first so that the `data/` folder contains the simulation outputs. You can change the datapath at the top of the `animate.py` script.

---

## Output

* `data/Hybrid_data.npz`: Hybrid simulation results including SSA, PDE, and combined grids.
* `data/Hybrid_data_parameters.json`: JSON file containing all simulation parameters.
* Additional outputs from stochastic and PDE simulations in the `data/` folder.

---