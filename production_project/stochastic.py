import numpy as np
from tqdm import tqdm
import os
import json
from copy import deepcopy
import ctypes

clibrary = ctypes.CDLL("c_class/clibrarySSA.so")  # import the c library

clibrary.CalculatePropensity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_float,
    ctypes.c_float,
    ctypes.c_float
]


class Stochastic:
    def __init__(self, domain_length, compartment_number, total_time, timestep, production_rate, degradation_rate, diffusion_rate, SSA_initial, use_c_functions):
        self.L = domain_length
        self.SSA_M = compartment_number
        self.production_rate = production_rate
        self.degradation_rate = degradation_rate
        self.total_time = total_time
        self.timestep = timestep

        self.use_c_functions = use_c_functions
        
        self.h = self.L / compartment_number
        self.diffusion_rate = diffusion_rate
        self.d = diffusion_rate / (self.h ** 2)  # Jump rate in SSA
        self.production_rate_per_compartment = production_rate * self.h
        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M)

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA initial is not a np array")
        elif not len(SSA_initial) == compartment_number:
            raise ValueError("The length of the SSA initial is not the same as compartment number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("The SSA initial is not an integer")
        else:
            self.SSA_initial = SSA_initial

        self.time_vector = np.arange(0, total_time, timestep)  # The time vector
        print("Successfully initialized the Stochastic model")

    def create_initial_dataframe(self):
        """Creates the initial dataframes to be used throughout the simulation
        Returns: the initial dataframes for discrete and continuous numbers of molecules. With initial conditions"""

        SSA_matrix = np.zeros((self.SSA_M, len(self.time_vector)))  # Discrete molecules grid
        SSA_matrix[:, 0] = self.SSA_initial
        return SSA_matrix

    def propensity_calculationPython(self, SSA_list):
        """Calculate the propensity functions for each reaction"""
        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]  # Left boundary condition (only move right)
        movement_propensity[-1] = self.d * SSA_list[-1]  # Right boundary condition (only move left)

        production_propensity = self.production_rate_per_compartment * SSA_list
        degradation_propensity = self.degradation_rate * SSA_list * (SSA_list - 1)

        combined_propensity = np.concatenate((movement_propensity, production_propensity, degradation_propensity))

        return combined_propensity

    def propensity_calculationC(self, SSA_list: np.ndarray) -> np.ndarray:
        """
        Calculates the propensity functions for each reaction using the C library.

        Args:
            SSA_list (np.ndarray): Discrete molecules list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        SSA_list = np.ascontiguousarray(SSA_list, dtype=np.int32)

        # Initialize propensity_list and other arrays
        propensity_list = np.zeros(3 * self.SSA_M, dtype=np.float32)

        # Call the C function to calculate propensities
        clibrary.CalculatePropensity(
            SSA_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            self.SSA_M,
            propensity_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.d,
            self.production_rate,
            self.degradation_rate,
        )

        # Convert propensity list back into numpy floats
        propensity_list = np.ctypeslib.as_array(propensity_list, shape=propensity_list.shape)
        return propensity_list

    def propensity_calculation(self, SSA_list):
        if self.use_c_functions:
            return self.propensity_calculationC(SSA_list)
        else:
            return self.propensity_calculationPython(SSA_list)

    def stochastic_simulation(self, SSA_grid):
        t = 0
        old_time = t
        SSA_list = deepcopy(SSA_grid[:, 0])  # Starting SSA_list
        while t < self.total_time:
            total_propensity = self.propensity_calculationPython(SSA_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:  # Stop if no reactions can occur
                break

            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)  # Time until next reaction
            alpha_cum = np.cumsum(total_propensity)  # Cumulative sum of propensities
            index = np.searchsorted(alpha_cum, r2 * alpha0)  # Determine which reaction occurs

            compartment_index = index % self.SSA_M  # The compartmental index is just the modulo of SSA.
            if index <= self.SSA_M - 2 and index >= 1:
                if r3 < 0.5:  # Move left
                    SSA_list[index] = SSA_list[index] - 1
                    SSA_list[index - 1] += 1
                else:  # Move right
                    SSA_list[index] = SSA_list[index] - 1
                    SSA_list[index + 1] += 1
            elif index == 0:  # Left boundary (can only move right)
                SSA_list[index] = SSA_list[index] - 1
                SSA_list[index + 1] += 1
            elif index == self.SSA_M - 1:  # Right boundary (can only move left)
                SSA_list[index] = SSA_list[index] - 1
                SSA_list[index - 1] += 1
            elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:  # Production reaction
                SSA_list[compartment_index] += 1
            elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:  # Degradation reaction
                SSA_list[compartment_index] = SSA_list[compartment_index] - 1

            ind_before = np.searchsorted(self.time_vector, old_time, 'right')
            ind_after = np.searchsorted(self.time_vector, t, 'left')
            for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                SSA_grid[:, time_index] = SSA_list

            old_time = t  # Update old_time
            t += tau  # Update time by the time step

        return SSA_grid

    def run_simulation(self, number_of_repeats):
        """This will run the simulation with a total number of repeats"""
        SSA_initial = self.create_initial_dataframe()
        SSA_average = np.zeros_like(SSA_initial)
        filled_SSA_grid = np.zeros_like(SSA_initial)

        for _ in tqdm(range(number_of_repeats), desc="Running the Stochastic simulations"):
            SSA_grid_initial = self.create_initial_dataframe()
            SSA_current = self.stochastic_simulation(SSA_grid_initial)
            SSA_average += SSA_current

        filled_SSA_grid = SSA_average / number_of_repeats

        print("Simulation completed")
        return filled_SSA_grid

    def save_simulation_data(self, filled_SSA_grid, datadirectory='data'):
        if not os.path.exists(datadirectory):
            os.makedirs(datadirectory)

        params = {
            'domain_length': self.L,
            'compartment_number': self.SSA_M,
            'total_time': self.total_time,
            'timestep': self.timestep,
            'production_rate': self.production_rate,
            'degradation_rate': self.degradation_rate,
            'diffusion_rate': self.diffusion_rate,
            'initial_SSA': self.SSA_initial.tolist(),
            'h': self.h,
        }

        np.savez(os.path.join(datadirectory, 'Pure_SSA_data'),
                 SSA_grid=filled_SSA_grid,
                 time_vector=self.time_vector,
                 SSA_X=self.SSA_X
                 )

        with open(os.path.join(datadirectory, "parameters_pure_SSA.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)

        print("Data saved successfully")