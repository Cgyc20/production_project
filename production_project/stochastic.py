import numpy as np
from tqdm import tqdm
import os
import json
from copy import deepcopy
import ctypes
from production_project.clibrary_argtypes import set_clibrary_argtypes  # Each data type for the C functions

clibrary = ctypes.CDLL("c_class/clibrary.so")  # Import the C library
set_clibrary_argtypes(clibrary)  # Set data types for C functions


class Stochastic:
    """
    Implements the stochastic Fisher-KPP system using Gillespie SSA.

    Attributes:
        L (float): Total domain length.
        SSA_M (int): Number of compartments.
        production_rate (float): Rate of production reactions.
        degradation_rate (float): Rate of degradation reactions.
        total_time (float): Total simulation time.
        timestep (float): Time increment for recording states.
        h (float): Spatial step size (domain_length / compartment_number).
        diffusion_rate (float): PDE diffusion coefficient.
        d (float): Jump rate in SSA (diffusion_rate / h^2).
        SSA_X (np.ndarray): Spatial positions of compartments.
        SSA_initial (np.ndarray): Initial discrete particle counts.
        time_vector (np.ndarray): Vector of time points for simulation.
    """

    def __init__(self, domain_length, compartment_number, total_time, timestep,
                 production_rate, degradation_rate, diffusion_rate, SSA_initial):
        """
        Initialize the stochastic Fisher-KPP model.

        Args:
            domain_length (float): Total spatial length of the domain.
            compartment_number (int): Number of discrete SSA compartments.
            total_time (float): Total simulation time.
            timestep (float): Time increment for recording SSA states.
            production_rate (float): Rate of production reactions.
            degradation_rate (float): Rate of degradation reactions.
            diffusion_rate (float): Diffusion coefficient for SSA jumps.
            SSA_initial (np.ndarray): Initial discrete particle counts (integer array).
        """
        self.L = domain_length
        self.SSA_M = compartment_number
        self.production_rate = production_rate
        self.degradation_rate = degradation_rate
        self.total_time = total_time
        self.timestep = timestep

        self.h = self.L / compartment_number
        print(f"h in the stochastic method is: {self.h}")
        self.diffusion_rate = diffusion_rate
        self.d = diffusion_rate / (self.h ** 2)  # Jump rate in SSA
        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M)

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA_initial must be a numpy array")
        elif len(SSA_initial) != compartment_number:
            raise ValueError("SSA_initial length must match compartment_number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("SSA_initial must contain integers")
        else:
            self.SSA_initial = SSA_initial

        self.time_vector = np.arange(0, total_time, timestep)
        print("Successfully initialized the Stochastic model")

    def create_initial_dataframe(self) -> np.ndarray:
        """
        Creates the initial SSA grid for the simulation.

        Returns:
            np.ndarray: A 2D array of shape (compartments x time_points) initialized with SSA_initial.
        """
        SSA_matrix = np.zeros((self.SSA_M, len(self.time_vector)))
        SSA_matrix[:, 0] = self.SSA_initial
        return SSA_matrix

    def propensity_calculationPython(self, SSA_list: np.ndarray) -> np.ndarray:
        """
        Calculate the propensity of each possible reaction in Python.

        Args:
            SSA_list (np.ndarray): Current SSA particle counts.

        Returns:
            np.ndarray: Combined propensities for movement, production, and degradation.
        """
        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]  # Left boundary
        movement_propensity[-1] = self.d * SSA_list[-1]  # Right boundary

        production_propensity = self.production_rate * SSA_list
        degradation_propensity = self.degradation_rate / self.h * SSA_list * (SSA_list - 1)

        return np.hstack((movement_propensity, production_propensity, degradation_propensity))

    def propensity_calculation(self, SSA_list: np.ndarray) -> np.ndarray:
        """
        Wrapper for propensity calculation (can switch to C implementation if desired).

        Args:
            SSA_list (np.ndarray): Current SSA particle counts.

        Returns:
            np.ndarray: Combined propensities for all reactions.
        """
        return self.propensity_calculationPython(SSA_list)

    def stochastic_simulation(self, SSA_grid: np.ndarray) -> np.ndarray:
        """
        Run the stochastic Gillespie simulation.

        Args:
            SSA_grid (np.ndarray): Initial SSA grid.

        Returns:
            np.ndarray: SSA grid updated over the simulation time.
        """
        t = 0
        old_time = t
        SSA_list = deepcopy(SSA_grid[:, 0])

        while t < self.total_time:
            total_propensity = self.propensity_calculationPython(SSA_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:
                break  # Stop if no reactions can occur

            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)
            alpha_cum = np.cumsum(total_propensity)
            index = np.searchsorted(alpha_cum, r2 * alpha0)
            compartment_index = index % self.SSA_M

            # Diffusion events
            if 1 <= index <= self.SSA_M - 2:
                if r3 < 0.5:
                    SSA_list[index] -= 1
                    SSA_list[index - 1] += 1
                else:
                    SSA_list[index] -= 1
                    SSA_list[index + 1] += 1
            elif index == 0:
                SSA_list[index] -= 1
                SSA_list[index + 1] += 1
            elif index == self.SSA_M - 1:
                SSA_list[index] -= 1
                SSA_list[index - 1] += 1
            # Production
            elif self.SSA_M <= index <= 2 * self.SSA_M - 1:
                SSA_list[compartment_index] += 1
            # Degradation
            elif 2 * self.SSA_M <= index <= 3 * self.SSA_M - 1:
                SSA_list[compartment_index] -= 1

            ind_before = np.searchsorted(self.time_vector, old_time, 'right')
            ind_after = np.searchsorted(self.time_vector, t, 'left')
            for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                SSA_grid[:, time_index] = SSA_list

            old_time = t
            t += tau

        return SSA_grid

    def run_simulation(self, number_of_repeats: int) -> np.ndarray:
        """
        Run multiple stochastic simulations and compute the average SSA grid.

        Args:
            number_of_repeats (int): Number of simulation repetitions.

        Returns:
            np.ndarray: Averaged SSA grid over all repeats.
        """
        SSA_average = np.zeros_like(self.create_initial_dataframe())

        for _ in tqdm(range(number_of_repeats), desc="Running the Stochastic simulations"):
            SSA_grid_initial = self.create_initial_dataframe()
            SSA_current = self.stochastic_simulation(SSA_grid_initial)
            SSA_average += SSA_current

        SSA_average /= number_of_repeats
        print("Simulation completed")
        return SSA_average

    def save_simulation_data(self, filled_SSA_grid: np.ndarray, datadirectory='data', filename='Pure_SSA_data'):
        """
        Save the SSA simulation results and parameters to disk.

        Args:
            filled_SSA_grid (np.ndarray): Averaged SSA grid.
            datadirectory (str): Directory to save data files.
            filename (str): Base filename (without extension) to use for saved files.
                            Defaults to 'Pure_SSA_data'.
        """
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

        np.savez(os.path.join(datadirectory, f'{filename}.npz'),
                 SSA_grid=filled_SSA_grid,
                 time_vector=self.time_vector,
                 SSA_X=self.SSA_X
                 )

        with open(os.path.join(datadirectory, f"{filename}_parameters.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)

        print("Data saved successfully")
