import numpy as np
from tqdm import tqdm
import os
import json

class Stochastic:
    def __init__(self, domain_length, compartment_number, total_time, timestep, production_rate, degredation_rate, diffusion_rate, SSA_initial):
        self.L = domain_length
        self.SSA_M = compartment_number
        self.production_rate = production_rate
        self.total_time = total_time
        self.timestep = timestep
        
        self.degredation_rate = degredation_rate
        
        self.h = self.L / compartment_number
        self.diffusion_rate = diffusion_rate*self.h
        self.production_rate_per_compartment = production_rate*self.h
        self.d = diffusion_rate / (self.h ** 2)  # Jump rate in SSA

        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M)

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA initial is not a np array")
        elif not len(SSA_initial) == compartment_number:
            raise ValueError("The length of the SSA initial is not the same as compartment number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("The SSA initial is not an integer")
        else:
            self.SSA_initial = SSA_initial

        self.steady_state = production_rate / degredation_rate
        self.time_vector = np.arange(0, total_time, timestep)  # The time vector
        print("Successfully initialized the Stochastic model")

        # print(f"initial = {SSA_initial}")
    
    def create_initial_dataframe(self):
        """Creates the intiial dataframes to be used throughout the simulation
        Returns: the initial dataframes for discrete and continuous numbers of molecules. With initial conditions"""

        SSA_matrix = np.zeros((self.SSA_M,len(self.time_vector) ))  # Discrete molecules grid
        SSA_matrix[:, 0] = self.SSA_initial
        # print(f"Shape of SSA_grid = {SSA_grid.shape}")  # Debugging statement
        # print(f"Initial SSA_grid: {SSA_grid}")  # Debugging statement

        return SSA_matrix

    def propensity_calculation(self,SSA_list):
        """Calculate the propensity functions for each reaction"""
        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]  # Left boundary condition (only move right)
        movement_propensity[-1] = self.d * SSA_list[-1]  # Right boundary condition (only move left)
        movement_propensity = np.maximum(movement_propensity, 0)

        R1_propensity = self.production_rate_per_compartment * np.ones_like(SSA_list)  
        R2_propensity = self.degredation_rate * SSA_list

        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity))

        return combined_propensity
    

    def stochastic_simulation(self,SSA_grid):
        t = 0
        old_time = t
        SSA_list = SSA_grid[:, 0]  # Starting SSA_list
        # print(f"SSA list = {SSA_list}")
        while t < self.total_time:
            total_propensity = self.propensity_calculation(SSA_list)
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
                    SSA_list[index] = max(SSA_list[index] - 1, 0)
                    SSA_list[index - 1] += 1
                else:  # Move right
                    SSA_list[index] = max(SSA_list[index] - 1, 0)
                    SSA_list[index + 1] += 1
            elif index == 0:  # Left boundary (can only move right)
                SSA_list[index] = max(SSA_list[index] - 1, 0)
                SSA_list[index + 1] += 1
            elif index == self.SSA_M - 1:  # Right boundary (can only move left)
                SSA_list[index] = max(SSA_list[index] - 1, 0)
                SSA_list[index - 1] += 1
            elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:  # Production reaction
                SSA_list[compartment_index] += 1
            elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:  # Degradation reaction
                SSA_list[compartment_index] = max(SSA_list[compartment_index] - 1, 0)

                
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
            # print(f"Initial SSA_grid: {SSA_grid_initial}")  # Debugging statement
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
            'degredation_rate': self.degredation_rate,
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




