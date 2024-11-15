import numpy as np

"""Created 1/10/2024 as new project file"""
from tqdm import tqdm
import os
import json
from copy import deepcopy

class Hybrid:
    
    def __init__(self, domain_length, compartment_number, PDE_multiple, total_time, timestep, threshold, gamma, production_rate, degredation_rate, diffusion_rate, SSA_initial):
        self.L = domain_length
        self.SSA_M = compartment_number
        self.PDE_multiple = PDE_multiple
        self.production_rate = production_rate
        self.PDE_M = compartment_number * PDE_multiple+1
        self.deltax = self.L / self.PDE_M
        self.total_time = total_time
        self.timestep = timestep
        self.threshold = threshold


        self.gamma = gamma
        
        self.degredation_rate = degredation_rate
        
        self.h = self.L / compartment_number
        # self.diffusion_rate = diffusion_rate*self.h
        self.diffusion_rate = diffusion_rate
        self.production_rate_per_compartment = production_rate*self.h
        self.d = diffusion_rate / (self.h ** 2)  # Jump rate in SSA

        self.threshold_conc = threshold/self.h
        
    


        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M)
        self.PDE_X = np.linspace(0, self.L, self.PDE_M)

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA initial is not a np array")
        elif not len(SSA_initial) == compartment_number:
            raise ValueError("The length of the SSA initial is not the same as compartment number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("The SSA initial is not an integer")
        else:
            self.SSA_initial = SSA_initial

        self.PDE_initial_conditions = np.zeros_like(self.PDE_X, dtype=np.float64)
        self.steady_state = production_rate / degredation_rate
        self.DX_NEW = self.create_finite_difference()  # Ensure DX_NEW is initialized here
        self.time_vector = np.arange(0, total_time, timestep)  # The time vector
        self.Crank_matrix, self.M1_inverse = self.create_crank_nicholson() #The crank method

        print("Successfully initialized the hybrid model")

        print(f"The threshold concentration is: {self.threshold_conc}")
    def create_crank_nicholson(self):
        """Creates the matrix used for crank nicholson method"""

        H = self.create_finite_difference()

        M1 = np.identity(H.shape[0]) *(1+0.5*self.timestep*self.degredation_rate) - 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H
        M2 = np.identity(H.shape[0]) *(1-0.5*self.timestep*self.degredation_rate) + 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H

        M1_inverse = np.linalg.inv(M1)
        Crank_matrix = M1_inverse@M2
        
        return Crank_matrix, M1_inverse


    def create_finite_difference(self):
        """Creates finite difference matrix"""

        self.DX = np.zeros((self.PDE_M, self.PDE_M), dtype=int)  # Initialize the finite difference matrix

        self.DX[0, 0], self.DX[-1, -1] = -1, -1  # This is the zero-flux boundary conditions
        self.DX[0, 1], self.DX[-1, -2] = 1, 1
        # Fill over all the elements with the finite difference operators
        for i in range(1, self.DX.shape[0] - 1):
            self.DX[i, i] = -2  # Fill diagonal with the -2 terms (spatial operator)
            self.DX[i, (i + 1)] = 1
            self.DX[i, (i - 1)] = 1
        
        return self.DX
    
    def create_initial_dataframe(self):
        """Creates the intiial dataframes to be used throughout the simulation
        Returns: the initial dataframes for discrete and continuous numbers of molecules. With initial conditions"""

        D_grid = np.zeros((self.SSA_M,len(self.time_vector) ))  # Discrete molecules grid
        D_grid[:, 0] = self.SSA_initial
        C_grid = np.zeros((self.PDE_M, len(self.time_vector)))  # Continuous mass
        C_grid[:, 0] = self.PDE_initial_conditions

        return C_grid,D_grid 
        
 
   
    def crank_nicholson(self,old_vector):
        """Returns the new vector using the crank nicholson matrix"""

        return self.Crank_matrix@old_vector #The new vector!

    def approximate_mass_left_hand(self, C_list):
        """Works out the approximate mass of the PDE domain over each grid point. Using the left hand rule"""
        approximation_number_cont = np.zeros(self.SSA_M)

        for i in range(self.SSA_M):
            start_index = self.PDE_multiple * i
            end_index = self.PDE_multiple * (i + 1)
            sum_value = 0.0
            sum_value = np.sum(C_list[start_index:end_index])*self.deltax
            approximation_number_cont[i] = sum_value

        return approximation_number_cont


    def propensity_calculation(self, D_list, C_list):
        """
        Calculates the propensity functions for each reaction.

        Args:
            D_list (np.ndarray): Discrete molecules list.
            C_list (np.ndarray): Continuous mass list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        movement_propensity = 2 * self.d * D_list
        movement_propensity[0] = self.d * D_list[0]
        movement_propensity[-1] = self.d * D_list[-1]
        # movement_propensity = np.maximum(movement_propensity, 0)

        R1_propensity = self.production_rate_per_compartment * np.ones_like(D_list)
        R2_propensity = self.degredation_rate * D_list

        approximate_mass_list = self.approximate_mass_left_hand(C_list)
        combined_list = np.add(D_list, approximate_mass_list)

        conversion_to_discrete = np.zeros_like(D_list)
        conversion_to_cont = np.zeros_like(approximate_mass_list)

        # Ensure the boolean index matches the array size
        if combined_list.shape != D_list.shape:
            raise ValueError("Shape mismatch between combined_list and D_list")
    
        conversion_to_cont[combined_list >= self.threshold] = D_list[combined_list >= self.threshold] * self.gamma
        conversion_to_discrete[combined_list < self.threshold] = approximate_mass_list[combined_list < self.threshold] * self.gamma

        # Debug: Print the propensity values for conversion
        # print(f"Conversion to Continuous Propensity: {conversion_to_cont}")
        # print(f"Conversion to Discrete Propensity: {conversion_to_discrete}")
        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity

    def hybrid_simulation(self,D_grid, C_grid):
        t = 0
        old_time = t
        td = self.timestep
        D_list = D_grid[:, 0]  # Starting D_list
        C_list = C_grid[:, 0]  # Starting C_list

        while t < self.total_time:
            total_propensity = self.propensity_calculation(D_list, C_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:  # Stop if no reactions can occur
                break

            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)  # Time until next reaction
            alpha_cum = np.cumsum(total_propensity)  # Cumulative sum of propensities
            # print(f"Length of alpha_cum: {len(alpha_cum)}")
            index = np.searchsorted(alpha_cum, r2 * alpha0)  # Determine which reaction occurs

            compartment_index = index%self.SSA_M #The compartmental index is just the modulo of SSA. 
            if t + tau <= td:  # Execute Gillespie
                if index <= self.SSA_M - 2 and index >= 1:
                    # print("Reaction diffusion")
                    if r3 < 0.5:  # Move left
                        D_list[index] = max(D_list[index] - 1, 0)
                        D_list[index - 1] += 1
                    else:  # Move right
                        D_list[index] = max(D_list[index] - 1, 0)
                        D_list[index + 1] += 1
                elif index == 0:  # Left boundary (can only move right)
                    D_list[index] = max(D_list[index] - 1, 0)
                    D_list[index + 1] += 1
                elif index == self.SSA_M - 1:  # Right boundary (can only move left)
                    D_list[index] = max(D_list[index] - 1, 0)
                    D_list[index - 1] += 1

                elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:  # Production reaction
                    # print("Production reaction")
                    D_list[compartment_index] += 1
                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:  # Degradation reaction
                    # print("Degradation reaction")
                    D_list[compartment_index] = max(D_list[compartment_index] - 1, 0)
                elif index >= 3 * self.SSA_M and index <= 4 * self.SSA_M - 1:  # Conversion from continuous to discrete
                    # print("Conversion from continuous to discrete")
                    D_list[compartment_index] += 1
                    C_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1/self.h
                    C_list = np.maximum(C_list, 0)  # Ensure non-negativity for continuous list
                    
                elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M-1:  # Conversion from discrete to continuous
                 
                    D_list[compartment_index] = max(D_list[compartment_index] - 1, 0)
                    C_list[self.PDE_multiple*compartment_index:self.PDE_multiple*(compartment_index+1)] += 1/self.h
                  
                    # approx_mass_after = self.approximate_mass_left_hand(D_list, C_list)
        
                    # print(f"Change in mass: {change}"
                # Store the results for the current time step
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    D_grid[:, time_index] = D_list
                    C_grid[:, time_index] = C_list

                old_time = t  # Update old_time
                t += tau  # Update time by the time step
            else:  # Else we run the ODE step

                #C_list = self.RK4(C_list)
                C_list = self.crank_nicholson(C_list)
                C_list = np.maximum(C_list, 0)  # Ensure non-negativity after RK4 step
                t = td
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    C_grid[:, time_index] = C_list
                    D_grid[:, time_index] = D_list
           
        return D_grid, C_grid


    def run_simulation(self,number_of_repeats):
        """This will run the simulation with a total number of repeats"""
        C_initial, D_initial = self.create_initial_dataframe()
        D_average = np.zeros_like(D_initial)
        C_average = np.zeros_like(C_initial)
        filled_C_grid = np.zeros_like(C_initial)
        filled_D_grid = np.zeros_like(D_initial)

       
        for _ in tqdm(range(number_of_repeats),desc="Running the Hybrid simulations"):
            D_current, C_current = self.hybrid_simulation(deepcopy(D_initial), deepcopy(C_initial))
            D_average += D_current
            C_average += C_current

        filled_D_grid = D_average/number_of_repeats
        filled_C_grid = C_average/number_of_repeats


        combined_grid = np.zeros_like(filled_C_grid)
    

        for i in range(filled_D_grid.shape[1]):
            for j in range(filled_D_grid.shape[0]):
                start_index = j*self.PDE_multiple
                end_index = (j+1)*self.PDE_multiple

                combined_grid[start_index:end_index,i] = filled_C_grid[start_index:end_index,i]+(1/self.h)*filled_D_grid[j,i]

        combined_grid[-1,:] = combined_grid[-2,:]
       

        print("Simulation completed")
        return filled_D_grid, filled_C_grid, combined_grid
    

    def save_simulation_data(self,D_grid,C_grid,combined_grid, datadirectory='data'):

        if not os.path.exists(datadirectory):
            os.makedirs(datadirectory)
        
        params = {
            'domain_length': self.L,
            'compartment_number': self.SSA_M,
            'PDE_multiple': self.PDE_multiple,
            'total_time': self.total_time,
            'timestep': self.timestep,
            'threshold': self.threshold,
            'gamma': self.gamma,
            'deltax': self.deltax,
            'production_rate': self.production_rate,
            'degredation_rate': self.degredation_rate,
            'diffusion_rate': self.diffusion_rate,
            'threshold_conc': self.threshold_conc,
            'initial_SSA': self.SSA_initial.tolist(),
            'h': self.h,
        }
        np.savez(os.path.join(datadirectory, 'Hybrid_data'),
                 D_grid = D_grid,
                 C_grid = C_grid,
                 combined_grid = combined_grid,
                 time_vector = self.time_vector,
                 SSA_X = self.SSA_X,
                 PDE_X = self.PDE_X)
        
        with open(os.path.join(datadirectory, "parameters.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)

        print("Data saved successfully")