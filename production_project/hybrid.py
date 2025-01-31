import numpy as np
from tqdm import tqdm
import os
import json
from copy import deepcopy, copy
import ctypes
from .base_function import UtilityFunctions
from production_project.clibrary_argtypes import set_clibrary_argtypes #Each data type for the c functions
clibrary = ctypes.CDLL("c_class/clibrary.so") #import the c library

set_clibrary_argtypes(clibrary) #Import the data types for each c function

class Hybrid:
    
    def __init__(self, domain_length, compartment_number, PDE_multiple, total_time, timestep, threshold, gamma, production_rate, degradation_rate, diffusion_rate, SSA_initial,use_c_functions):
        self.L = domain_length
        self.SSA_M = compartment_number
        self.PDE_multiple = PDE_multiple
        self.production_rate = production_rate
        self.PDE_M = compartment_number * PDE_multiple
        self.deltax = self.L / self.PDE_M

        self.use_c_functions = use_c_functions #Whether use c_function or not 
        if self.use_c_functions:
            print("Using c functions")
        else: 
            print(f"Using python function")
        
        self.total_time = total_time
        self.timestep = timestep
        self.threshold = threshold
        self.gamma = gamma
        self.degradation_rate = degradation_rate
        self.h = self.L / compartment_number
        self.diffusion_rate = diffusion_rate
        self.d = diffusion_rate / (self.h**2)
        self.threshold_conc = threshold / self.h
        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M)
        self.PDE_X = np.linspace(0, self.L, self.PDE_M)

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA initial is not a np array")
        elif not len(SSA_initial) == compartment_number:
            raise ValueError("The length of the SSA initial is not the same as compartment number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("The SSA initial is not an integer")
        else:
            self.SSA_initial = SSA_initial.astype(int)

        self.PDE_initial_conditions = np.zeros_like(self.PDE_X, dtype=np.float64)
        self.steady_state = production_rate / degradation_rate
        self.DX_NEW = self.create_finite_difference()
        self.time_vector = np.arange(0, total_time, timestep)
        print("Successfully initialized the hybrid model")
        print(f"The threshold concentration is: {self.threshold_conc}")

    def create_finite_difference(self) -> np.ndarray:
        self.DX = np.zeros((self.PDE_M, self.PDE_M), dtype=int)
        self.DX[0, 0], self.DX[-1, -1] = -1, -1
        self.DX[0, 1], self.DX[-1, -2] = 1, 1
        for i in range(1, self.DX.shape[0] - 1):
            self.DX[i, i] = -2
            self.DX[i, (i + 1)] = 1
            self.DX[i, (i - 1)] = 1
        return self.DX
    
    def create_initial_dataframe(self) -> np.ndarray:
        SSA_grid = np.zeros((self.SSA_M, len(self.time_vector)), dtype=int)
        SSA_grid[:, 0] = self.SSA_initial
        PDE_grid = np.zeros((self.PDE_M, len(self.time_vector)), dtype=float)
        PDE_grid[:, 0] = self.PDE_initial_conditions
        return PDE_grid, SSA_grid 
    
    def calculate_total_mass(self, PDE_list: np.ndarray, SSA_list: np.ndarray) -> np.ndarray:
        """This will calculate the total mass of discrete + continuous"""

        return UtilityFunctions.calculate_total_mass(PDE_list, SSA_list, self.use_c_functions,self.PDE_multiple, self.deltax, self.SSA_M )
      
    def threshold_boolean(self, combined_list: np.ndarray) -> np.ndarray:
        """Generate a boolean list based on the threshold"""

        compartment_bool_list, PDE_bool_list =  UtilityFunctions.threshold_boolean(combined_list, self.threshold, self.PDE_multiple ,self.SSA_M)

        return compartment_bool_list, PDE_bool_list

    def boolean_if_less_mass(self, PDE_list: np.ndarray) -> np.ndarray: 

        return UtilityFunctions.boolean_if_less_mass(PDE_list, self.h, self.PDE_multiple, self.SSA_M)
        
    def RHS_derivative(self, old_vector, boolean_threshold, SSA_fine_mass):
        dudt = np.zeros_like(old_vector)
        nabla = self.DX_NEW
        
        bool_production = self.production_rate * boolean_threshold 
        dudt = self.diffusion_rate * (1 / self.deltax)**2 * nabla @ old_vector - self.degradation_rate * (old_vector ** 2) + bool_production * (old_vector + SSA_fine_mass)
        return dudt

    def fine_grid_SSA_mass(self, SSA_mass):
        """Convert the SSA_mass to the same fine resolution as the PDE"""
        return UtilityFunctions.fine_grid_SSA_mass(SSA_mass, self.PDE_X, self.SSA_M, self.PDE_multiple)
    

    def RK4(self, old_vector, boolean_threshold, SSA_fine_mass):
        k1 = self.RHS_derivative(old_vector, boolean_threshold, SSA_fine_mass)
        k2 = self.RHS_derivative(old_vector + 0.5 * self.timestep * k1, boolean_threshold, SSA_fine_mass)
        k3 = self.RHS_derivative(old_vector + 0.5 * self.timestep * k2, boolean_threshold, SSA_fine_mass)
        k4 = self.RHS_derivative(old_vector + self.timestep * k3, boolean_threshold, SSA_fine_mass)
        return old_vector + self.timestep * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
 
    def propensity_calculation(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        SSA_list = SSA_list.astype(int)
        PDE_list = PDE_list.astype(float)

        combined_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)
        boolean_SSA_threshold, boolean_PDE_threshold = self.threshold_boolean(combined_list)

        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]
        movement_propensity[-1] = self.d * SSA_list[-1]

        R1_propensity = self.production_rate * combined_list * boolean_SSA_threshold
        R2_propensity = self.degradation_rate * (1 / self.h) * SSA_list * (SSA_list - 1)
        R3_propensity = 2 * self.degradation_rate * (1 / self.h) * approximate_PDE_mass * SSA_list

        conversion_to_discrete = np.zeros_like(SSA_list)
        conversion_to_cont = np.zeros_like(approximate_PDE_mass)
        boolean_SSA_threshold = self.boolean_if_less_mass(PDE_list).astype(int)
        conversion_to_discrete[combined_list < self.threshold] = approximate_PDE_mass[combined_list < self.threshold] * self.gamma
        conversion_to_discrete *= boolean_SSA_threshold
        conversion_to_cont[combined_list >= self.threshold] = SSA_list[combined_list >= self.threshold] * self.gamma
        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity, R3_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity

    def check_negative_values(self, vector: np.ndarray, vector_name: str):
        """
        Checks if a vector has negative values and raises an error if any are found.

        Args:
            vector (np.ndarray): The input vector to check.
            vector_name (str): The name of the vector (for error message context).

        Raises:
            ValueError: If the vector contains negative values below the machine error threshold.
        """
        machine_error = 10e-5
        negative_indices = np.where(vector < -machine_error)[0]
        if negative_indices.size > 0:
            print(f"Negative values found at indices: {negative_indices}")
            print(vector)
            raise ValueError(f"The vector named '{vector_name}' has negative values.")
        return None

    def hybrid_simulation(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, approx_mass: np.ndarray) -> np.ndarray:
        t = 0
        old_time = t
        td = self.timestep
        PDE_particles = np.zeros_like(approx_mass)
        SSA_list = SSA_grid[:, 0].astype(int)
        PDE_list = PDE_grid[:, 0].astype(float)
        ind_after = 0
        while t < self.total_time:
            total_propensity = self.propensity_calculation(SSA_list, PDE_list)
            fine_SSA_mass = self.fine_grid_SSA_mass(SSA_list)
            combined_mass = self.calculate_total_mass(PDE_list, SSA_list)[0]
            SSA_boolean_threshold, PDE_boolean_threshold = self.threshold_boolean(combined_mass)
            
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:
                PDE_list = self.RK4(PDE_list, PDE_boolean_threshold, fine_SSA_mass)
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    self.check_negative_values(PDE_list, "PDE_list")
                    self.check_negative_values(SSA_grid, "SSA_list")
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)
                old_time = t 
                continue 

            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)
            alpha_cum = np.cumsum(total_propensity)
            index = np.searchsorted(alpha_cum, r2 * alpha0)
            compartment_index = index % self.SSA_M
            if t + tau <= td:
                if index <= self.SSA_M - 2 and index >= 1:
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
                elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:
                    SSA_list[compartment_index] += 1
                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:
                    SSA_list[compartment_index] -= 1
                elif index >= 3 * self.SSA_M and index <= 4 * self.SSA_M - 1:
                    SSA_list[compartment_index] -= 1
                elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M - 1:
                    SSA_list[compartment_index] += 1
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1 / self.h
                else:
                    SSA_list[compartment_index] -= 1 
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] += 1 / self.h
                t += tau
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    SSA_grid[:, time_index] = SSA_list
                    PDE_grid[:, time_index] = PDE_list
                    self.check_negative_values(PDE_list, "PDE_list")
                    self.check_negative_values(SSA_list, "SSA_list")
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)
                old_time = t  
            else:
                PDE_list = self.RK4(PDE_list, PDE_boolean_threshold, fine_SSA_mass)
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    self.check_negative_values(PDE_list, "PDE_list")
                    self.check_negative_values(SSA_list, "SSA_list")
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)
                old_time = t 
        return SSA_grid, PDE_grid, approx_mass

    def run_simulation(self, number_of_repeats: int) -> np.ndarray:
        PDE_initial, SSA_initial = self.create_initial_dataframe()
        approx_mass_initial = np.zeros_like(SSA_initial)
        approx_mass_initial[:, 0] = self.calculate_total_mass(PDE_initial[:, 0], SSA_initial[:, 0])[0]
        SSA_sum = np.zeros_like(SSA_initial)
        PDE_sum = np.zeros_like(PDE_initial)
        approx_mass_sum = np.zeros_like(approx_mass_initial)
        for _ in tqdm(range(number_of_repeats), desc="Running the Hybrid simulations"):
            SSA_current, PDE_current, approx_mass_current = self.hybrid_simulation(deepcopy(SSA_initial), deepcopy(PDE_initial), deepcopy(approx_mass_initial))
            SSA_sum += SSA_current
            PDE_sum += PDE_current
            approx_mass_sum += approx_mass_current
        SSA_average = SSA_sum / number_of_repeats
        PDE_average = PDE_sum / number_of_repeats
        approx_sum_average = approx_mass_sum / number_of_repeats
        combined_grid = np.zeros_like(PDE_average)
        for i in range(SSA_average.shape[1]):
            for j in range(SSA_average.shape[0]):
                start_index = j * self.PDE_multiple
                end_index = (j + 1) * self.PDE_multiple
                combined_grid[start_index:end_index, i] = PDE_average[start_index:end_index, i] + (1 / self.h) * SSA_average[j, i]
        combined_grid[-1, :] = combined_grid[-2, :]
        print("Simulation completed")
        return SSA_average, PDE_average, combined_grid

    def save_simulation_data(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, combined_grid: np.ndarray, datadirectory='data'):
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
            'degradation_rate': self.degradation_rate,
            'diffusion_rate': self.diffusion_rate,
            'threshold_conc': self.threshold_conc,
            'initial_SSA': self.SSA_initial.tolist(),
            'h': self.h,
        }
        np.savez(os.path.join(datadirectory, 'Hybrid_data'),
                 SSA_grid=SSA_grid,
                 PDE_grid=PDE_grid,
                 combined_grid=combined_grid,
                 time_vector=self.time_vector,
                 SSA_X=self.SSA_X,
                 PDE_X=self.PDE_X)
        with open(os.path.join(datadirectory, "parameters.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)
        print("Data saved successfully")