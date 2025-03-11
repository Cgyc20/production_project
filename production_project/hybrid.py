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
        
    # def RHS_derivative(self, old_vector, boolean_threshold, SSA_fine_mass):
    #     dudt = np.zeros_like(old_vector)
    #     nabla = self.DX_NEW
        
    #     bool_production = self.production_rate * boolean_threshold 
    #     dudt = self.diffusion_rate * (1 / self.deltax)**2 * nabla @ old_vector - self.degradation_rate * (old_vector ** 2) + bool_production * (old_vector+SSA_fine_mass)
    #     #dudt = self.diffusion_rate * (1 / self.deltax)**2 * nabla @ old_vector - self.degradation_rate * (old_vector ** 2) + self.production_rate* (old_vector)
    #     return dudt

    def RHS_derivative(self, old_vector, boolean_threshold, SSA_fine_mass):
        nabla = self.DX_NEW
        diff_coeff = self.diffusion_rate * (1 / self.deltax) ** 2

        # Precompute terms
        diffusion_term = diff_coeff * (nabla @ old_vector)
        degradation_term = self.degradation_rate * boolean_threshold*((old_vector+SSA_fine_mass) ** 2)
        production_term = self.production_rate * boolean_threshold * (old_vector + SSA_fine_mass)
        degradation_term = self.degradation_rate * (old_vector) ** 2
        # Combine all terms
        dudt = diffusion_term - degradation_term + production_term

        return dudt

    def fine_grid_SSA_mass(self, SSA_mass):
        """Convert the SSA_mass to the same fine resolution as the PDE"""
        return UtilityFunctions.fine_grid_SSA_mass(SSA_mass, self.PDE_X, self.SSA_M, self.PDE_multiple, self.h)
    
    def RK4(self, old_vector, boolean_threshold, SSA_fine_mass, dt=None):

        if dt == None:
            dt = self.timestep
        
        k1 = self.RHS_derivative(old_vector, boolean_threshold, SSA_fine_mass)
        k2 = self.RHS_derivative(old_vector + 0.5 * dt * k1, boolean_threshold, SSA_fine_mass)
        k3 = self.RHS_derivative(old_vector + 0.5 * dt * k2, boolean_threshold, SSA_fine_mass)
        k4 = self.RHS_derivative(old_vector + dt * k3, boolean_threshold, SSA_fine_mass)
        return old_vector + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    

 
    def propensity_calculation(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        SSA_list = SSA_list.astype(int)
        PDE_list = PDE_list.astype(float)

        combined_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)
        boolean_SSA_threshold, _ = self.threshold_boolean(combined_list)

        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]
        movement_propensity[-1] = self.d * SSA_list[-1]

        R1_propensity = self.production_rate * combined_list * boolean_SSA_threshold #D-> 2D
        #R1_propensity = self.production_rate * SSA_list #D-> 2A
        R2_propensity = self.degradation_rate * (1 / self.h) * SSA_list * (SSA_list - 1) #2D -> D
        R3_propensity = 2 * self.degradation_rate * (1 / self.h) * approximate_PDE_mass * SSA_list # D+C -> C

        conversion_to_discrete = np.zeros_like(SSA_list)
        conversion_to_cont = np.zeros_like(approximate_PDE_mass)
        boolean_SSA_threshold = self.boolean_if_less_mass(PDE_list).astype(int)
        conversion_to_discrete[combined_list < self.threshold] = approximate_PDE_mass[combined_list < self.threshold] * self.gamma
        conversion_to_discrete *= boolean_SSA_threshold
        conversion_to_cont[combined_list >= self.threshold] = SSA_list[combined_list >= self.threshold] * self.gamma
        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity, R3_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity
    

    
    def propensity_calculationC(self, SSA_list, PDE_list, combined_mass_list, Approximate_PDE_Mass, boolean_SSA_list):

        # Convert inputs to appropriate data types
        SSA_list = np.ascontiguousarray(SSA_list, dtype=np.int32)
        PDE_list = np.ascontiguousarray(PDE_list, dtype=np.float32)
        combined_mass_list = np.ascontiguousarray(combined_mass_list, dtype=np.float32)
        Approximate_PDE_Mass = np.ascontiguousarray(Approximate_PDE_Mass, dtype=np.float32)
        boolean_SSA_list = np.ascontiguousarray(boolean_SSA_list, dtype=np.int32)
        boolean_if_less_mass = self.boolean_if_less_mass(PDE_list).astype(int)
        boolean_if_less_mass = np.ascontiguousarray(boolean_if_less_mass, dtype=np.int32)
                # Constants

       

        SSA_M = ctypes.c_int(len(SSA_list))
        degradation_rate_h = ctypes.c_float(self.degradation_rate)
        threshold = ctypes.c_float(self.threshold_conc)
        production_rate = ctypes.c_float(self.production_rate)
        gamma = ctypes.c_float(self.gamma)
        jump_rate = ctypes.c_float(self.d)


        # Output array for propensities
        propensity_list = np.zeros(6 * SSA_M.value, dtype=np.float32)

        # Convert NumPy arrays to ctypes pointers
        PDE_list_ptr = PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        SSA_list_ptr = SSA_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        propensity_list_ptr = propensity_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        boolean_SSA_list_ptr = boolean_SSA_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        combined_mass_list_ptr = combined_mass_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        Approximate_PDE_Mass_ptr = Approximate_PDE_Mass.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        boolean_if_less_mass_ptr = boolean_if_less_mass.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        
        # Call the C function
        clibrary.CalculatePropensity(
            SSA_M,
            PDE_list_ptr,
            SSA_list_ptr,
            propensity_list_ptr,
            boolean_SSA_list_ptr,
            combined_mass_list_ptr,
            Approximate_PDE_Mass_ptr,
            boolean_if_less_mass_ptr,
            degradation_rate_h,
            threshold,
            production_rate,
            gamma,
            jump_rate
        )
        return propensity_list



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

    def test_boolean(self, combined_mass: np.ndarray, compartment_boolean_input: np.ndarray, PDE_boolean_input: np.ndarray):
        """THis will test if the Boolean is doing as expected."""

        compartment_boolean_internal = combined_mass[combined_mass<= self.threshold]
        PDE_boolean_internal = np.zeros_like(PDE_boolean_input)
        # print(f"The compartment_boolean shape: {len(compartment_boolean_input)}")
        for i in range(min(self.SSA_M, len(compartment_boolean_internal))):
            value = compartment_boolean_internal[i]
            new_value = 0 if value == 1 else 1
            start_index = i * self.PDE_multiple
            PDE_boolean_internal[start_index:start_index + self.PDE_multiple] = new_value
                    

        if compartment_boolean_internal.any() != compartment_boolean_input.any():
            raise ValueError("The compartment booleaan input does not match current combined_mass")
        if PDE_boolean_input.any() != PDE_boolean_input.any():
            raise ValueError("The PDE_boolean_threshold input does not match current combined mass")
        
        return None

    def hybrid_simulation(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, approx_mass: np.ndarray) -> np.ndarray:
        t = 0
        old_time = t
        td = self.timestep
        PDE_particles = np.zeros_like(approx_mass)
        SSA_list = SSA_grid[:, 0].astype(int)
        PDE_list = PDE_grid[:, 0].astype(float)
        ind_after = 0

        # Arrays to track SSA events and PDE updates
        SSA_events_log = []  # Format: (time, compartment_index, reaction_type)
        PDE_update_times = []  # List of times when PDE is updated

        while t < self.total_time:
            #total_propensity = self.propensity_calculation(SSA_list, PDE_list)
            
            fine_SSA_mass = self.fine_grid_SSA_mass(SSA_list)
            combined_mass, approx_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)
            SSA_boolean_list, PDE_boolean_threshold = self.threshold_boolean(combined_mass)
            total_propensity = self.propensity_calculation(SSA_list, PDE_list)
            #total_propensity = self.propensity_calculationC(SSA_list, PDE_list, combined_mass, approx_PDE_mass,SSA_boolean_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:
                PDE_list = self.RK4(PDE_list, PDE_boolean_threshold, fine_SSA_mass)
                PDE_update_times.append(t)

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
                continue

            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)
            alpha_cum = np.cumsum(total_propensity)
            index = np.searchsorted(alpha_cum, r2 * alpha0)
            compartment_index = index % self.SSA_M

            if t + tau <= td:
                reaction_type = ""
                if index <= self.SSA_M - 2 and index >= 1:  # Diffusion
                    if r3 < 0.5:
                        SSA_list[index] -= 1
                        SSA_list[index - 1] += 1
                    
                        reaction_type = "diffusion"

                    else:
                        SSA_list[index] -= 1
                        SSA_list[index + 1] += 1
                        reaction_type = "diffusion"
                elif index == 0:
                    SSA_list[index] -= 1
                    SSA_list[index + 1] += 1
                    reaction_type = "diffusion"
                elif index == self.SSA_M - 1:
                    SSA_list[index] -= 1
                    SSA_list[index - 1] += 1
                    reaction_type = "diffusion"
                elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:
                    SSA_list[compartment_index] += 1  # D -> 2D
                    reaction_type = "D duplication"
                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:
                    SSA_list[compartment_index] -= 1  # 2D -> D
                    reaction_type = "D degradation"
                elif index >= 3 * self.SSA_M and index <= 4 * self.SSA_M - 1:
                    SSA_list[compartment_index] -= 1  # D + C -> C
                    reaction_type = "J degredation"
                elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M - 1:  # C -> D
                    SSA_list[compartment_index] += 1
                    PDE_list[self.PDE_multiple * compartment_index: self.PDE_multiple * (compartment_index + 1)] -= 1 / self.h
                    reaction_type = "conversion_C_to_D"
                else:
                    SSA_list[compartment_index] -= 1  # D -> C
                    PDE_list[self.PDE_multiple * compartment_index: self.PDE_multiple * (compartment_index + 1)] += 1 / self.h
                    reaction_type = "conversion_D_to_C"
                
                SSA_events_log.append((t, compartment_index, reaction_type))

                # PDE_list = self.RK4(PDE_list, PDE_boolean_threshold, fine_SSA_mass,tau)
                t += tau
                old_time = t


            else:
                PDE_list = self.RK4(PDE_list, PDE_boolean_threshold, fine_SSA_mass)
                PDE_update_times.append(t)

                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                # print(f"Ind before in PDE {ind_before}")
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                # print(f"ind_after - ind_after in PDE {ind_after-ind_before}")
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    self.check_negative_values(PDE_list, "PDE_list")
                    self.check_negative_values(SSA_list, "SSA_list")
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)

                old_time = t

        return SSA_grid, PDE_grid, approx_mass, SSA_events_log, PDE_update_times

    def run_simulation(self, number_of_repeats: int) -> np.ndarray:
        PDE_initial, SSA_initial = self.create_initial_dataframe()
        approx_mass_initial = np.zeros_like(SSA_initial)
        approx_mass_initial[:, 0] = self.calculate_total_mass(PDE_initial[:, 0], SSA_initial[:, 0])[0]
        SSA_sum = np.zeros_like(SSA_initial)
        PDE_sum = np.zeros_like(PDE_initial)
        approx_mass_sum = np.zeros_like(approx_mass_initial)

        # Arrays to track all SSA events and PDE update times across repeats
        all_SSA_events_logs = []
        all_PDE_update_times = []

        for _ in tqdm(range(number_of_repeats), desc="Running the Hybrid simulations"):
            SSA_current, PDE_current, approx_mass_current, SSA_events_log, PDE_update_times = self.hybrid_simulation(
                deepcopy(SSA_initial), deepcopy(PDE_initial), deepcopy(approx_mass_initial))
            SSA_sum += SSA_current
            PDE_sum += PDE_current
            approx_mass_sum += approx_mass_current

            all_SSA_events_logs.append(SSA_events_log)
            all_PDE_update_times.append(PDE_update_times)

        SSA_average = SSA_sum / number_of_repeats
        PDE_average = PDE_sum / number_of_repeats

        combined_grid = np.zeros_like(PDE_average)
        for i in range(SSA_average.shape[1]):
            for j in range(SSA_average.shape[0]):
                start_index = j * self.PDE_multiple
                end_index = (j + 1) * self.PDE_multiple
                combined_grid[start_index:end_index, i] = PDE_average[start_index:end_index, i] + (1 / self.h) * SSA_average[j, i]
        combined_grid[-1, :] = combined_grid[-2, :]

        print("Simulation completed")

        # Save simulation data including SSA events and PDE update times

        return SSA_average, PDE_average, combined_grid, all_SSA_events_logs, all_PDE_update_times

     
    def save_simulation_data(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, combined_grid: np.ndarray, all_SSA_events_logs: list, all_PDE_update_times: list, datadirectory='data'):
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
        
        # Save the lists of SSA events and PDE update times separately
        np.save(os.path.join(datadirectory, 'SSA_events_logs.npy'), np.array(all_SSA_events_logs, dtype=object))
        np.save(os.path.join(datadirectory, 'PDE_update_times.npy'), np.array(all_PDE_update_times, dtype=object))
        
        with open(os.path.join(datadirectory, "parameters.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)
        print("Data saved successfully")