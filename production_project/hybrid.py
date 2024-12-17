import numpy as np
from tqdm import tqdm
import os
import json
from copy import deepcopy, copy

class Hybrid:
    
    def __init__(self, domain_length, compartment_number, PDE_multiple, total_time, timestep, threshold, gamma, production_rate, degradation_rate, diffusion_rate, SSA_initial, use_c_functions=False):
        self.L = domain_length
        self.SSA_M = compartment_number
        self.PDE_multiple = PDE_multiple
        self.production_rate = production_rate
        self.PDE_M = compartment_number * PDE_multiple
        self.deltax = self.L / (self.PDE_M)
        
        self.total_time = total_time
        self.timestep = timestep
        self.threshold = threshold
        self.gamma = gamma
        self.degradation_rate = degradation_rate
        self.h = self.L / compartment_number
        self.diffusion_rate = diffusion_rate
        self.production_rate_per_compartment = production_rate * self.h
        self.degradation_rate_per_compartment = degradation_rate * self.h
        # self.d = diffusion_rate / (self.h ** 2)
        self.d = diffusion_rate / (self.deltax**2*self.h ** 2)
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
        self.use_c_functions = use_c_functions
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
        
    def RHS_derivative(self, old_vector):
        dudt = np.zeros_like(old_vector)
        nabla = self.DX_NEW
        dudt = self.diffusion_rate*(1/self.deltax)**2 * nabla @ old_vector - self.degradation_rate * old_vector ** 2 + self.production_rate * old_vector
        return dudt
    
    def RK4(self, old_vector):
        k1 = self.RHS_derivative(old_vector)
        k2 = self.RHS_derivative(old_vector + 0.5 * self.timestep * k1)
        k3 = self.RHS_derivative(old_vector + 0.5 * self.timestep * k2)
        k4 = self.RHS_derivative(old_vector + self.timestep * k3)
        return old_vector + self.timestep * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    def ApproximateLeftHandPython(self, PDE_list: np.ndarray) -> np.ndarray:
        PDE_list = PDE_list.astype(float)
        approximation_number_cont = np.zeros(self.SSA_M)
        for i in range(self.SSA_M):
            start_index = self.PDE_multiple * i
            end_index = self.PDE_multiple * (i + 1)
            sum_value = np.sum(PDE_list[start_index:end_index]) * self.deltax
            approximation_number_cont[i] = sum_value
        return approximation_number_cont

    def calculate_total_mass(self, PDE_list: np.ndarray, SSA_list: np.ndarray) -> np.ndarray:
        PDE_list = PDE_list.astype(float)
        SSA_list = SSA_list.astype(int)
        approximate_PDE_mass = self.ApproximateLeftHandPython(PDE_list)
        combined_list = np.add(SSA_list, approximate_PDE_mass)
        return combined_list, approximate_PDE_mass

    def booleanMassPython(self, PDE_list: np.ndarray) -> np.ndarray:
        PDE_list = PDE_list.astype(float)
        boolean_PDE_list = np.zeros_like(PDE_list)
        boolean_PDE_list[PDE_list > 1 / self.h] = 1
        boolean_threshold_SSA = np.zeros(self.SSA_M)
        for i in range(self.SSA_M):
            start_index = i * self.PDE_multiple
            BOOL_VALUE = True
            for j in range(self.PDE_multiple):
                current_index = start_index + j
                if boolean_PDE_list[j] == 0:
                    BOOL_VALUE = False
            if BOOL_VALUE:
                boolean_threshold_SSA[i] = 1 
            else:
                boolean_threshold_SSA[i] = 0
        return boolean_threshold_SSA

    def boolean_if_less_mass(self, PDE_list: np.ndarray) -> np.ndarray:
        return self.booleanMassPython(PDE_list)
        
    def propensity_calculationPython(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        SSA_list = SSA_list.astype(int)
        PDE_list = PDE_list.astype(float)
        movement_propensity = 2 * self.d * SSA_list
        movement_propensity[0] = self.d * SSA_list[0]
        movement_propensity[-1] = self.d * SSA_list[-1]


        R1_propensity = self.production_rate_per_compartment * SSA_list
        R2_propensity = self.degradation_rate * SSA_list * (SSA_list - 1)


        combined_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)

        R3_propensity = self.degradation_rate * approximate_PDE_mass * SSA_list

        conversion_to_discrete = np.zeros_like(SSA_list)
        conversion_to_cont = np.zeros_like(approximate_PDE_mass)
        boolean_SSA_threshold = self.boolean_if_less_mass(PDE_list).astype(int)
        conversion_to_discrete[combined_list < self.threshold] = approximate_PDE_mass[combined_list < self.threshold] * self.gamma
        conversion_to_discrete *= boolean_SSA_threshold
        conversion_to_cont[combined_list >= self.threshold] = SSA_list[combined_list >= self.threshold] * self.gamma
        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity, R3_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity

    def propensity_calculation(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        return self.propensity_calculationPython(SSA_list, PDE_list)

    def hybrid_simulation(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, approx_mass: np.ndarray) -> np.ndarray:
        t = 0
        old_time = t
        td = self.timestep
        PDE_particles = np.zeros_like(approx_mass)
        SSA_list = SSA_grid[:, 0].astype(int)
        PDE_list = PDE_grid[:, 0].astype(float)
        ind_after = 0
        while t < self.total_time:
            total_propensity = self.propensity_calculationPython(SSA_list, PDE_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:
                PDE_list = self.RK4(PDE_list)
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
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
                        SSA_list[index] = SSA_list[index] - 1
                        SSA_list[index - 1] += 1
                    else:
                        SSA_list[index] = SSA_list[index] - 1
                        SSA_list[index + 1] += 1
                elif index == 0:
                    SSA_list[index] = SSA_list[index] - 1
                    SSA_list[index + 1] += 1
                elif index == self.SSA_M - 1:
                    SSA_list[index] = SSA_list[index] - 1
                    SSA_list[index - 1] += 1
                elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1: #D -> D+D
                    SSA_list[compartment_index] += 1
                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1: # D+D -> D
                    SSA_list[compartment_index] -= 1

                elif index >= 3 * self.SSA_M and index <= 4 * self.SSA_M - 1: #D+C -> D
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1 / self.h
                elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M - 1: # C -> D
                    SSA_list[compartment_index] += 1
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1 / self.h
                else: #D-> C
                    SSA_list[compartment_index] -= 1
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] += 1 / self.h
                t += tau 
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    SSA_grid[:, time_index] = SSA_list
                    PDE_grid[:, time_index] = PDE_list
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)
                old_time = t  
            else:
                PDE_list = self.RK4(PDE_list)
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
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