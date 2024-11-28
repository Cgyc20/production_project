import numpy as np
"""Created 1/10/2024 as new project file"""
from tqdm import tqdm
import os
import ctypes
import json
from copy import deepcopy, copy
from production_project.clibrary_argtypes import set_clibrary_argtypes #Each data type for the c functions


clibrary = ctypes.CDLL("c_class/clibrary.so") #import the c library

set_clibrary_argtypes(clibrary) #Import the data types for each c function


class Hybrid:
    
    def __init__(self, domain_length, compartment_number, PDE_multiple, total_time, timestep, threshold, gamma, production_rate, degradation_rate, diffusion_rate, SSA_initial, use_c_functions=False):
        self.L = domain_length #The domain length
        self.SSA_M = compartment_number #The number of compartments
        self.PDE_multiple = PDE_multiple #The number of PDE points per compartment
        self.production_rate = production_rate #The producti
        self.PDE_M = compartment_number * PDE_multiple+1
        self.deltax = self.L / (self.PDE_M-1) #The PDE stepsize
        self.total_time = total_time #total simulation time
        self.timestep = timestep #The timestep size
        self.threshold = threshold #THe threshold (which is per compartment, the number of cells 
        self.gamma = gamma #The rate of conversion
        self.degradation_rate = degradation_rate
        
        self.h = self.L / compartment_number #The size of each compartment
        self.diffusion_rate = diffusion_rate #Rate of diffusion
        self.production_rate_per_compartment = production_rate*self.h #THe 
        self.d = diffusion_rate / (self.h ** 2)  # The jump rate

        self.threshold_conc = threshold/self.h #The threshold concentration per cell
        
        self.SSA_X = np.linspace(0, self.L - self.h, self.SSA_M) #The spatial domain in which the SSA is defined on
        self.PDE_X = np.linspace(0, self.L, self.PDE_M) #The spatial domain that the PDE is defined on

        if not isinstance(SSA_initial, np.ndarray):
            raise ValueError("SSA initial is not a np array")
        elif not len(SSA_initial) == compartment_number:
            raise ValueError("The length of the SSA initial is not the same as compartment number")
        elif not np.issubdtype(SSA_initial.dtype, np.integer):
            raise ValueError("The SSA initial is not an integer")
        else:
            self.SSA_initial = SSA_initial.astype(int)

        # self.PDE_initial_conditions = np.zeros_like(self.PDE_X, dtype=np.float64) #Initially zero states for PDE
        self.PDE_initial_conditions = np.zeros_like(self.PDE_X, dtype=np.float64) #Initially zero states for PDE
        self.steady_state = production_rate / degradation_rate #Steady state of the system
        self.DX_NEW = self.create_finite_difference()  # Ensure DX_NEW is initialized here
        self.time_vector = np.arange(0, total_time, timestep)  # The time vector
        self.Crank_matrix, self.M1_inverse = self.create_crank_nicholson() #The crank method respective matrices

        self.use_c_functions = use_c_functions
        print("Successfully initialized the hybrid model")

        print(f"The threshold concentration is: {self.threshold_conc}")

    def create_crank_nicholson(self) -> np.ndarray:
        """Creates the matrix used for crank nicholson method """

        H = self.create_finite_difference()
        H = H/(self.deltax**2)
        
        M1  = np.identity(H.shape[0])-(0.5*self.timestep*self.diffusion_rate)*H
        M2 = np.identity(H.shape[0])+(0.5*self.timestep*self.diffusion_rate)*H
        # M1 = np.identity(H.shape[0]) *(1+0.5*self.timestep*self.degradation_rate) - 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H
        # M2 = np.identity(H.shape[0]) *(1-0.5*self.timestep*self.degradation_rate) + 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H
        M1_inverse = np.linalg.inv(M1)
        Crank_matrix = M1_inverse@M2
        
        return Crank_matrix, M1_inverse

    def create_finite_difference(self) -> np.ndarray:
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
    
    def create_initial_dataframe(self) -> np.ndarray:
        """Creates the intiial dataframes to be used throughout the simulation
        Returns: the initial dataframes for discrete and continuous numbers of molecules. With initial conditions"""
        SSA_grid = np.zeros((self.SSA_M,len(self.time_vector) ), dtype=int)  # Discrete molecules grid
        SSA_grid[:, 0] = self.SSA_initial
        PDE_grid = np.zeros((self.PDE_M, len(self.time_vector)), dtype=float)  # Continuous mass
        PDE_grid[:, 0] = self.PDE_initial_conditions

        return PDE_grid, SSA_grid 
        
    def crank_nicholson(self, old_vector: np.ndarray) -> np.ndarray:
        """Returns the new vector using the crank nicholson matrix"""
        old_vector = old_vector.astype(float)
        return self.Crank_matrix @ old_vector  # The new vector!

    def ApproximateLeftHandPython(self, PDE_list: np.ndarray) -> np.ndarray:
        """Works out the approximate mass of the PDE domain over each grid point. Using the left hand rule"""
        PDE_list = PDE_list.astype(float)
        approximation_number_cont = np.zeros(self.SSA_M)  # set empty list the length of number of compartments

        for i in range(self.SSA_M):  # Range over each compartment
            start_index = self.PDE_multiple * i  # start index of PDE point
            end_index = self.PDE_multiple * (i + 1)  # end index of PDE point
            sum_value = 0.0  # sum value as stands
            sum_value = np.sum(PDE_list[start_index:end_index]) * self.deltax  # Summing the list, left hand rule
            approximation_number_cont[i] = sum_value  # Adding to the list

        return approximation_number_cont
    
    def ApproximateLeftHandC(self, PDE_list):
        # Assuming the C library has been loaded and has a function `ApproxMassLeftHand`
        # which expects the arguments as C pointers to the arrays
        approximate_PDE_mass = np.zeros(self.SSA_M)
        PDE_list = np.array(PDE_list, dtype=np.float32)
        PDE_list_Ctypes = PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

        approximate_PDE_mass_Ctypes = approximate_PDE_mass.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        # Call the C function
        clibrary.ApproximateMassLeftHand(self.SSA_M, self.PDE_multiple, PDE_list_Ctypes, approximate_PDE_mass_Ctypes, self.deltax)

        # Convert the result back into a NumPy array
        approximate_PDE_mass = np.ctypeslib.as_array(approximate_PDE_mass_Ctypes, shape=approximate_PDE_mass.shape)
        return approximate_PDE_mass

    def calculate_total_mass(self, PDE_list: np.ndarray, SSA_list: np.ndarray) -> np.ndarray:
        """This will calculate the total mass of discrete + continuous"""
        PDE_list = PDE_list.astype(float)
        SSA_list = SSA_list.astype(int)
        if self.use_c_functions:
            mass_solver = self.ApproximateLeftHandC
        else:
            mass_solver = self.ApproximateLeftHandPython
        
        approximate_PDE_mass = mass_solver(PDE_list)
        combined_list = np.add(SSA_list, approximate_PDE_mass) 
        return combined_list, approximate_PDE_mass
    

    def booleanMassC(self, PDE_list: np.ndarray) -> np.ndarray:
        """Calculate boolean mass using c-type function"""
        PDE_list = PDE_list.astype(np.float32)
        boolean_PDE_list = np.zeros_like(PDE_list, dtype=np.int32)
        boolean_SSA_list = np.zeros(self.SSA_M, dtype=np.int32)

        PDE_list_Ctypes = PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        boolean_PDE_list_Ctypes = boolean_PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        boolean_SSA_list_Ctypes = boolean_SSA_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        clibrary.BooleanMass(self.SSA_M, self.PDE_M, self.PDE_multiple, PDE_list_Ctypes, boolean_PDE_list_Ctypes, boolean_SSA_list_Ctypes, self.h)

        boolean_SSA_list = np.ctypeslib.as_array(boolean_SSA_list_Ctypes, shape=boolean_SSA_list.shape)
        return boolean_SSA_list

    def booleanMassPython(self, PDE_list: np.ndarray) -> np.ndarray:
        """Calculate mass using python function"""
        PDE_list = PDE_list.astype(float)
        boolean_PDE_list = np.zeros_like(PDE_list)  # Set the boolean list to zero
        boolean_PDE_list[PDE_list > 1 / self.h] = 1  # s

        boolean_threshold_SSA = np.zeros(self.SSA_M)

        for i in range(self.SSA_M):
            """Run over the compartments"""
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
        """This takes in the PDE_list. If there is any instance in which a point is less than 1/h, then the boolean value will be 0
        Input: the PDE_list
        Returns: A boolean list with length SSA_M. 1 if all points are above 1/h (in that compartment), else 0 if there exists at least one point less than 1/h
        """
        if self.use_c_functions:
            return self.booleanMassC(PDE_list)
        else:
            return self.booleanMassPython(PDE_list)
        
            

    def propensity_calculationPython(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        """
        Calculates the propensity functions for each reaction.

        Args:
            SSA_list (np.ndarray): Discrete molecules list.
            PDE_list (np.ndarray): Continuous mass list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        SSA_list = SSA_list.astype(int)
        PDE_list = PDE_list.astype(float)
        
        movement_propensity = 2 * self.d * SSA_list  # The diffusion rates
        movement_propensity[0] = self.d * SSA_list[0]
        movement_propensity[-1] = self.d * SSA_list[-1]

        approximate_PDE_mass = np.zeros_like(SSA_list)
        combined_list = np.zeros_like(SSA_list)

        #approximate_PDE_mass = self.ApproximateLeftHandC(PDE_list)
        combined_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)
        
        conversion_to_discrete = np.zeros_like(SSA_list)  # length of SSA_m
        conversion_to_cont = np.zeros_like(approximate_PDE_mass)  # length as SSA_m 

        # Ensure the boolean index matches the array size
        if combined_list.shape != SSA_list.shape:
            raise ValueError("Shape mismatch between combined_list and SSA_list")

        boolean_SSA_threshold = self.boolean_if_less_mass(PDE_list).astype(int)

        conversion_to_discrete[combined_list < self.threshold] = approximate_PDE_mass[combined_list < self.threshold] * self.gamma
        conversion_to_discrete *= boolean_SSA_threshold
        conversion_to_cont[combined_list >= self.threshold] = SSA_list[combined_list >= self.threshold] * self.gamma
        
        combined_propensity = np.concatenate((movement_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity
    


        # Assuming the C library is already loaded (e.g., clibrary = ctypes.CDLL('./path_to_your_c_library.so'))

    def propensity_calculationC(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        """
        Calculates the propensity functions for each reaction using the C library.

        Args:
            SSA_list (np.ndarray): Discrete molecules list.
            PDE_list (np.ndarray): Continuous mass list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        SSA_list = np.ascontiguousarray(SSA_list, dtype=np.int32)
        PDE_list = np.ascontiguousarray(PDE_list, dtype=np.float32)
        

        # Initialize propensity_list and other arrays
        propensity_list = np.zeros(5 * self.SSA_M, dtype=np.float32)
        combined_mass_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)
        boolean_SSA_threshold = self.boolean_if_less_mass(PDE_list).astype(int)
        boolean_SSA_threshold = np.ascontiguousarray(boolean_SSA_threshold, dtype=np.int32)
        combined_mass_list = np.ascontiguousarray(combined_mass_list, dtype=np.float32)
        approximate_PDE_mass = np.ascontiguousarray(approximate_PDE_mass, dtype=np.float32)

      
        # Call the C function to calculate propensities
        clibrary.CalculatePropensity(
            self.SSA_M,
            PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            SSA_list.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            propensity_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            combined_mass_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            approximate_PDE_mass.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            boolean_SSA_threshold.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            self.degradation_rate,
            self.threshold,
            self.production_rate_per_compartment,
            self.gamma,
            self.d
        )

        # Convert propensity list back into numpy floats
        propensity_list = np.ctypeslib.as_array(propensity_list, shape=propensity_list.shape)
        return propensity_list
        
    

    def propensity_calculation(self, SSA_list: np.ndarray, PDE_list: np.ndarray) -> np.ndarray:
        """
        Wrapper function to choose between Python and C implementation of propensity calculation.

        Args:
            SSA_list (np.ndarray): Discrete molecules list.
            PDE_list (np.ndarray): Continuous mass list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        if self.use_c_functions:
            return self.propensity_calculationC(SSA_list, PDE_list)
        else:
            return self.propensity_calculationPython(SSA_list, PDE_list)

        
        # print(f"Propensity C: {propensity_c_list}")
        # print(f"Propensity Python: {propensity_python_list}")

        return propensity_python_list

    def hybrid_simulation(self, SSA_grid: np.ndarray, PDE_grid: np.ndarray, approx_mass: np.ndarray) -> np.ndarray:
        t = 0
        old_time = t
        td = self.timestep
        PDE_particles = np.zeros_like(approx_mass)
        SSA_list = SSA_grid[:, 0].astype(int)  # Starting SSA_list
        PDE_list = PDE_grid[:, 0].astype(float)  # Starting PDE_list
        ind_after = 0
        while t < self.total_time:
            total_propensity = self.propensity_calculation(SSA_list, PDE_list)
            alpha0 = np.sum(total_propensity)
            if alpha0 == 0:  # If no reactions can occur, execute PDE step
                PDE_list = self.crank_nicholson(PDE_list)
                PDE_list = np.maximum(PDE_list, 0)  # Ensure non-negativity after RK4 step
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list)
                old_time = t  # Update old_time
                continue  # Skip the rest of the loop and continue with the next iteration
    
            r1, r2, r3 = np.random.rand(3)
            tau = (1 / alpha0) * np.log(1 / r1)  # Time until next reaction
            alpha_cum = np.cumsum(total_propensity)  # Cumulative sum of propensities
            index = np.searchsorted(alpha_cum, r2 * alpha0)  # Determine which reaction occurs
    
            compartment_index = index % self.SSA_M  # The compartmental index is just the modulo of SSA. 
            if t + tau <= td:  # Execute Gillespie
                """The diffusion reactions are executed here"""
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
    
                
                    """Finally the conversion reactions here"""
                elif index >=  self.SSA_M and index <= 2 * self.SSA_M - 1:  # Conversion from continuous to discrete
                    SSA_list[compartment_index] += 1
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1 / self.h
                    PDE_list = np.maximum(PDE_list, 0)  # Ensure non-negativity for continuous list (probably don't need)
                
                #elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M-1:  # Conversion from discrete to continuous
                    # print(f"*"*30)
                    # print(f"Checking conversion to PDE, given this occurs")
                    # print(f"  {SSA_list}")
                    # print(f"Continuous mass at time {t:.1f}:")
                    # print(f"  {PDE_list.round(1)}")
                    # print(f"Number of particles continuous")
                    # print(f" {PDE_particles[:,min(ind_after+1, len(self.time_vector))-1]}")
                    # print(f"*"*30)
    
                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:  # Conversion from discrete to continuous
    
                    SSA_list[compartment_index] = max(SSA_list[compartment_index] - 1, 0)
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] += 1 / self.h
                    
                t += tau 
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    SSA_grid[:, time_index] = SSA_list
                    PDE_grid[:, time_index] = PDE_list
                    approx_mass[:, time_index], PDE_particles[:, time_index] = self.calculate_total_mass(PDE_list, SSA_list,)
    
                old_time = t  # Update old_time
                    # Update time by the time step
                #Printing a barrier
    
                
                # print(f"{'Simulation Step':^30}")  # Centered title within the asterisks
                # print("*" * 30)
    
                # # Time and mass information
                # print(f"Time: {t:.2f}")
                # print(f"Mass conversion threshold: {self.threshold}")
                # print("-" * 30)  # Separator line
    
                # # Particle and mass details
                # print(f"Stochastic particles in each box at time {t}:")
                # print(f"  {SSA_list}")
                # print(f"Continuous mass at time {t:.1f}:")
                # print(f"  {PDE_list.round(1)}")
                # print(f"Number of particles continuous")
                # print(f" {PDE_particles[:,min(ind_after+1, len(self.time_vector))-1]}")
                # print(f"Approximate mass at time {t:.1f}:")
                # print(f"  {approx_mass[:, min(ind_after+1, len(self.time_vector))-1]}")
                # print("-" * 30)
    
                # # Propensity information
                # print(f"{'Propensity Details':^30}")
                # print(f"Index of reaction chosen: {index}")
                # print("-" * 30)
                # print(f"Movement propensity:           {total_propensity[:self.SSA_M]}")
                # print(f"Production propensity:         {total_propensity[self.SSA_M:2*self.SSA_M]}")
                # print(f"Degradation propensity:        {total_propensity[2*self.SSA_M:3*self.SSA_M]}")
                # print(f"Conversion to discrete prop.:  {total_propensity[3*self.SSA_M:4*self.SSA_M]}")
                # print(f"Conversion to continuous prop.: {total_propensity[4*self.SSA_M:]}")
                # print("*" * 30)
                # print("\n")  # Extra blank line for space between steps
    
    
    
            else:  # Else we run the ODE step
                PDE_list = self.crank_nicholson(PDE_list)
                PDE_list = np.maximum(PDE_list, 0)  # Ensure non-negativity after RK4 step
                t = copy(td)
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    approx_mass[:, time_index], PDE_particles[:, time_index]= self.calculate_total_mass(PDE_list, SSA_list)
            
        return SSA_grid, PDE_grid, approx_mass
    

    def run_simulation(self, number_of_repeats: int) -> np.ndarray:
        """This will run the simulation with a total number of repeats"""
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
            approx_mass_sum += approx_mass_sum
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

    def save_simulation_data(self,SSA_grid:np.ndarray,PDE_grid:np.ndarray,combined_grid:np.ndarray, datadirectory='data'):

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
                 SSA_grid = SSA_grid,
                 PDE_grid = PDE_grid,
                 combined_grid = combined_grid,
                 time_vector = self.time_vector,
                 SSA_X = self.SSA_X,
                 PDE_X = self.PDE_X)
        
        with open(os.path.join(datadirectory, "parameters.json"), 'w') as params_file:
            json.dump(params, params_file, indent=4)

        print("Data saved successfully")