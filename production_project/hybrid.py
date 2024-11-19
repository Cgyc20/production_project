import numpy as np
"""Created 1/10/2024 as new project file"""
from tqdm import tqdm
import os
import json
from copy import deepcopy

class Hybrid:
    
    def __init__(self, domain_length, compartment_number, PDE_multiple, total_time, timestep, threshold, gamma, production_rate, degredation_rate, diffusion_rate, SSA_initial):
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
        self.degredation_rate = degredation_rate
        
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
            self.SSA_initial = SSA_initial

        # self.PDE_initial_conditions = np.zeros_like(self.PDE_X, dtype=np.float64) #Initially zero states for PDE
        self.PDE_initial_conditions = np.ones_like(self.PDE_X, dtype=np.float64) #Initially zero states for PDE
        self.steady_state = production_rate / degredation_rate #Steady state of the system
        self.DX_NEW = self.create_finite_difference()  # Ensure DX_NEW is initialized here
        self.time_vector = np.arange(0, total_time, timestep)  # The time vector
        self.Crank_matrix, self.M1_inverse = self.create_crank_nicholson() #The crank method respective matrices

        print("Successfully initialized the hybrid model")

        print(f"The threshold concentration is: {self.threshold_conc}")
    def create_crank_nicholson(self) -> np.ndarray:
        """Creates the matrix used for crank nicholson method """

        H = self.create_finite_difference()
        M1 = np.identity(H.shape[0]) *(1+0.5*self.timestep*self.degredation_rate) - 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H
        M2 = np.identity(H.shape[0]) *(1-0.5*self.timestep*self.degredation_rate) + 0.5*(self.timestep*self.diffusion_rate/self.deltax**2)*H
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
        SSA_grid = np.zeros((self.SSA_M,len(self.time_vector) ))  # Discrete molecules grid
        SSA_grid[:, 0] = self.SSA_initial
        PDE_grid = np.zeros((self.PDE_M, len(self.time_vector)))  # Continuous mass
        PDE_grid[:, 0] = self.PDE_initial_conditions

        return PDE_grid,SSA_grid 
        
    def crank_nicholson(self,old_vector: np.ndarray) -> np.ndarray:
        """Returns the new vector using the crank nicholson matrix"""

        return self.Crank_matrix@old_vector #The new vector!

    def approximate_mass_left_hand(self, PDE_list: np.ndarray) -> np.ndarray:
        """Works out the approximate mass of the PDE domain over each grid point. Using the left hand rule"""
        approximation_number_cont = np.zeros(self.SSA_M) #set empty list the length of number of compartments

        for i in range(self.SSA_M): #Range over each compartment
            start_index = self.PDE_multiple * i #start index of PDE point
            end_index = self.PDE_multiple * (i + 1) #end index of PDE point
            sum_value = 0.0 #sum value as stands
            sum_value = np.sum(PDE_list[start_index:end_index])*self.deltax  #Summing the list, left hand rule
            approximation_number_cont[i] = sum_value #Adding to the list

        return approximation_number_cont


    def calculate_total_mass(self,PDE_list:np.ndarray, SSA_list: np.ndarray) -> np.ndarray:
        """This will calculate the total mass of discrete + continuous"""

        approximate_PDE_mass = self.approximate_mass_left_hand(PDE_list)
        combined_list = np.add(SSA_list, approximate_PDE_mass) 
        return combined_list, approximate_PDE_mass
    

    def propensity_calculation(self, SSA_list: np.ndarray, PDE_list : np.ndarray) -> np.ndarray:
        """
        Calculates the propensity functions for each reaction.

        Args:
            SSA_list (np.ndarray): Discrete molecules list.
            PDE_list (np.ndarray): Continuous mass list.

        Returns:
            np.ndarray: Combined propensity list.
        """
        movement_propensity = 2 * self.d * SSA_list #The diffusion rates
        movement_propensity[0] = self.d * SSA_list[0]
        movement_propensity[-1] = self.d * SSA_list[-1]
    

        R1_propensity = self.production_rate_per_compartment * np.ones_like(SSA_list) #The production propensity
        R2_propensity = self.degredation_rate * SSA_list #degredation propensity


        combined_list, approximate_PDE_mass = self.calculate_total_mass(PDE_list, SSA_list)  #Add with the discrete mass to gather 
    
        conversion_to_discrete = np.zeros_like(SSA_list)
        conversion_to_cont = np.zeros_like(approximate_PDE_mass)

        # Ensure the boolean index matches the array size
        if combined_list.shape != SSA_list.shape:
            raise ValueError("Shape mismatch between combined_list and SSA_list")
    
        conversion_to_cont[combined_list >= self.threshold] = SSA_list[combined_list >= self.threshold] * self.gamma
        conversion_to_discrete[combined_list < self.threshold] = approximate_PDE_mass[combined_list < self.threshold] * self.gamma

        combined_propensity = np.concatenate((movement_propensity, R1_propensity, R2_propensity, conversion_to_discrete, conversion_to_cont))
        return combined_propensity

    def hybrid_simulation(self, SSA_grid : np.ndarray, PDE_grid : np.ndarray, approx_mass: np.ndarray) -> np.ndarray:
        t = 0
        old_time = t
        td = self.timestep
        PDE_particles = np.zeros_like(approx_mass)
        SSA_list = SSA_grid[:, 0]  # Starting SSA_list
        PDE_list = PDE_grid[:, 0]  # Starting PDE_list
        
        while t < self.total_time:
            total_propensity = self.propensity_calculation(SSA_list, PDE_list)
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
                """The diffusion reactions are executed here"""
                if index <= self.SSA_M - 2 and index >= 1:
                    # print("Reaction diffusion")
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



                    """Now the reaction kinetics"""

                elif index >= self.SSA_M and index <= 2 * self.SSA_M - 1:  # Production reaction
                    SSA_list[compartment_index] += 1

                elif index >= 2 * self.SSA_M and index <= 3 * self.SSA_M - 1:  # Degradation reaction
                    SSA_list[compartment_index] = max(SSA_list[compartment_index] - 1, 0)


                    """Finally the conversion reactions here"""
                elif index >= 3 * self.SSA_M and index <= 4 * self.SSA_M - 1:  # Conversion from continuous to discrete
                    # print("Conversion from continuous to discrete")
                    SSA_list[compartment_index] += 1
                    PDE_list[self.PDE_multiple * compartment_index : self.PDE_multiple * (compartment_index + 1)] -= 1/self.h
                    PDE_list = np.maximum(PDE_list, 0)  # Ensure non-negativity for continuous list
                    
                elif index >= 4 * self.SSA_M and index <= 5 * self.SSA_M-1:  # Conversion from discrete to continuous
                 
                    SSA_list[compartment_index] = max(SSA_list[compartment_index] - 1, 0)
                    PDE_list[self.PDE_multiple*compartment_index:self.PDE_multiple*(compartment_index+1)] += 1/self.h
                 
                
                t += tau 
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    SSA_grid[:, time_index] = SSA_list
                    PDE_grid[:, time_index] = PDE_list
                    approx_mass[:, time_index], PDE_particles[:,time_index] = self.calculate_total_mass(PDE_list, SSA_list)

                old_time = t  # Update old_time
                 # Update time by the time step
                #Printing a barrier

               
                print(f"{'Simulation Step':^30}")  # Centered title within the asterisks
                print("*" * 30)

                # Time and mass information
                print(f"Time: {t:.2f}")
                print(f"Mass conversion threshold: {self.threshold}")
                print("-" * 30)  # Separator line

                # Particle and mass details
                print(f"Stochastic particles in each box at time {t}:")
                print(f"  {SSA_list}")
                print(f"Continuous mass at time {t:.1f}:")
                print(f"  {PDE_list.round(1)}")
                print(f"Number of particles continuous")
                print(f" {PDE_particles[:,min(ind_after+1, len(self.time_vector))-1]}")
                print(f"Approximate mass at time {t:.1f}:")
                print(f"  {approx_mass[:, min(ind_after+1, len(self.time_vector))-1]}")
                print("-" * 30)

                # Propensity information
                print(f"{'Propensity Details':^30}")
                print(f"Index of reaction chosen: {index}")
                print("-" * 30)
                print(f"Movement propensity:           {total_propensity[:self.SSA_M]}")
                print(f"Production propensity:         {total_propensity[self.SSA_M:2*self.SSA_M]}")
                print(f"Degradation propensity:        {total_propensity[2*self.SSA_M:3*self.SSA_M]}")
                print(f"Conversion to discrete prop.:  {total_propensity[3*self.SSA_M:4*self.SSA_M]}")
                print(f"Conversion to continuous prop.: {total_propensity[4*self.SSA_M:]}")
                print("*" * 30)
                print("\n")  # Extra blank line for space between steps


            else:  # Else we run the ODE step

                #PDE_list = self.RK4(PDE_list)
                PDE_list = self.crank_nicholson(PDE_list)
                PDE_list = np.maximum(PDE_list, 0)  # Ensure non-negativity after RK4 step
                t = td
                td += self.timestep
                ind_before = np.searchsorted(self.time_vector, old_time, 'right')
                ind_after = np.searchsorted(self.time_vector, t, 'left')
                for time_index in range(ind_before, min(ind_after + 1, len(self.time_vector))):
                    PDE_grid[:, time_index] = PDE_list
                    SSA_grid[:, time_index] = SSA_list
                    approx_mass[:, time_index], PDE_particles[:,time_index]  = self.calculate_total_mass(PDE_list, SSA_list)
           
        return SSA_grid, PDE_grid, approx_mass


    def run_simulation(self,number_of_repeats:int) -> np.ndarray:
        """This will run the simulation with a total number of repeats"""
        PDE_initial, SSA_initial = self.create_initial_dataframe()
        approx_mass_initial = np.zeros_like(SSA_initial)
        approx_mass_initial[:,0] = self.calculate_total_mass(PDE_initial[:,0], SSA_initial[:,0])[0]
        SSA_sum = np.zeros_like(SSA_initial)
        PDE_sum = np.zeros_like(PDE_initial)
        approx_mass_sum = np.zeros_like(approx_mass_initial)

        for _ in tqdm(range(number_of_repeats),desc="Running the Hybrid simulations"):
            SSA_current, PDE_current, approx_mass_current = self.hybrid_simulation(deepcopy(SSA_initial), deepcopy(PDE_initial),deepcopy(approx_mass_initial))
            SSA_sum+= SSA_current
            PDE_sum += PDE_current
            approx_mass_sum += approx_mass_sum
        SSA_average = SSA_sum/number_of_repeats
        PDE_average = PDE_sum/number_of_repeats
        approx_sum_average = approx_mass_sum/number_of_repeats


        combined_grid = np.zeros_like(PDE_average)


        for i in range(SSA_average.shape[1]):
            for j in range(SSA_average.shape[0]):
                start_index = j*self.PDE_multiple
                end_index = (j+1)*self.PDE_multiple

                combined_grid[start_index:end_index,i] = PDE_average[start_index:end_index,i]+(1/self.h)*SSA_average[j,i]

        combined_grid[-1,:] = combined_grid[-2,:]
        

        # combined_grid = np.zeros_like(PDE_average)


        # for i in range(SSA_average.shape[1]):
        #     for j in range(SSA_average.shape[0]):
        #         start_index = j*self.PDE_multiple
        #         end_index = (j+1)*self.PDE_multiple

        #         combined_grid[start_index:end_index,i] = PDE_average[start_index:end_index,i]+(1/self.h)*SSA_average[j,i]

        # combined_grid[-1,:] = combined_grid[-2,:]
        

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
            'degredation_rate': self.degredation_rate,
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