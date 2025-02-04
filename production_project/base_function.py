import numpy as np
from tqdm import tqdm
from copy import deepcopy, copy
import ctypes

from production_project.clibrary_argtypes import set_clibrary_argtypes # Each data type for the c functions
clibrary = ctypes.CDLL("c_class/clibrary.so") # Import the C library
set_clibrary_argtypes(clibrary) # Import the data types for each C function

class UtilityFunctions:
    
    @staticmethod
    def calculate_total_mass(PDE_list: np.ndarray, SSA_list: np.ndarray, use_c_functions: bool, PDE_multiple: int, deltax: float, SSA_M: int) -> np.ndarray:
        """This will calculate the total mass of discrete + continuous"""
        PDE_list = PDE_list.astype(float)
        SSA_list = SSA_list.astype(int)
        if use_c_functions:
            mass_solver = UtilityFunctions.ApproximateLeftHandC
        else:
            mass_solver = UtilityFunctions.ApproximateLeftHandPython
        
        approximate_PDE_mass = mass_solver(PDE_list, PDE_multiple, deltax, SSA_M)
        combined_list = np.add(SSA_list, approximate_PDE_mass) 
        return combined_list, approximate_PDE_mass

    @staticmethod
    def ApproximateLeftHandPython(PDE_list: np.ndarray, PDE_multiple: int, deltax: float, SSA_M: int) -> np.ndarray:
        PDE_list = PDE_list.astype(float)
        approximation_number_cont = np.zeros(SSA_M)
        for i in range(SSA_M):
            start_index = PDE_multiple * i
            end_index = PDE_multiple * (i + 1)
            sum_value = np.sum(PDE_list[start_index:end_index]) * deltax
            approximation_number_cont[i] = sum_value
        return approximation_number_cont

    @staticmethod
    def ApproximateLeftHandC(PDE_list, PDE_multiple: int, deltax: float, SSA_M: int):
        approximate_PDE_mass = np.zeros(SSA_M)
        PDE_list = np.array(PDE_list, dtype=np.float32)
        PDE_list_Ctypes = PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

        approximate_PDE_mass_Ctypes = approximate_PDE_mass.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        clibrary.ApproximateMassLeftHand(SSA_M, PDE_multiple, PDE_list_Ctypes, approximate_PDE_mass_Ctypes, deltax)

        approximate_PDE_mass = np.ctypeslib.as_array(approximate_PDE_mass_Ctypes, shape=approximate_PDE_mass.shape)
        return approximate_PDE_mass

    @staticmethod
    def boolean_if_less_mass(PDE_list: np.ndarray, h: float, PDE_multiple: int, SSA_M: int) -> np.ndarray: 
        PDE_list = PDE_list.astype(float)
        boolean_PDE_list = np.zeros_like(PDE_list)
        boolean_PDE_list[PDE_list > 1 / h] = 1
        boolean_threshold_SSA = np.zeros(SSA_M)
        for i in range(SSA_M):
            start_index = i * PDE_multiple
            BOOL_VALUE = True
            for j in range(PDE_multiple):
                current_index = start_index + j
                if boolean_PDE_list[current_index] == 0:
                    BOOL_VALUE = False
                    break
            boolean_threshold_SSA[i] = 1 if BOOL_VALUE else 0
        return boolean_threshold_SSA

    @staticmethod
    def threshold_boolean(combined_list: np.ndarray, threshold: float, PDE_multiple: int, SSA_M: int) -> np.ndarray:
        """Generate a boolean list based on the threshold"""
        compartment_bool_list = np.array([0 if total_mass > threshold else 1 for total_mass in combined_list])
        PDE_bool_list = np.zeros(SSA_M * PDE_multiple)
        for i in range(SSA_M):
            value = compartment_bool_list[i]
            new_value = 0 if value == 1 else 1
            start_index = i * PDE_multiple 
            PDE_bool_list[start_index:start_index + PDE_multiple] = new_value
        return compartment_bool_list, PDE_bool_list


    @staticmethod
    def fine_grid_SSA_mass(SSA_mass: np.ndarray, PDE_X: np.ndarray, SSA_M: int, PDE_multiple: int, h: float):
        """Convert the SSA_mass to the same fine resolution as the PDE"""
        
        fine_SSA_mass = np.zeros_like(PDE_X, dtype=float)  # Ensure float type for PDE consistency
        
        for i in range(SSA_M):
            start_index = i * PDE_multiple
            end_index = (i + 1) * PDE_multiple
            fine_SSA_mass[start_index:end_index] = float(SSA_mass[i])/h

        return fine_SSA_mass
            
