import numpy as np
import ctypes


clibrary = ctypes.CDLL("c_class/clibrary.so") #import the c library

clibrary.ApproxMassLeftHand.argtypes = [
    ctypes.c_int,  # SSA_M (int)
    ctypes.c_int,  # PDE_M (int)
    ctypes.POINTER(ctypes.c_float),  # PDE_list (pointer to float)
    ctypes.POINTER(ctypes.c_float),  # approxMass (pointer to float)
    ctypes.c_float  # deltax (float)
]

# Assuming clibrary is already loaded and ApproxMassLeftHand is defined in it

def calculate_mass(PDE_list, SSA_list, SSA_M, PDE_M, deltax):
    approximate_PDE_mass = np.zeros_like(SSA_list)
    combined_list = np.zeros_like(SSA_list)

    # Convert numpy arrays to ctypes pointers
    PDE_list_Ctypes = PDE_list.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    approximate_PDE_mass_Ctypes = approximate_PDE_mass.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    # Call the C function
    clibrary.ApproxMassLeftHand(SSA_M, PDE_M, PDE_list_Ctypes, approximate_PDE_mass_Ctypes, deltax)

    # Convert the ctypes pointer back to numpy array
    approximate_PDE_mass = np.ctypeslib.as_array(approximate_PDE_mass_Ctypes, shape=approximate_PDE_mass.shape)

    # Combine the lists
    combined_list = np.add(SSA_list, approximate_PDE_mass)

    # Debug prints
    print("PDE_list:", PDE_list)
    print("approximate_PDE_mass:", approximate_PDE_mass)
    print("combined_list:", combined_list)

    return combined_list, approximate_PDE_mass

# Example usage
PDE_list = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
SSA_list = np.array([0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5], dtype=np.float32)
SSA_M = len(SSA_list)
PDE_M = 2
deltax = 0.1

combined_list, approximate_PDE_mass = calculate_mass(PDE_list, SSA_list, SSA_M, PDE_M, deltax)