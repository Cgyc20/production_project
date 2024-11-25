import ctypes


def set_clibrary_argtypes(clibrary):
    clibrary.ApproximateMassLeftHand.argtypes = [
        ctypes.c_int,  # SSA_M (int)
        ctypes.c_int,  # PDE_M (int)
        ctypes.POINTER(ctypes.c_float),  # PDE_list (pointer to float)
        ctypes.POINTER(ctypes.c_float),  # approxMass (pointer to float)
        ctypes.c_float  # deltax (float)
    ]

    clibrary.BooleanMass.argtypes = [
        ctypes.c_int,                  # SSA_m
        ctypes.c_int,                  # PDE_m
        ctypes.c_int,                  # PDE_multiple
        ctypes.POINTER(ctypes.c_float),  # PDE_list (float array)
        ctypes.POINTER(ctypes.c_int),  # boolean_PDE_list (int array)
        ctypes.POINTER(ctypes.c_int),  # boolean_SSA_list (int array)
        ctypes.c_float                 # h
    ]

    clibrary.CalculatePropensity.argtypes = [
        ctypes.c_int,  # SSA_M (int)
        ctypes.POINTER(ctypes.c_float),  # PDE_list (pointer to float)
        ctypes.POINTER(ctypes.c_int),  # SSA_list (pointer to int)
        ctypes.POINTER(ctypes.c_float),  # propensity_list (pointer to float)
        ctypes.POINTER(ctypes.c_float),  # combined_mass_list (pointer to float)
        ctypes.POINTER(ctypes.c_float),  # Approximate_PDE_Mass (pointer to float)
        ctypes.POINTER(ctypes.c_int),  # boolean_mass_list (pointer to int)
        ctypes.c_float,  # degradation_rate (float)
        ctypes.c_float,  # threshold 
        ctypes.c_float,  # Production_rate_PC (float)
        ctypes.c_float,  # gamma (float)
        ctypes.c_float  # jump_rate (float)
    ]