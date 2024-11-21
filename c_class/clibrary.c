#include <stdio.h>
void ApproximateMassLeftHand(int SSA_M, int PDE_multiple, float *PDE_list, float *approxMass, float deltax) {
    // Works out the approximate mass of the PDE domain over each grid point using the left-hand rule
    int start_index, end_index;
    float sum_value;

    // Iterate over each compartment
    for (int i = 0; i < SSA_M; i++) {
        start_index = PDE_multiple * i;  // Start index for this compartment
        end_index = PDE_multiple * (i + 1);  // End index for this compartment
        sum_value = 0.0;  // Initialize sum for the current compartment

        // Sum the PDE values over the range corresponding to this compartment
        for (int j = start_index; j < end_index; j++) {  // Sum over the range
            sum_value += PDE_list[j];  // Add the value at each grid point
        }

        // Multiply by the step size (deltax) and store the result
        approxMass[i] = sum_value * deltax;  
    }
}

void BooleanMass(int SSA_m, int PDE_m, int PDE_multiple, float *PDE_list, int *boolean_PDE_list, int *boolean_SSA_list, float h) {
    // This will create boolean lists based on PDE and SSA conditions
    int start_index;
    int BOOL_value;
    int current_index;

    // Create the boolean_PDE_list based on the condition PDE_list[i] > 1/h
    for (int i = 0; i < PDE_m; i++) {  // Use `< PDE_m` to avoid out-of-bounds
        if (PDE_list[i] > 1 / h) {
            boolean_PDE_list[i] = 1;
        } else {
            boolean_PDE_list[i] = 0;
        }
    }

    // Create the boolean_SSA_list based on boolean_PDE_list
    for (int i = 0; i < SSA_m; i++) {  // Use `< SSA_m` to avoid out-of-bounds
        start_index = i * PDE_multiple;
        BOOL_value = 1; // Default to 1 (true)

        for (int j = 0; j < PDE_multiple; j++) {  // Iterate through the block of PDE_multiple
            current_index = start_index + j;

            // Ensure we do not access out-of-bounds in boolean_PDE_list
            if (current_index >= PDE_m || boolean_PDE_list[current_index] == 0) {
                BOOL_value = 0; // If any element is 0, the SSA block is 0
                break;          // Break early for efficiency
            }
        }

        boolean_SSA_list[i] = BOOL_value; // Assign the result to the SSA list
    }
}