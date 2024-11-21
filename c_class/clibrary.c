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

// Function to calculate the propensities


void CalculatePropensity(int SSA_M, float *PDE_list, int *SSA_list, float *propensity_list, 
                         float *combined_mass_list, float *Approximate_PDE_Mass, 
                         int *boolean_mass_list, float degradation_rate, float threshold, 
                         float Production_rate_PC, float gamma, float jump_rate) {
    
    // Initialize all propensity values to zero
    for (int i = 0; i < 5 * SSA_M; i++) {
        propensity_list[i] = 0.0f;
    }

    // Diffusion propensities (boundary conditions)
    propensity_list[0] = jump_rate * (float)SSA_list[0];
    for (int i = 1; i < SSA_M - 1; i++) {
        propensity_list[i] = 2 * jump_rate * (float)SSA_list[i];
    }
    propensity_list[SSA_M - 1] = jump_rate * (float)SSA_list[SSA_M - 1];

    // Production rates (constant for each compartment)
    for (int i = SSA_M; i < 2 * SSA_M; i++) {
        propensity_list[i] = Production_rate_PC;
    }

    // Degradation rates (depends on SSA_list)
    for (int i = 2 * SSA_M; i < 3 * SSA_M; i++) {
        propensity_list[i] = degradation_rate * (float)SSA_list[i - 2 * SSA_M];
    }

    // Conversion from continuous to discrete (below threshold)
    for (int i = 3 * SSA_M; i < 4 * SSA_M; i++) {
        if (combined_mass_list[i - 3 * SSA_M] < threshold) {
            propensity_list[i] = gamma * Approximate_PDE_Mass[i - 3 * SSA_M] * (float)boolean_mass_list[i - 3 * SSA_M];
        }
    }

    // Conversion from discrete to continuous (above threshold)
    for (int i = 4 * SSA_M; i < 5 * SSA_M; i++) {
        if (combined_mass_list[i - 4 * SSA_M] >= threshold) {
            propensity_list[i] = gamma * (float)SSA_list[i - 4 * SSA_M];
        }
    }
}