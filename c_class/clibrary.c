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
                         int *boolean_SSA_list, float *combined_mass_list, 
                         float *Approximate_PDE_Mass, int *boolean_mass_list, 
                         float degradation_rate_h, float threshold, 
                         float production_rate, float gamma, float jump_rate) {
    
    // Precompute commonly used values
    int two_SSA_M = 2 * SSA_M;
    int three_SSA_M = 3 * SSA_M;
    int four_SSA_M = 4 * SSA_M;
    int five_SSA_M = 5 * SSA_M;
    int six_SSA_M = 6 * SSA_M;

    // Initialize all propensity values to zero
    for (int i = 0; i < six_SSA_M; i++) {
        propensity_list[i] = 0.0f;
    }

    // 1. Diffusion Propensities (Boundary Conditions)
    propensity_list[0] = jump_rate * SSA_list[0];
    for (int i = 1; i < SSA_M - 1; i++) {
        propensity_list[i] = 2.0f * jump_rate * SSA_list[i];
    }
    propensity_list[SSA_M - 1] = jump_rate * SSA_list[SSA_M - 1];

    // 2. Production Rates (D → 2D)
    for (int i = SSA_M; i < two_SSA_M; i++) {
        int index = i - SSA_M;
        propensity_list[i] = production_rate * combined_mass_list[index] * boolean_SSA_list[index];
    }

    // 3. First Degradation: D + D → D
    for (int i = two_SSA_M; i < three_SSA_M; i++) {
        int index = i - two_SSA_M;
        propensity_list[i] = degradation_rate_h * SSA_list[index] * (SSA_list[index] - 1);
    }

    // 4. Second Degradation: D + C → C (Fixed missing factor of 2.0f)
    for (int i = three_SSA_M; i < four_SSA_M; i++) {
        int index = i - three_SSA_M;
        propensity_list[i] = 2.0f * degradation_rate_h * SSA_list[index] * Approximate_PDE_Mass[index];
    }

    // 5. Conversion from Continuous to Discrete (Below Threshold)
    for (int i = four_SSA_M; i < five_SSA_M; i++) {
        int index = i - four_SSA_M;
        float combined_mass = combined_mass_list[index];
        float approx_mass = Approximate_PDE_Mass[index];
        int boolean_mass = boolean_mass_list[index];

        propensity_list[i] = (combined_mass < threshold) 
                              ? gamma * approx_mass * boolean_mass 
                              : 0.0f;
    }

    // 6. Conversion from Discrete to Continuous (Above Threshold)
    for (int i = five_SSA_M; i < six_SSA_M; i++) {
        int index = i - five_SSA_M;
        float combined_mass = combined_mass_list[index];
        int SSA_mass = SSA_list[index];

        propensity_list[i] = (combined_mass >= threshold) 
                              ? gamma * SSA_mass 
                              : 0.0f;
    }
}