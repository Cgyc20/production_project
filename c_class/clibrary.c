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

void CalculatePropensitySub(int SSA_M, float *PDE_list, int *SSA_list, float *propensity_list, float degradation_rate, float Production_rate_PC, float jump_rate, int *boolean_mass_list){
    // works out the propensity functions
    // first we will in the diffusion rates
    float boolean_mass_list;

    propensity_list[0] = 2*jump_rate*SSA_list[0];
    propensity_list[SSA_M-1] = 2*jump_rate*SSA_list[SSA_M-1];
    for (int i=1; i<=SSA_M-2; i++){
        propensity_list[i] = SSA_list[i]*jump_rate;
    }

    // now we fill in the production_rates

    for (int i =SSA_M; i<=2*SSA_M-1,i++){
        propensity_list[i] = Production_rate_PC;
    }
    
    // now the degradation
    for (int i = 2*SSA_M; i<=3*SSA_M-1,i++){
        propensity_list[i] = degradation_rate*SSA_list[i];
    }

    boolean_mas_list

}