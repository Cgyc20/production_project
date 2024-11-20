#include <stdio.h>

void ApproxMassLeftHand(int SSA_M, int PDE_M, float *PDE_list, float *approxMass, float deltax) {
    // works out the total mass of the PDE at each compartment using left hand rule
    float sum_value;
    int start_index, end_index;

    // Iterate over each compartment
    for (int i = 0; i < SSA_M; i++) {
        start_index = i * PDE_M;  // Start index for this compartment
        end_index = (i + 1) * PDE_M;  // End index for this compartment
        sum_value = 0.0;  // Initialize sum to 0.0 for each compartment
        
        // Sum the PDE values over the range corresponding to this compartment
        for (int j = start_index; j < end_index; j++) {  // Use '<' for correct range
            sum_value += PDE_list[j];  // Access values using pointer (PDE_list is a reference)
        }
        
        // Store the result in the approxMass array for this compartment
        approxMass[i] = sum_value * deltax;
    }
}

