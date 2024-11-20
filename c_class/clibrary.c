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