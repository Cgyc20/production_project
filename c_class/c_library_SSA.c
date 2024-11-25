#include <stdio.h>

void CalculatePropensity(int *SSA_list, int SSA_M, float *propensity_list, float jump_rate, float Production_rate_PC, float degradation_rate){

    int two_SSA_M = 2 * SSA_M;
    int three_SSA_M = 3 * SSA_M;

    // set the terms in the list to zero
    for (int i = 0; i < three_SSA_M; i++) {
        propensity_list[i] = 0.0f;
    }

    float jump_rate_f = jump_rate;

    propensity_list[0] = jump_rate_f * SSA_list[0];
    for (int i = 1; i < SSA_M - 1; i++) {
        propensity_list[i] = 2.0f * jump_rate_f * SSA_list[i];
    }
    propensity_list[SSA_M - 1] = jump_rate_f * SSA_list[SSA_M - 1];

    // Production rates (constant for each compartment)
    for (int i = SSA_M; i < two_SSA_M; i++) {
        propensity_list[i] = Production_rate_PC;
    }

    // Degradation rates (depends on SSA_list)
    float degradation_rate_f = degradation_rate;
    for (int i = two_SSA_M; i < three_SSA_M; i++) {
        propensity_list[i] = degradation_rate_f * SSA_list[i - two_SSA_M];
    }

}