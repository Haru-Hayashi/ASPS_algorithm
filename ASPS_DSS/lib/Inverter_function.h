#ifndef INVERTER_FUNCTION_H_
#define INVERTER_FUNCTION_H_

// *** header file include *** //
#include <math.h>
#include <stdbool.h>

// *** macro definition *** //

// *** type definition *** //
/**
 * @brief インバータ電圧積分モデル
 * 
 * @param Pout
 * 
 */

typedef struct{
    float Pout_a[3], Pout_b[3]; 
    float Pref_a[100000];
    float Pref_b[100000];
    float Perr[3];
    float Perr_sum;
    float theta_out;
    float theta_ref;
} Integral_model;

#endif