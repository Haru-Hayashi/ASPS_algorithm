#ifndef ASPS_H_
#define ASPS_H_

#include <math.h>


void VectorSearch();
void OutputVoltage(int out_swp, double Vdc, double *Vout_u, double *Vout_v, double *Vout_w);
void VectorVoltage(int i, double Vdc, double *Vest_a, double *Vest_b);
int PhaseCalc(int *sector, double Step_Level, int x, double omega_ref, double Ts, double *phase_int);
int PhaseCalc2(double A, double B, double theta, int x, int *sector, int  *jump, double *phase);

#endif
