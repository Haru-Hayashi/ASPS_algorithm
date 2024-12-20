#ifndef INVERTER_FUNCTION_H_
#define INVERTER_FUNCTION_H_

// *** header file include *** //
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include "ParameterDefined.h"
#include "Common.h"

// *** macro definition *** //

// *** type definition *** //

/**
 * @struct Vector_parameter
 * @brief 三相交流用の変数
 */
typedef struct{
    double uvw[3];
    double ab[2];
    double dq[2];
    double nrm;   
    double phase;
} Vector_parameter;

/**
 * @struct Inverter_parameter
 * @brief インバータの変数
 * 
 * @param Vdc [V]       DCリンク電圧
 * @param out_swp       スイッチング状態
 * @param vout[8] [V]   ベクトル番号ごとの出力電圧
 */
typedef struct{
    double Vdc;
    uint8_t out_swp;
    Vector_parameter vout[8];
} Inverter_parameter;

/**
 * @brief インバータ電圧積分モデル
 * 
 * @param VV_Num[3]   電圧ベクトル番号 [0]:最新 [1]:前回 [2]:前々回 
 * @param Pref_ab     電圧積分指令  
 * @param Pout        電圧積分出力   
 * @param Perr        電圧積分誤差
 * @param Perr_nrm    電圧積分誤差ノルム
 * @param Perr_sum    電圧積分誤差面積
 * @param theta_out   電圧積分出力位相
 * @param theta_ref   電圧積分指令位相
 */
typedef struct{
    uint8_t VV_Num[3];
    double Pref_a[1000];
    double Pref_b[1000];
    double Pout_nrm;
    double Pref_nrm;
    double Perr_nrm;
    double Perr_sum;
    double Perr_sum_min;
    double theta_out;
    double theta_ref[1000];
    Vector_parameter Vout;
    Vector_parameter Pout;
    Vector_parameter Perr;
} Integral_model;


static bool SwitchingState[8][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {1, 1, 0},
    {0, 1, 0},
    {0, 1, 1},
    {0, 0, 1},
    {1, 0, 1},
    {1, 1, 1}
};

static uint8_t SwCntTable[8][8] = {
	{0,1,2,1,2,1,2,3},
	{1,0,1,2,3,2,1,2},
	{2,1,0,1,2,3,2,1},
	{1,2,1,0,1,2,3,2},
	{2,3,2,1,0,1,2,1},
	{1,2,3,2,1,0,1,2},
	{2,1,2,3,2,1,0,1},
	{3,2,1,2,1,2,1,0}
};


// *** prototype definition *** //
// 座標変換
void uvw2ab(double u, double v, double w, double* a, double*b);
void ab2dq(double a, double b, double theta, double* d, double* q);
void dq2ab(double d, double q, double theta, double* a, double*b);
void ab2uvw( double a, double b, double* u, double* v, double* w);
//　関数
bool Inv_UpdateOutputVoltage(Inverter_parameter* ip);
bool AdjoinVector(double Step, int VV, double omega_ref, double Ts);
bool Endpoint_adjust(Integral_model* VI, Inverter_parameter* Vector, uint32_t* SwCnt, uint32_t Layer);

#endif