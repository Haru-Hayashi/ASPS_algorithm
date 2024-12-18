#ifndef INVERTER_FUNCTION_H_
#define INVERTER_FUNCTION_H_

// *** header file include *** //
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include "ParameterDefined.h"


// *** macro definition *** //

// *** type definition *** //

/**
 * @struct Vector_parameter
 * @brief 三相交流用の変数
 */
typedef struct{
    float uvw[3];
    float ab[2];
    float dq[2];
    float nrm;   
    float phase;
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
    float Vdc;
    uint8_t out_swp;
    Vector_parameter vout[8];
} Inverter_parameter;

/**
 * @brief インバータ電圧積分モデル
 * 
 * @param Pref_ab     電圧積分指令  
 * @param Pout        電圧積分出力   
 * @param Perr        電圧積分誤差
 * @param Perr_nrm    電圧積分誤差ノルム
 * @param Perr_sum    電圧積分誤差面積
 * @param theta_out   電圧積分出力位相
 * @param theta_ref   電圧積分指令位相
 */
typedef struct{
    float Pref_a[100000];
    float Pref_b[100000];
    float Perr_nrm;
    float Perr_sum;
    float theta_out;
    float theta_ref;
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


// *** prototype definition *** //
// 座標変換
void uvw2ab(float u, float v, float w, float* a, float*b);
void ab2dq(float a, float b, float theta, float* d, float* q);
void dq2ab(float d, float q, float theta, float* a, float*b);
void ab2uvw( float a, float b, float* u, float* v, float* w);
//　関数
bool Inv_UpdateOutputVoltage(Inverter_parameter* ip);

#endif