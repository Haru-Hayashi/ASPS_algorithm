#ifndef ASPS_DSS_H_
#define ASPS_DSS_H_

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include "Inverter_function.h"

/**
 * @brief ASPSのためのDFS変数
 * 
 * @param Layer                 DFSのレイヤ(ステップ数)
 * @param SwCnt                 スイッチング回数
 * @param SwCnt_min             最小だったスイッチング回数
 * @param VV_Num[3]             電圧ベクトル番号 [0]:最新 [1]:前回 [2]:前々回 
 * @param VV_seq[100]           電圧ベクトルシーケンス (探索用)
 * @param VV_seq_record[100]  　電圧ベクトルシーケンス (最適シーケンスの保存用)
 */
typedef struct{
    uint32_t Layer;
    uint8_t SwCnt;
    uint8_t SwCnt_min;
    uint8_t VV_Num[3];
    uint8_t VV_seq[100];
    uint8_t VV_seq_record[100];
    Integral_model VI;
} ASPS_DFS;

void DFSbasedASPS();

#endif
