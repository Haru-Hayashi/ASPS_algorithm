#ifndef USERDEFINED_H_
#define USERDEFINED_H_

#include <math.h>

#define FileNum 3   //　CSVデータで保存するファイル数

#define Mod_Max 1.15470053837925  // 最大変調率(六角形の内接円)
#define Omega_rated 360  // 定格角周波数

#define Mul_Reso 100  // 指令値の演算周期 -> 演算周期 Ts * Mul_Reso
#define VV_NUM 7  // 探索するベクトルの本数 -> VV_NUM-1
#define V_DC 282  // DCリンク電圧

#endif