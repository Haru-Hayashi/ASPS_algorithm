#include "Inverter_function.h"

/**
 * @brief Initialize VCT_tVect3State
 * 
 * @param vs3 pointer of VCT_tVect3State
 */
void InitVect3(Vector_parameter* vs3){
    int i;
    for(i=0;i<2;++i){
        vs3->ab[i] = 0;
        vs3->dq[i] = 0;
        vs3->uvw[i] = 0;
    }
    vs3->uvw[2] = 0;
    vs3->phase = 0;
    vs3->nrm = 0;
}

/**
 * @brief Initialize the parameter of an inverter
 * 
 * @param ip initialized object
 * @param Vdc [V] DC link voltage
 */
bool Inv_UpdateOutputVoltage(Inverter_parameter* ip){
    int i_vv; int j_axis;
 
    // インバータの出力可能な電圧ベクトルを計算
    ip->Vdc = V_DC;  // 直流リンク電圧
    for(i_vv=0;i_vv<8;++i_vv){
        // 各出力電圧ベクトルの初期化
        InitVect3(&(ip->vout[i_vv]));
        // SWPから3相電圧出力を計算（中性点電位考慮なし）
        for(j_axis=0;j_axis<3;++j_axis){
            if(SwitchingState[i_vv][j_axis]==0){
                ip->vout[i_vv].uvw[j_axis] = -ip->Vdc/2;
            }else{
                ip->vout[i_vv].uvw[j_axis] = ip->Vdc/2;
            }
        }
        // 3相電圧出力から2相電圧出力を計算
        uvw2ab(ip->vout[i_vv].uvw[0], ip->vout[i_vv].uvw[1], ip->vout[i_vv].uvw[2], &(ip->vout[i_vv].ab[0]), &(ip->vout[i_vv].ab[1]));
        // 2相電圧出力から3相電圧出力を再計算
        ab2uvw(ip->vout[i_vv].ab[0], ip->vout[i_vv].ab[1], &(ip->vout[i_vv].uvw[0]), &(ip->vout[i_vv].uvw[1]), &(ip->vout[i_vv].uvw[2]));
    }
    return 0;
}

void uvw2ab(float u, float v, float w, float* a, float*b){
    *a = (float)(sqrt(2.0/3.0)*(u - 0.5*v - 0.5*w));
    *b = (float)(sqrt(2.0/3.0)*sqrt(3)/2.0*(v-w));
}

void ab2dq(float a, float b, float theta, float* d, float* q){
    float c = (float)cos(theta);
    float s = (float)sin(theta);
    *d = c*a + s*b;
    *q = -s*a + c*b;
}

void dq2ab(float d, float q, float theta, float* a, float*b){
    float c = (float)cos(theta);
    float s = (float)sin(theta);
    *a = c*d - s*q;
    *b = s*d + c*q;
}

void ab2uvw( float a, float b, float* u, float* v, float* w){
    *u = (float)(sqrt(2.0/3.0)*a);
    *v = (float)(sqrt(2.0/3.0)*(-0.5*a + sqrt(3)/2.0*b));
    *w = (float)(sqrt(2.0/3.0)*(-0.5*a - sqrt(3)/2.0*b));
}