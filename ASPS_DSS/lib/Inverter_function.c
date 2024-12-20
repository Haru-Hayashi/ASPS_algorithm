#include "Inverter_function.h"

/**
 * @brief 3相ベクトル変数の初期化
 * 
 * @param vs3 3相ベクトル変数のポインタ
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
 * @brief インバータの出力電圧を計算
 * 
 * @param ip オブジェクト
 * @param Vdc DCリンク電圧
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

/**
 * @brief 隣接ベクトルの判定
 * 
 * @param Step レイヤ
 * @param VV         電圧ベクトル番号
 * @param omega_ref  角周波数指令
 * @param Ts         演算周期
 */
bool AdjoinVector(double Step, int VV, double omega_ref, double Ts){
    double theta;
    
	theta = omega_ref*Ts*Step - pi*0.5;
	if(theta > pi){
		theta -= 2*pi;
	}

    if((theta >= -PI16) && (theta < PI16)){
		return(!((VV ==2) || (VV ==3) || (VV ==0))) ? 1 :0;
    }
    else if((theta >= PI16) && (theta < PI12)){
		return(!((VV ==3) || (VV ==4) || (VV ==0))) ? 1 :0;
    }
    else if((theta >= PI12) && (theta < PI56)){
		return(!((VV ==4) || (VV ==5) || (VV ==0))) ? 1 :0;
    }
    else if((pi >= theta && theta >= PI56) || (-pi <= theta && theta < -PI56)){
		return(!((VV ==5) || (VV ==6) || (VV ==0))) ? 1 :0;
    }
    else if((theta >= -PI56) && (theta < -PI12)){
		return(!((VV ==6) || (VV ==1) || (VV ==0))) ? 1 :0;
    }
    else if((theta >= -PI12) && (theta < -PI16)){
		return(!((VV ==1) || (VV ==2) || (VV ==0))) ? 1 :0;
    }
}

/**
 * @brief 電圧積分ベクトルの最終値を(0,0)に調整する
 * 
 * @param VI       オブジェクトｒ
 * @param Vector   オブジェクト
 * @param SwCnt    スイッチング回数
 */
bool Endpoint_adjust(Integral_model* VI, Inverter_parameter* Vector, uint32_t* SwCnt, uint32_t Layer){

    // VV=6を出力し続ける処理
    if(VI->theta_ref[Layer] > pi && VI->Pout.ab[1] < 1.5e-2){
        VI->Vout.ab[0] = Vector->vout[6].ab[0];
        VI->Vout.ab[1] = Vector->vout[6].ab[1];
        VI->VV_Num[0] = 6;
        if(VI->VV_Num[1] != VI->VV_Num[0]){
            *SwCnt++;
        }
    }

    // VV=1を出力し続ける処理
    if(VI->theta_ref[Layer] > pi && VI->VV_Num[0] == 6 && fabs(VI->Pout.ab[1]) < 1.0e-6){
        VI->Vout.ab[0] = Vector->vout[1].ab[0];
        VI->Vout.ab[1] = Vector->vout[1].ab[1];
        VI->VV_Num[0] = 1;
        if(VI->VV_Num[1] != VI->VV_Num[0]){
            *SwCnt++;
        }
    }

    // VV=0を出力し続ける処理
    if(VI->theta_ref[Layer] > pi && fabs(VI->Pout.ab[0]) < 1.0e-6){
        VI->Vout.ab[0] = Vector->vout[0].ab[0];
        VI->Vout.ab[1] = Vector->vout[0].ab[1];
        VI->VV_Num[0] = 0;
        if(VI->VV_Num[1] != VI->VV_Num[0]){
            *SwCnt++;
        }
    }

    return 0;
}

// **************** 座標変換の関数 **************** // 
void uvw2ab(double u, double v, double w, double* a, double*b){
    *a = (double)(sqrt(2.0/3.0)*(u - 0.5*v - 0.5*w));
    *b = (double)(sqrt(2.0/3.0)*sqrt(3)/2.0*(v-w));
}

void ab2dq(double a, double b, double theta, double* d, double* q){
    double c = (double)cos(theta);
    double s = (double)sin(theta);
    *d = c*a + s*b;
    *q = -s*a + c*b;
}

void dq2ab(double d, double q, double theta, double* a, double*b){
    double c = (double)cos(theta);
    double s = (double)sin(theta);
    *a = c*d - s*q;
    *b = s*d + c*q;
}

void ab2uvw( double a, double b, double* u, double* v, double* w){
    *u = (double)(sqrt(2.0/3.0)*a);
    *v = (double)(sqrt(2.0/3.0)*(-0.5*a + sqrt(3)/2.0*b));
    *w = (double)(sqrt(2.0/3.0)*(-0.5*a - sqrt(3)/2.0*b));
}