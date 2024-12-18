#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#include "ASPS.h"
#include "lib/Common.h"
#include "lib/ParameterDefined.h"
#include "lib/Inverter_function.h"

// オブジェクト
Integral_model Integral_param;
int f_ref;
double omega_ref;
double Ts = 10e-6;

int i;
int j;
int k;
double t;
double A;
double theta;

// VIモデル関連
float Pref_a_tmp, Pref_b_tmp;

static double Vout_u[3], Vout_v[3], Vout_w[3];
double Vout_a[3], Vout_b[3];
double Vref_a[3], Vref_b[3];
static double Pout_a[3], Pout_b[3];
static double Perr_a[3], Perr_b[3];
double Vref_nrm = 0.0;
double Perr_nrm[3];
double once_Perr;
double Vest_a[8], Vest_b[8];
static double Pest_a[8], Pest_b[8];
double Vtheta_ref_tmp;
double Vtheta2;
double min_Perr = 1e6;
int min_vector;

int u = 0, v = 0, w = 0;
static int U[2];
double Vu_inv = 0.0, Vv_inv = 0.0, Vw_inv = 0.0;
double Vwu = 0.0, Vuv = 0.0, Vvw = 0.0;
double Vu = 0.0, Vv = 0.0, Vw = 0.0;

static double Vz;
double V_amp;
double mod;
int k = 0;
int n = 0;
// int N[2] = {0,1};
int N_act;
int N = 0;
long long Pattern_Count;
int P = 0, P2 = 0, P3 = 0;
int tree[100][10000];
int SwCnt[100];
int vector[3];
int Vector_tmp;
int vector0_tmp[2];
int StepLevel[2] = {0,0};
int StepVIRTM;
//volatile int Line_flag[10] = 0x00;
int Line_flag[1000][8];
volatile double Boundary = 0.0;
double nothing;
int Flag;
double SwFreq;
int temp;
double t;
int SwCnt_min;
int SwCnt_min_total;
int optimal_pattern;
int optimal_count;
int jump;
int sector;
double phase;
int renew;
int renew_max;
int Sw_rep[100000];
int C_eva;
int C_eva_max;
int lastplot;
double Mod;

double Err_a, Err_b;
double Grad_a, Grad_b;
double N_lim_a, N_lim_b;

int DFS_Loop;
int DFS_MAX;
int Step_MAX;

double Perr_nrm_eva[2];
double Perr_nrm_eva_min;
double Perr_nrm_eva_limit;
double Perr_nrm_MAX;
double Perr_nrm_tmp[100];
double Pref_nrm;
double Pout_nrm;
double Perr_nrm_total;
int flag_zero;
int flag_Pa_zero;
int threshold;
long long search_num_MAX;
double phase_int;
double StepLevel_cont;
int DFS_MAX_SUM;
int SwCnt_limit;
double Kvf;
int Flag_DFS_MAX;

double W;

int i_temp[1000];

clock_t start_clock, end_clock;
clock_t start_clock2, end_clock2;

/*確認用変数*/////////
int Vector[3];
int Vector_enc[2];
int Count_enc;
int flag_s;
int Loop_Time;
long long Node_Counter;
int flag_OT;
int Perr_skip;
int Perr_sum;
///////////////////

int SwCntTable[8][8] = 
{
	{0,1,2,1,2,1,2,3},
	{1,0,1,2,3,2,1,2},
	{2,1,0,1,2,3,2,1},
	{1,2,1,0,1,2,3,2},
	{2,3,2,1,0,1,2,1},
	{1,2,3,2,1,0,1,2},
	{2,1,2,3,2,1,0,1},
	{3,2,1,2,1,2,1,0}
};



int main(void){
	// ファイルの保存先フォルダとファイル名を指定
    const char *folderPath = "C:\\Haru_Hayashi\\ASPS_algorithm\\data\\"; // フォルダのパス
	const char *fileNames[] = {"test.csv", "test2.csv", "test3.csv"};	

	// FILE * ポインタの配列を宣言
    FILE *fps[FileNum]; 

	for (int i = 0; i < FileNum; i++) {
		// ファイルパスを結合
		char filePath[255]; // ファイルパスのバッファ
		snprintf(filePath, sizeof(filePath), "%s%s", folderPath, fileNames[i]);
		// ファイルを書き込みモード("w")で開く
        fps[i] = fopen(filePath, "w");
	}

	// ***** 変調率と角周波数の設定 ***** //
	Mod = Mod_Max;
	omega_ref = Omega_rated*Mod/Mod_Max;
	f_ref =omega_ref/(2*pi);

	// ******* 演算周期とDFSの設定 ******* //
	Ts = 25e-6;  // 演算周期
	Step_MAX = floor((double)1/f_ref/Ts);    // 1周期のステップ数
	DFS_MAX = floor(pi/(6*omega_ref*Ts));  // 1周期のステップ数を12分割 
	DFS_Loop = floor(Step_MAX/DFS_MAX) + 1;  // DFSを実行する回数
	DFS_Loop = 2; 

	// ********** 許容値の設定 ********** //
	Boundary = 12.0e-3;       // 境界幅 (境界付きで探索するときに使用)
	SwCnt_limit = 6;          // スイッチング回数の上限値
	SwCnt_min = SwCnt_limit;  // 現在のスイッチング回数の最小値
	Perr_nrm_eva_limit = 0.5; // DSSの上限値

	Kvf = V_DC*1/sqrt(2)/360;
	Vref_nrm = Kvf*omega_ref;

// **************** ここからアルゴリズムの記述 **************** //
	start_clock2 = clock();

	// ********** DFSの基づくASPSの処理開始(基本12回ループ) ********** //
	for(k = 0; k < DFS_Loop; ++k){
		Loop_Time++;
		start_clock = clock();

		// VI指令値の計算とデータ格納 
		for(i = 0; i < DFS_MAX * Mul_Reso; ++i){
			Pref_a_tmp += Vref_nrm*cos(Integral_param.theta_ref)*Ts/Mul_Reso;
			Pref_b_tmp += Vref_nrm*cos(Integral_param.theta_ref-pi*(0.5))*Ts/Mul_Reso;
			Integral_param.Pref_a[i] = Pref_a_tmp;
			Integral_param.Pref_b[i] = Pref_b_tmp;
			Integral_param.theta_ref += omega_ref*Ts/Mul_Reso;
		}

		loop:
		Node_Counter = 0;
		VectorSearch();
		
		DFS_MAX_SUM += DFS_MAX;
		N = N_act;

		// VV_sequenceからVIを計算してデータを保存
		for(i = 1; i <= DFS_MAX; ++i){

			VectorVoltage(tree[i][N], V_DC, &Vest_a[0], &Vest_b[0]);
			Vector[0] = tree[i][N];

			/******************最終値と初期値をつなげる処理******************/
			//ベクトル6を打ち続ける処理
			if(Integral_param.theta_ref > 5.8 && Pout_b[0] < 1.5e-2){
				VectorVoltage(6, V_DC, &Vest_a[0], &Vest_b[0]);
				Vector[0] = 6;
				if(Vector[1] != Vector[0] && flag_zero != 1 && flag_zero != 2){
					SwCnt_min_total++;
				}
			}
			//出力値が初期値に戻ったら零ベクトルを打ち続けるフラグ
			if(Integral_param.theta_ref > 5.8 && fabs(Pout_a[0]) < 1.0e-6){
				//Count_enc = 1;
				flag_zero = 2;
			}
			//ベクトル1を打ち続ける処理
			if(flag_zero == 1){
				VectorVoltage(1, V_DC, &Vest_a[0], &Vest_b[0]);
				Vector[0] = 1;
				if(Vector[1] != Vector[0]){
					SwCnt_min_total++;
				}
			//}
			//零ベクトルを打ち続ける
			}else if(flag_zero == 2){
				VectorVoltage(0, V_DC, &Vest_a[0], &Vest_b[0]);
				Vector[0] = 0;
				if(Vector[1] != Vector[0]){
					SwCnt_min_total++;
				}
			}
			/***************************************************************/

			Pout_a[0] += Vest_a[0]*Ts;
			Pout_b[0] += Vest_b[0]*Ts;

			/******************最終地と初期値をつなげる処理******************/
			//ベクトル6を打ち続けてα軸にあたるときのフラグ
			if(Integral_param.theta_ref > 5.8 && Vector[0] == 6 && fabs(Pout_b[0]) < 1.0e-6){
				flag_zero = 1;
			}
			/**************************************************************/

			Perr_a[0] = Pout_a[0] - Integral_param.Pref_a[i*Mul_Reso-Mul_Reso];
			Perr_b[0] = Pout_b[0] - Integral_param.Pref_b[i*Mul_Reso-Mul_Reso];
			Perr_nrm[0] = sqrt(Perr_a[0]*Perr_a[0]+Perr_b[0]*Perr_b[0]);
			Perr_nrm_total += Perr_nrm[0];
			Pout_nrm = sqrt(Pout_a[0]*Pout_a[0] + (Pout_b[0]-Kvf)*(Pout_b[0]-Kvf));
			Pref_nrm = sqrt(Integral_param.Pref_a[i*Mul_Reso-Mul_Reso]*Integral_param.Pref_a[i*Mul_Reso-Mul_Reso]+(Integral_param.Pref_b[i*Mul_Reso-Mul_Reso]-Kvf)*(Integral_param.Pref_b[i*Mul_Reso-Mul_Reso]-Kvf));

			//ベクトル番号をエンコードして保存
			if(Vector[0] != Vector[1]){
				fprintf(fps[1], "%d, %d\n", Vector[1], Count_enc);
				Count_enc = 1;
			}else {
				Count_enc++;
			}
			Vector[1] = Vector[0];

			// t = Ts*(i-1)+Ts*DFS_MAX*k;
			t += Ts;

			//指令値が初期値にきたら終わる
			if(Integral_param.theta_ref > 5.8 && Integral_param.Pref_a[i*Mul_Reso-Mul_Reso] > Vref_nrm*Ts){
				Count_enc -= 1;
				fprintf(fps[1], "%d, %d\n", Vector[1], Count_enc);
				Count_enc = 1;
				goto stop;
			}
			// fprintf(fps[0], "%f, %f, %f, %f, %f, %f, %f, %f, %d\n", t, Pout_a[0], Pout_b[0], Integral_param.Pref_a[i-1], Integral_param.Pref_b[i-1], Boundary, Perr_nrm[0], Vtheta_tmp[i-1], Vector[0]);
			fprintf(fps[0], "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d\n", t, Pout_a[0], Pout_b[0], Integral_param.Pref_a[i*Mul_Reso-Mul_Reso], Integral_param.Pref_b[i*Mul_Reso-Mul_Reso], Boundary, Perr_nrm[0], Pout_nrm, Kvf+Boundary, Kvf-Boundary, Pref_nrm, Integral_param.theta_ref, Vector[0]);
		}

		Pest_a[1] = Pout_a[0];
		Pest_b[1] = Pout_b[0];

		vector[0] = tree[DFS_MAX][N];
		Vector_tmp = tree[DFS_MAX-1][N];
		vector0_tmp[0] = tree[DFS_MAX][N];
		vector0_tmp[1] = tree[DFS_MAX-1][N];
		 

		for(j = 0; j < 2; ++j){
			SwCnt[j] = 0;
			Perr_nrm_eva[j] = 0;
		}
		for(i = 0; i <= DFS_MAX; ++i){
			for(j = 0; j < 2; ++j){
				tree[i][j] = 0;
			}
		}

		end_clock = clock();
		// printf("Loop_TIme=%d, Switching_min=%d\n",Loop_Time, SwCnt_min);
		printf("Loop_TIme/Loop_MAX=%d/%d, Clock_Time=%f\n",Loop_Time, DFS_Loop, (double)(end_clock-start_clock)/CLOCKS_PER_SEC);
		printf("Search_Number=%d, DFS_MAX=%d\n",Node_Counter, DFS_MAX);
		// fprintf(fps[2], "%d, %d, %d, %f\n", Loop_Time, Node_Counter, DFS_MAX, (double)(end_clock-start_clock)/CLOCKS_PER_SEC);
		Node_Counter = 0;

		SwCnt_min_total += SwCnt_min;
		SwCnt_min = SwCnt_limit;

		StepLevel[0] = 0;

		tree[0][0] = vector[0]; 
		tree[0][1] = vector[0]; 

	}
	stop:

	end_clock2 = clock();
///////////////////////////////////////////////////////////
		
	SwFreq = (SwCnt_min_total/3)/(DFS_MAX*DFS_Loop*Ts);

	printf("DFS_MAX=%d\n", DFS_MAX);
	printf("DFS_Loop=%d\n", DFS_Loop);
	printf("Step_Max=%d\n", Step_MAX);
	printf("Ts[us]=%f\n", Ts*1e6);
	printf("Boundary=%f\n", Boundary);
    printf("f_ref=%d\n", f_ref);
	printf("flag_zero=%d\n", flag_Pa_zero);
	printf("Minimum SwCnt=%d\n", SwCnt_min_total);
	printf("Minimum SwFreq=%f\n", SwFreq);
	printf("Optimal Count=%d\n", Sw_rep[renew]+1);
	printf("SwCnt_renew=%d\n", renew);
	printf("How many patterns=%d\n", Pattern_Count);
	printf("last_plot=%d\n", lastplot);
	printf("CPU_clock=%f\n", (double)(end_clock2-start_clock2)/CLOCKS_PER_SEC);
	printf("Perr_skip=%d\n", Perr_skip);
	printf("Perr_nrm_limit=%f\n", Perr_nrm_eva_limit);
	printf("Perr_nrm_total=%f\n", Perr_nrm_total);
	printf("Perr_sum=%d\n", Perr_sum);

    for(i = 0; i < FileNum; i++){
		// ファイルを閉じる
        fclose(fps[i]);
	}
	
}


void VectorSearch(){

	for(i = 0; i < VV_NUM; ++i){

		VectorVoltage(i, V_DC, &Vest_a[i], &Vest_b[i]);
		Pest_a[1] += Ts*Vest_a[i];
		Pest_b[1] += Ts*Vest_b[i];
		Perr_a[1] = Pest_a[1] - Integral_param.Pref_a[StepLevel[0]*Mul_Reso];
		Perr_b[1] = Pest_b[1] - Integral_param.Pref_b[StepLevel[0]*Mul_Reso];
		Perr_nrm[1] = sqrt(Perr_a[1]*Perr_a[1]+Perr_b[1]*Perr_b[1]);

		//探索ベクトル制限(非零ベクトル2,零ベクトル2,合計4本)
		StepLevel_cont = StepLevel[0]+DFS_MAX_SUM;
		if(PhaseCalc(&sector, StepLevel_cont, i, omega_ref, Ts, &phase_int) == 1){
			goto END;
		}
	
		//許容値に入っているか比較する。入ってなければfor文がまわる
		// if(fabs(Perr_nrm[1]) < fabs(Boundary)){ 
			i_temp[StepLevel[0]] = i;
			vector[2] = vector[1];
			vector[1] = vector[0];  //前回のベクトル番号を保存しておく
			StepLevel[0]++;  //ステップ数を更新
			vector[0] = i;  //格納するベクトル番号
			if(vector[0] == 0 && SwCntTable[7][vector[1]] < SwCntTable[0][vector[1]]){
				vector[0] = 7;
			}
			tree[StepLevel[0]][N] = vector[0];  //ベクトル番号格納
			SwCnt[N] += SwCntTable[vector[0]][vector[1]];  //ベクトルパターンの合計スイッチング回数を更新
	
			//スイッチング回数による枝切り
			if(SwCnt[N] > SwCnt_min){
				SwCnt[N] -= SwCntTable[vector[0]][vector[1]];
				StepLevel[0]--;
				vector[0] = vector[1];
				vector[1] = vector[2];
				goto END;
			}

			Perr_nrm_eva[N] += Perr_nrm[1]; //スイッチングの枝切りがなければ評価エラー加算
			Perr_nrm_tmp[StepLevel[0]] = Perr_nrm[1]; 

			//電圧積分誤差による枝切り
			if(Perr_nrm_eva[N] > Perr_nrm_eva_limit){
				Perr_skip++; //確認用
				Perr_nrm_eva[N] -= Perr_nrm[1];
				SwCnt[N] -= SwCntTable[vector[0]][vector[1]];
				StepLevel[0]--;
				vector[0] = vector[1];
				vector[1] = vector[2];
				goto END;
			}

			// Perr_nrm_tmp[StepLevel[0]] = Perr_nrm[1]; 

			/****************************任意のステップ数になった時の処理****************************/
			if(StepLevel[0] == DFS_MAX){
				Flag_DFS_MAX = 1;
				//スイッチング回数の最小値が更新されたらベクトル保存
				if(SwCnt[N] < SwCnt_min){
					SwCnt_min = SwCnt[N]; //スイッチング回数の最小値を保存
					Perr_nrm_eva_min = Perr_nrm_eva[N]; //トータルエラーの最小値として保存しておく(これを基準にする)
					N_act = N;
					N = !N;
					SwCnt[N] = SwCnt[!N];
					Perr_nrm_eva[N] = Perr_nrm_eva[!N];
					for(j = 0; j <= DFS_MAX; ++j){
						tree[j][N] = tree[j][!N];
					}
					renew++;                 
				}

				//スイッチング回数が最小値と同じだったら,トータルエラーが最小のものを選ぶ
				/***************************************************************************/
				if(SwCnt[N] == SwCnt_min && Perr_nrm_eva[N] < Perr_nrm_eva_min){
					N_act = N;
					N = !N;
					SwCnt[N] = SwCnt[!N];
					Perr_nrm_eva[N] = Perr_nrm_eva[!N];
					for(j = 0; j <= DFS_MAX; ++j){
						tree[j][N] = tree[j][!N];
					}
					Perr_nrm_eva_min = Perr_nrm_eva[N]; //トータルエラーの最小値を保存
				}
				/***************************************************************************/
				
				SwCnt[N] -= SwCntTable[vector[0]][vector[1]]; //ひとつ前の階層に戻るのでスイッチングカウントをひとつ前に戻す
				Perr_nrm_eva[N] -= Perr_nrm[1]; //エラーもひとつ前に戻す
				Pattern_Count++; //ベクトルパターン数更新(確認用)
				vector[0] = vector[1];
				vector[1] = tree[StepLevel[0]-2][N];
				StepLevel[0]--;
				goto END;  //再起呼び出しをスキップ					
			}
			/****************************任意のステップ数になった時の処理終わり****************************/

		/***********再帰呼び出し***********/
		VectorSearch();
		/*********************************/

		i = i_temp[StepLevel[0]]; //関数から戻った時にその階層のfor文の場所を代入

		/******探索数制限を超えていたらreturn******/
		if(flag_OT == 1){
			return;
		}
		/************************************/

		END: //枝切時とゴール後の再起呼び出しスキップ場所

		VectorVoltage(i, V_DC, &Vest_a[i], &Vest_b[i]);
		Pest_a[1] -= Ts*Vest_a[i];  //前回の積分値に戻る
		Pest_b[1] -= Ts*Vest_b[i];	//同上
	}	
	/* -----for文がすべて回った場合の処理---- */
	SwCnt[N] -= SwCntTable[vector[1]][vector[0]];
	Perr_nrm_eva[N] -= Perr_nrm_tmp[StepLevel[0]]; //進んだ先に見つからなかったので戻る
	vector[0] = vector[1];
	if(StepLevel[0]-2 < 0){
		vector[1] = Vector_tmp;
	}else{
		vector[1] = tree[StepLevel[0]-2][N];
	}
	StepLevel[0]--;  //ステップ数更新

	return;	
}


void VectorVoltage(int i, double V_dc, double *Vest_a, double *Vest_b){
	if(i == 0){u = 0; v = 0; w = 0;}
	if(i == 1){u = 1; v = 0; w = 0;}
	if(i == 2){u = 1; v = 1; w = 0;}
	if(i == 3){u = 0; v = 1; w = 0;}
	if(i == 4){u = 0; v = 1; w = 1;}
	if(i == 5){u = 0; v = 0; w = 1;}
	if(i == 6){u = 1; v = 0; w = 1;}
	if(i == 7){u = 1; v = 1; w = 1;}
	//電圧ベクトル計算
	Vu_inv = V_DC*1*u;
	Vv_inv = V_DC*1*v;
	Vw_inv = V_DC*1*w;
	Vwu = Vw_inv-Vu_inv;
	Vuv = Vu_inv-Vv_inv;
	Vvw = Vv_inv-Vw_inv;
	Vu = -(Vwu-Vuv)/3;
	Vv = -(Vuv-Vvw)/3;
	Vw = -(Vvw-Vwu)/3;
	*Vest_a = (double)(sqrt(2.0/3.0)*(Vu - 0.5*Vv - 0.5*Vw));
	*Vest_b = (double)(sqrt(2.0/3.0)*sqrt(3)/2.0*(Vv-Vw));

	return;
}

void OutputVoltage(int out_swp, double V_dc, double *Vout_u, double *Vout_v, double *Vout_w){
	if(out_swp == 0){u = 0, v = 0, w = 0;}
	if(out_swp == 1){u = 1, v = 0, w = 0;}
	if(out_swp == 2){u = 1, v = 1, w = 0;}
	if(out_swp == 3){u = 0, v = 1, w = 0;}
	if(out_swp == 4){u = 0, v = 1, w = 1;}
	if(out_swp == 5){u = 0, v = 0, w = 1;}
	if(out_swp == 6){u = 1, v = 0, w = 1;}
	if(out_swp == 7){u = 1, v = 1, w = 1;}
	//電圧ベクトル計算
	Vu_inv = V_DC/1*u;
	Vv_inv = V_DC/1*v;
	Vw_inv = V_DC/1*w;
	Vwu = Vw_inv-Vu_inv;
	Vuv = Vu_inv-Vv_inv;
	Vvw = Vv_inv-Vw_inv;
	*Vout_u = -(Vwu-Vuv)/3;
	*Vout_v = -(Vuv-Vvw)/3;
	*Vout_w = -(Vvw-Vwu)/3;

	return;
}

int PhaseCalc(int *sector, double Step_Level, int x, double omega_ref, double Ts, double *phase_int){
	// double *phase_int;
	// *phase = atan2(B,A);
	// *phase_int = *phase - pi/2;
	// *phase_int = theta - pi/2;
	// *phase_int = *phase;

	*phase_int = omega_ref*Ts*Step_Level - pi*0.5;

	if(*phase_int > pi){
		*phase_int -= 2*pi;
	}

    if((*phase_int >= -pi/6) && (*phase_int < pi/6))
    {
        *sector = 1; 
		if(!((x ==2) || (x ==3) || (x ==0))){
			// *jump = 1;
			return 1;
		}   
			return 0;
    }
    else if((*phase_int >= pi/6) && (*phase_int < pi/2))
    {
        *sector = 2;
		if(!((x ==3) || (x ==4) || (x ==0))){
			// *jump = 1;
			return 1;
		} 
			return 0;
    }
    else if((*phase_int >= pi/2) && (*phase_int < pi*5/6))
    {
        *sector = 3;
		if(!((x ==4) || (x ==5) || (x ==0))){
			// *jump = 1;
			return 1;
		} 
			return 0;
    }
    else if((pi >= *phase_int && *phase_int >= pi*5/6) || (-pi <= *phase_int && *phase_int < -pi*5/6))
    {
        *sector = 4;
		if(!((x ==5) || (x ==6) || (x ==0))){
			// *jump = 1;
			return 1;
		} 
			return 0;
    }
    else if((*phase_int >= -pi*5/6) && (*phase_int < -pi/2))
    {
        *sector = 5;
		if(!((x ==6) || (x ==1) || (x ==0))){
			// *jump = 1;
			return 1;
		} 
			return 0;
    }
    else if((*phase_int >= -pi/2) && (*phase_int < -pi/6))
    {
        *sector = 6;
		if(!((x ==1) || (x ==2) || (x ==0))){
			// *jump = 1;
			return 1;
		} 
			return 0;
    }
}

