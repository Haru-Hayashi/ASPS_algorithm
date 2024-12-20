#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#include "ASPS_DSS.h"
#include "lib/Common.h"
#include "lib/ParameterDefined.h"
#include "lib/Inverter_function.h"

// オブジェクト
Integral_model VI_data;
Inverter_parameter Inv_param;
ASPS_DFS dfs;

int f_ref;
double omega_ref;
double Ts = 10e-6;

int i, j, k;
double time_data;

// インバータモデル
float Pref_a_tmp, Pref_b_tmp;
double Vref_nrm = 0.0;
double Mod;

// DFS
int DFS_Loop;
int DFS_MAX;
int Step_MAX;
int Vector_tmp;
int vector0_tmp[2];
volatile double Boundary = 0.0;
double SwFreq;
double time_data;
int SwCnt_total;
long long Pattern_Update;

double Perr_sum_limit;
double Perr_nrm_tmp[100];
double Pref_nrm;
double Pout_nrm;
double Perr_nrm_total;
int flag_zero;
int flag_Pa_zero;
int threshold;
long long search_num_MAX;
double phase_int;
double Total_Step;
int DFS_MAX_SUM;
int SwCnt_limit;
double Kvf;
double W;

int VV_node_tmp[1000];

double Timer_clock[2];
clock_t start_clock, end_clock;
clock_t start_clock2, end_clock2;

/*確認用変数*/////////
int Count_enc;
int Loop_Time;
long long Node_Counter;
///////////////////

int main(void){
	// ファイルの保存先フォルダとファイル名を指定
	const char *folderPath = "D:\\Remoto\\ASPS_algorithm\\ASPS_DSS\\data\\";  // フォルダのパス
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
	SwCnt_limit = 5;          // スイッチング回数の上限値
	dfs.SwCnt_min = SwCnt_limit;  // 現在のスイッチング回数の最小値
	Perr_sum_limit = 0.7; // DSSの上限値

	Kvf = V_DC*1/sqrt(2)/360;
	Vref_nrm = Kvf*omega_ref;

// **************** ここからアルゴリズムの記述 **************** //
	start_clock2 = clock();

	// ********** DFSの基づくASPSの処理開始(基本12回ループ) ********** //

	//　インバータの出力可能な電圧ベクトルを計算
	Inv_UpdateOutputVoltage(&Inv_param);

	for(k = 0; k < DFS_Loop; ++k){
		Loop_Time++;
		start_clock = clock();

		// VI指令値の計算とデータ格納 
		for(j = 0; j < DFS_MAX * Mul_Reso; ++j){
			Pref_a_tmp += Vref_nrm*cos(VI_data.theta_ref)*Ts/Mul_Reso;
			Pref_b_tmp += Vref_nrm*cos(VI_data.theta_ref-pi*(0.5))*Ts/Mul_Reso;
			VI_data.Pref_a[j] = Pref_a_tmp;
			VI_data.Pref_b[j] = Pref_b_tmp;
			VI_data.theta_ref += omega_ref*Ts/Mul_Reso;
		}

		loop:
		Node_Counter = 0;
		DFSbasedASPS();
		
		DFS_MAX_SUM += DFS_MAX;

		// VV_sequenceからVIを計算してデータファイルに保存
		for(j = 1; j <= DFS_MAX; ++j){
			VI_data.Vout.ab[0] = Inv_param.vout[dfs.VV_seq_record[j]].ab[0];
			VI_data.Vout.ab[1] = Inv_param.vout[dfs.VV_seq_record[j]].ab[1];
			VI_data.VV_Num[0] = dfs.VV_seq_record[j];

			/******************最終値と初期値をつなげる処理******************/
			// ベクトル6を打ち続ける処理
			if(VI_data.theta_ref > 5.8 && VI_data.Pout.ab[1] < 1.5e-2){
				VI_data.Vout.ab[0] = Inv_param.vout[6].ab[0];
				VI_data.Vout.ab[1] = Inv_param.vout[6].ab[1];
				VI_data.VV_Num[0] = 6;
				if(VI_data.VV_Num[1] != VI_data.VV_Num[0] && flag_zero != 1 && flag_zero != 2){
					SwCnt_total++;
				}
			}
			// 出力値が初期値に戻ったら零ベクトルを打ち続けるフラグ
			if(VI_data.theta_ref > 5.8 && fabs(VI_data.Pout.ab[0]) < 1.0e-6){
				//Count_enc = 1;
				flag_zero = 2;
			}
			// ベクトル1を打ち続ける処理
			if(flag_zero == 1){
				VI_data.Vout.ab[0] = Inv_param.vout[1].ab[0];
				VI_data.Vout.ab[1] = Inv_param.vout[1].ab[1];
				// VI_data.VV_Num[0] = 1;
				if(VI_data.VV_Num[1] != VI_data.VV_Num[0]){
					SwCnt_total++;
				}
			//}
			// 零ベクトルを打ち続ける
			}else if(flag_zero == 2){
				VI_data.Vout.ab[0] = Inv_param.vout[0].ab[0];
				VI_data.Vout.ab[1] = Inv_param.vout[0].ab[1];
				VI_data.VV_Num[0] = 0;
				if(VI_data.VV_Num[1] != VI_data.VV_Num[0]){
					SwCnt_total++;
				}
			}
			/***************************************************************/

			VI_data.Pout.ab[0] += VI_data.Vout.ab[0]*Ts;
			VI_data.Pout.ab[1] += VI_data.Vout.ab[1]*Ts;

			/******************最終地と初期値をつなげる処理******************/
			// ベクトル6を打ち続けてα軸にあたるときのフラグ
			if(VI_data.theta_ref > 5.8 && VI_data.VV_Num[0] == 6 && fabs(VI_data.Pout.ab[1]) < 1.0e-6){
				flag_zero = 1;
			}
			/**************************************************************/

			// 電圧積分モデルの計算
			VI_data.Perr.ab[0] = VI_data.Pout.ab[0] - VI_data.Pref_a[j*Mul_Reso-Mul_Reso];
			VI_data.Perr.ab[1] = VI_data.Pout.ab[1] - VI_data.Pref_b[j*Mul_Reso-Mul_Reso];
			VI_data.Perr_nrm = sqrt(VI_data.Perr.ab[0]*VI_data.Perr.ab[0]+VI_data.Perr.ab[1]*VI_data.Perr.ab[1]);
			VI_data.Pout_nrm = sqrt(VI_data.Pout.ab[0]*VI_data.Pout.ab[0] + (VI_data.Pout.ab[1]-Kvf)*(VI_data.Pout.ab[1]-Kvf));
			VI_data.Pref_nrm = sqrt(VI_data.Pref_a[j*Mul_Reso-Mul_Reso]*VI_data.Pref_a[j*Mul_Reso-Mul_Reso]+(VI_data.Pref_b[j*Mul_Reso-Mul_Reso]-Kvf)*(VI_data.Pref_b[j*Mul_Reso-Mul_Reso]-Kvf));
			Perr_nrm_total += VI_data.Perr_nrm;

			// VVシーケンスをVV番号と持続回数の情報にエンコードしてファイルに保存
			if(VI_data.VV_Num[0] != VI_data.VV_Num[1]){
				fprintf(fps[1], "%d, %d\n", VI_data.VV_Num[1], Count_enc);
				Count_enc = 1;
			}else {
				Count_enc++;
			}
			VI_data.VV_Num[1] = VI_data.VV_Num[0];

			time_data += Ts;

			//指令値が初期値にきたら終わる
			if(VI_data.theta_ref > 5.8 && VI_data.Pref_a[j*Mul_Reso-Mul_Reso] > Vref_nrm*Ts){
				Count_enc -= 1;
				fprintf(fps[1], "%d, %d\n", VI_data.VV_Num[1], Count_enc);
				Count_enc = 1;
				goto stop;
			}
			fprintf(fps[0], "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d\n", time_data, VI_data.Pout.ab[0], VI_data.Pout.ab[1], VI_data.Pref_a[j*Mul_Reso-Mul_Reso], VI_data.Pref_b[j*Mul_Reso-Mul_Reso], Boundary, VI_data.Perr_nrm, Pout_nrm, Kvf+Boundary, Kvf-Boundary, Pref_nrm, VI_data.theta_ref, VI_data.VV_Num[0]);
		}

		// 次のDFSで使用する変数の情報を渡す
		dfs.VI.Pout.ab[0] = VI_data.Pout.ab[0];
		dfs.VI.Pout.ab[1] = VI_data.Pout.ab[1];
		dfs.VV_Num[0] = dfs.VV_seq_record[DFS_MAX];
		Vector_tmp = dfs.VV_seq_record[DFS_MAX-1];
		vector0_tmp[0] = dfs.VV_seq_record[DFS_MAX];
		vector0_tmp[1] = dfs.VV_seq_record[DFS_MAX-1];
		 
		dfs.SwCnt = 0;
		dfs.VI.Perr_sum = 0;

		for(j = 0; j <= DFS_MAX; ++j){
			dfs.VV_seq[j] = 0;
			dfs.VV_seq_record[j] = 0;
		}
	
		Node_Counter = 0;

		SwCnt_total += dfs.SwCnt_min;
		dfs.SwCnt_min = SwCnt_limit;

		dfs.Layer = 0;

		dfs.VV_seq[0] = dfs.VV_Num[0]; 
		dfs.VV_seq[1] = dfs.VV_Num[0]; 

		// セクター毎にDFSの処理時間を計測して出力
		end_clock = clock();
		Timer_clock[0]=(double)(end_clock-start_clock)/CLOCKS_PER_SEC;
		printf("Loop Time/Loop MAX : %d/%d, Clock_Time : %f\n",Loop_Time, DFS_Loop, Timer_clock[0]);
		printf("Search Number : %d, DFS_MAX : %d\n\n",Node_Counter, DFS_MAX);

	}
	stop:

	SwFreq = (SwCnt_total/6)/(Step_MAX*Ts);

	end_clock2 = clock();
	Timer_clock[1]=(double)(end_clock2-start_clock2)/CLOCKS_PER_SEC;

	printf("DFS Depth           : %d\n", DFS_MAX);
	printf("DFS Loop Times      : %d\n", DFS_Loop);
	printf("Total Steps         : %d\n", Step_MAX);
	printf("Comp. cycle Ts      : %f us\n", Ts * 1e6);
	// printf("Boundary           : %f\n", Boundary);
	printf("Fund. Freq.         : %d Hz\n", f_ref);
	printf("Switching times     : %d\n", SwCnt_total);
	printf("Switching frequency : %f Hz\n", SwFreq);
	printf("VI Norm Limit       : %f Vs\n", Perr_sum_limit);
	printf("VI Error Area       : %f Vs\n", Perr_nrm_total);
	printf("Pattern Update      : %d\n", Pattern_Update);
	printf("CPU Comp. Time      : %f s\n", Timer_clock[1]);


    for(i = 0; i < FileNum; i++){
		// ファイルを閉じる
        fclose(fps[i]);
	}
}


// ********************************************************************** //
// ************************* DFSによるASPSの処理 ************************* //
// ********************************************************************** //
void DFSbasedASPS(){
	uint8_t VV_node;

	// VVノード0～6を探索
	for(VV_node = 0; VV_node < VV_NUM; ++VV_node){

		// i番目の電圧ベクトルノードを参照して出力VIを計算
		dfs.VI.Vout.ab[0] = Inv_param.vout[VV_node].ab[0];
		dfs.VI.Vout.ab[1] = Inv_param.vout[VV_node].ab[1];
		dfs.VI.Pout.ab[0] += Ts*dfs.VI.Vout.ab[0];
		dfs.VI.Pout.ab[1] += Ts*dfs.VI.Vout.ab[1];
		// 現在のレイヤーの指令VIを参照してVI誤差を計算
		dfs.VI.Perr.ab[0] = dfs.VI.Pout.ab[0] - VI_data.Pref_a[dfs.Layer*Mul_Reso];
		dfs.VI.Perr.ab[1] = dfs.VI.Pout.ab[1] - VI_data.Pref_b[dfs.Layer*Mul_Reso];
		// VI誤差ノルムを計算
		dfs.VI.Perr_nrm = sqrt(dfs.VI.Perr.ab[0]*dfs.VI.Perr.ab[0]+dfs.VI.Perr.ab[1]*dfs.VI.Perr.ab[1]);

		// 探索を隣接ベクトルに限定
		Total_Step = dfs.Layer+DFS_MAX_SUM;
		if(AdjoinVector(Total_Step, VV_node, omega_ref, Ts)){
			goto NonRecursive;
		}

		VV_node_tmp[dfs.Layer] = VV_node;  // レイヤごとのノード探索状況を保存
		dfs.Layer++;   				       // レイヤを更新
		dfs.VV_Num[2] = dfs.VV_Num[1];     // 前々回のベクトル番号として保存
		dfs.VV_Num[1] = dfs.VV_Num[0];     // 前回のベクトル番号として保存
		dfs.VV_Num[0] = VV_node;           // 現在のVV番号を更新

		// ZVVの時の処理 (VV : 0 or 7を判定)
		if(dfs.VV_Num[0] == 0 && SwCntTable[7][dfs.VV_Num[1]] < SwCntTable[0][dfs.VV_Num[1]]){
			dfs.VV_Num[0] = 7;
		}
		dfs.VV_seq[dfs.Layer] = dfs.VV_Num[0];  		        // ベクトル番号格納
		dfs.SwCnt += SwCntTable[dfs.VV_Num[0]][dfs.VV_Num[1]];  // スイッチング回数を更新

		// スイッチング回数による枝刈(現在の最小値より大きければ枝刈)
		if(dfs.SwCnt > dfs.SwCnt_min){
			dfs.SwCnt -= SwCntTable[dfs.VV_Num[0]][dfs.VV_Num[1]]; // スイッチング回数を戻す
			dfs.Layer--;										   // レイヤを戻す													
			dfs.VV_Num[0] = dfs.VV_Num[1];						   
			dfs.VV_Num[1] = dfs.VV_Num[2];
			goto NonRecursive;
		}

		dfs.VI.Perr_sum += dfs.VI.Perr_nrm; 		 // VI誤差面積を更新
		Perr_nrm_tmp[dfs.Layer] = dfs.VI.Perr_nrm;   // レイヤごとのVI誤差を保存しておく 

		// 電圧積分誤差面積による枝刈(設定した許容値を超えれば枝刈)
		if(dfs.VI.Perr_sum > Perr_sum_limit){
			dfs.VI.Perr_sum -= dfs.VI.Perr_nrm;					   // 誤差面積を戻す
			dfs.SwCnt -= SwCntTable[dfs.VV_Num[0]][dfs.VV_Num[1]]; // スイッチング回数を戻す		
			dfs.Layer--;										   // レイヤを戻す	
			dfs.VV_Num[0] = dfs.VV_Num[1];
			dfs.VV_Num[1] = dfs.VV_Num[2];
			goto NonRecursive;
		}

		// ****************** 最終レイヤに到達したときの処理 ****************** //
		if(dfs.Layer == DFS_MAX){
			// 最小スイッチング回数が更新されたとき
			if(dfs.SwCnt < dfs.SwCnt_min){
				dfs.SwCnt_min = dfs.SwCnt; 			      // スイッチング回数の最小値を保存
				dfs.VI.Perr_sum_min = dfs.VI.Perr_sum;    // 誤差面積の最小値を更新(これを基準に判定)
				for(j = 0; j <= DFS_MAX; ++j){
					dfs.VV_seq_record[j] = dfs.VV_seq[j]; // VVシーケンスの保存
				}
			}

			// 最小スイッチング回数と同じ場合は誤差面積の大小で評価
			if(dfs.SwCnt == dfs.SwCnt_min && dfs.VI.Perr_sum < dfs.VI.Perr_sum_min){
				for(j = 0; j <= DFS_MAX; ++j){
					dfs.VV_seq_record[j] = dfs.VV_seq[j];  // VVシーケンスの保存
				}
				dfs.VI.Perr_sum_min = dfs.VI.Perr_sum;     // 誤差面積の最小値を更新
			}
			dfs.SwCnt -= SwCntTable[dfs.VV_Num[0]][dfs.VV_Num[1]];  // スイッチング回数を戻す
			dfs.VI.Perr_sum -= dfs.VI.Perr_nrm; 					// 誤差面積を戻す
			Pattern_Update++; 										// 更新回数記録
			dfs.VV_Num[0] = dfs.VV_Num[1];
			dfs.VV_Num[1] = dfs.VV_seq[dfs.Layer-2];
			dfs.Layer--;											// レイヤを戻す
			goto NonRecursive;  			
		}

		// ** 再帰呼び出し ** //
		DFSbasedASPS();

		VV_node = VV_node_tmp[dfs.Layer]; // バックトラックした後のレイヤのノードを代入

		NonRecursive: // 枝刈と最終レイヤの再起呼び出しスキップ場所

		// 前回の積分値に戻す
		dfs.VI.Vout.ab[0] = Inv_param.vout[VV_node].ab[0];
		dfs.VI.Vout.ab[1] = Inv_param.vout[VV_node].ab[1];
		dfs.VI.Pout.ab[0] -= Ts*dfs.VI.Vout.ab[0];  
		dfs.VI.Pout.ab[1] -= Ts*dfs.VI.Vout.ab[1];
	}	

	// **** レイヤの全ノードが処理済み(for文を完了) **** // 
	dfs.SwCnt -= SwCntTable[dfs.VV_Num[1]][dfs.VV_Num[0]];  // スイッチング回数を戻す
	dfs.VI.Perr_sum -= Perr_nrm_tmp[dfs.Layer]; 			// 誤差面積を戻す 
	dfs.VV_Num[0] = dfs.VV_Num[1];
	if(dfs.Layer-2 < 0){
		dfs.VV_Num[1] = Vector_tmp;
	}else{
		dfs.VV_Num[1] = dfs.VV_seq[dfs.Layer-2];
	}
	dfs.Layer--;  // レイヤを戻す

	return;	
}

bool InitDFS(){

}



