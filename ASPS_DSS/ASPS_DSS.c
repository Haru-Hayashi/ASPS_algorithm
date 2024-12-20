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

// 指令値データ生成
double Mod;
double omega_ref;
double f_ref;
double Vref_nrm;
double VI_nrm;
double Pref_a_tmp, Pref_b_tmp;

// DFS
int DFS_Loop;
int DFS_MAX;
int Step_MAX;
int Step_Current;
int Step_Total;
int VV_seq_tmp;
int SwCnt_total;
int SwCnt_limit;
int VV_node_tmp[100];
long Pattern_Update;
double Ts;
double SwFreq;
double Boundary;
double Perr_sum_limit;
double Perr_sum_total;
double Perr_sum_ave;
double Perr_nrm_tmp[100];
double VI_theta_tmp;

// データ管理
double time_data;
int header_written;
int VV_Counter;
int Loop_Time;

// 汎用イテレータ
int i, j, k;

// 処理時間
double Timer_clock[2];
clock_t start_clock, end_clock;
clock_t start_clock2, end_clock2;


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
	Mod = Mod_Max;							 // 変調率 
	omega_ref = Omega_rated*Mod/Mod_Max;     // 角周波数指令
	f_ref =omega_ref/(2*pi);			     // 周波数指令

	// ******* 演算周期とDFSの設定 ******* //
	Ts = 25e-6;  // 演算周期
	Step_MAX = floor((double)1/f_ref/Ts);    // 1周期のステップ数
	DFS_MAX = floor(pi/(6*omega_ref*Ts));    // 1周期のステップ数を12分割 
	DFS_Loop = floor(Step_MAX/DFS_MAX) + 1;  // DFSを実行する回数
	// DFS_Loop = 2; 

	// ********* 許容値の設定 ********* //
	Boundary = 12.0e-3;                      // 境界幅 (境界付きで探索するときに使用)
	SwCnt_limit = 8;                         // スイッチング回数の上限値
	dfs.SwCnt_min = SwCnt_limit;             // 現在のスイッチング回数の最小値
	Perr_sum_limit = 0.7;                    // 誤差面積の上限値
          
	VI_nrm = V_DC/sqrt(2)/360;               // 電圧積分ノルム(V/w[Vs/(rad/s)])
	Vref_nrm = VI_nrm*omega_ref;			 // 電圧ノルム(V[V])


// ********************************************************************** //
// ************************** アルゴリズムの記述 ************************** //
// ********************************************************************** //
	start_clock2 = clock();

	// ********** DFSの基づくASPSの処理開始(基本12回ループ) ********** //

	//　インバータの出力可能な電圧ベクトルを計算
	Inv_UpdateOutputVoltage(&Inv_param);

	for(k = 0; k < DFS_Loop; ++k){
		Loop_Time++;
		start_clock = clock();

		// VI指令値の計算とデータ格納 
		// 指令値はTs/Mul_Resoの分解能で計算(Mul_Reso∈Z)
		for(j = 0; j < DFS_MAX * Mul_Reso; ++j){
			Pref_a_tmp += Vref_nrm*cos(VI_theta_tmp)*Ts/Mul_Reso;
			Pref_b_tmp += Vref_nrm*cos(VI_theta_tmp-pi*(0.5))*Ts/Mul_Reso;
			VI_data.Pref_a[j] = Pref_a_tmp;
			VI_data.Pref_b[j] = Pref_b_tmp;
			VI_theta_tmp += omega_ref*Ts/Mul_Reso;
			VI_data.theta_ref[j] = VI_theta_tmp;
		}

		// DFSに基づくASPSの処理関数
		DFSbasedASPS();

		// 現在のトータルステップ数を更新
		Step_Total += DFS_MAX;

		// ******* VVシーケンスからVIを計算してデータファイルに保存 ******* //
		for(j = 1; j <= DFS_MAX; ++j){
			// 出力電圧をデータ構造体に渡す
			VI_data.Vout.ab[0] = Inv_param.vout[dfs.VV_seq_record[j]].ab[0];
			VI_data.Vout.ab[1] = Inv_param.vout[dfs.VV_seq_record[j]].ab[1];
			// VV番号をデータ構造体に渡す
			VI_data.VV_Num[0] = dfs.VV_seq_record[j];

			// 最終値と初期値を(0,0)でつなげる処理(最後のベクトル3本VV:0,1,6に適用)
			Endpoint_adjust(&VI_data, &Inv_param, &SwCnt_total, j);

			// *** 電圧積分モデルの計算 *** //
			// VIの計算結果をデータ構造体に渡す
			VI_data.Pout.ab[0] += VI_data.Vout.ab[0]*Ts;
			VI_data.Pout.ab[1] += VI_data.Vout.ab[1]*Ts;
			VI_data.Perr.ab[0] = VI_data.Pout.ab[0] - VI_data.Pref_a[(j-1)*Mul_Reso];
			VI_data.Perr.ab[1] = VI_data.Pout.ab[1] - VI_data.Pref_b[(j-1)*Mul_Reso];
			VI_data.Perr_nrm = sqrt(VI_data.Perr.ab[0]*VI_data.Perr.ab[0]+VI_data.Perr.ab[1]*VI_data.Perr.ab[1]);
			VI_data.Pout_nrm = sqrt(VI_data.Pout.ab[0]*VI_data.Pout.ab[0] + (VI_data.Pout.ab[1]-VI_nrm)*(VI_data.Pout.ab[1]-VI_nrm));
			VI_data.Pref_nrm = sqrt(VI_data.Pref_a[(j-1)*Mul_Reso]*VI_data.Pref_a[(j-1)*Mul_Reso]+(VI_data.Pref_b[(j-1)*Mul_Reso]-VI_nrm)*(VI_data.Pref_b[(j-1)*Mul_Reso]-VI_nrm));
			VI_data.Perr_sum += VI_data.Perr_nrm;
			Perr_sum_total += VI_data.Perr_nrm;

			// **** VV番号と持続回数の情報をファイルに保存 **** //
			// VV番号が変わったら保存
			if(VI_data.VV_Num[0] != VI_data.VV_Num[1]){
				fprintf(fps[1], "%d, %d\n", VI_data.VV_Num[1], VV_Counter);
				VV_Counter = 1;
			}else {
				VV_Counter++;  // VV持続回数をカウント
			}
			// 前回のVV番号を保存しておく
			VI_data.VV_Num[1] = VI_data.VV_Num[0];

			time_data += Ts;

			// 指令値が初期値(0,0)にきたら終わる
			if(VI_data.theta_ref[(j-1)*Mul_Reso] > pi && VI_data.Pref_a[(j-1)*Mul_Reso] > 1.0e-6){
				VV_Counter -= 1;
				fprintf(fps[1], "%d, %d\n", VI_data.VV_Num[1], VV_Counter);
				VV_Counter = 1;
				goto stop;
			}

			// 計算結果をファイル出力
			if(!header_written){
				fprintf(fps[0],"time, Pout_a, Pout_b, Pref_a, Pref_b, Perr_nrm, Pout_nrm, Pref_nrm, "
								"theta_ref, VV\n");
				header_written = 1;
			}
			fprintf(fps[0], "%f, %f, %f, %f, %f, %f, %f, %f, %f, %d\n", 
			time_data, VI_data.Pout.ab[0], VI_data.Pout.ab[1], VI_data.Pref_a[(j-1)*Mul_Reso], 
			VI_data.Pref_b[(j-1)*Mul_Reso], VI_data.Perr_nrm, VI_data.Pout_nrm, VI_data.Pref_nrm, 
			VI_data.theta_ref[(j-1)*Mul_Reso], VI_data.VV_Num[0]);
		}

		// スイッチング回数の更新
		SwCnt_total += dfs.SwCnt_min;

		// セクター毎に処理時間と探索状況を出力
		end_clock = clock();
		Timer_clock[0]=(double)(end_clock-start_clock)/CLOCKS_PER_SEC;
		printf("Current Sector/Max. Sector : %d/%d, Clock_Time : %f sec\n",Loop_Time, DFS_Loop-1, Timer_clock[0]);
		printf(" DFS Bott. : %d, SwCnt : %d, Error Area: %f Vs\n\n", Pattern_Update, dfs.SwCnt_min, VI_data.Perr_sum);

		// DFS更新回数リセット
		Pattern_Update = 0.0;
		// 誤差面積(セクタ毎)をリセット
		VI_data.Perr_sum = 0.0;

		// 次のセクタのDFS初期条件をセット
		dfs.VI.Pout.ab[0] = VI_data.Pout.ab[0];
		dfs.VI.Pout.ab[1] = VI_data.Pout.ab[1];
		dfs.VV_Num[0] = dfs.VV_seq_record[DFS_MAX];
		dfs.SwCnt_min = SwCnt_limit;
		dfs.SwCnt = 0;
		dfs.VI.Perr_sum = 0;
		dfs.Layer = 0;
		dfs.VV_seq[0] = dfs.VV_Num[0]; 
		dfs.VV_seq[1] = dfs.VV_Num[0]; 
		VV_seq_tmp = dfs.VV_seq_record[DFS_MAX-1];
		
	}
	stop:

	// スイッチング周波数を計算
	SwFreq = (SwCnt_total/6)/(Step_MAX*Ts);

	// セクタの平均誤差面積を計算
	Perr_sum_ave = Perr_sum_total/(DFS_Loop-1);

	// 合計の処理時間を計算
	end_clock2 = clock();
	Timer_clock[1]=(double)(end_clock2-start_clock2)/CLOCKS_PER_SEC;

    // 演算条件を出力
	printf("DFS Depth         : %d\n", DFS_MAX);
	printf("DFS Loop Times    : %d\n", DFS_Loop-1);
	printf("Total Steps       : %d\n", Step_MAX);
	printf("Comp. Cycle Ts    : %f us\n", Ts * 1e6);
	printf("Fund. Freq.       : %f Hz\n", f_ref);
	printf("Switching Count   : %d\n", SwCnt_total);
	printf("Switching Freq.   : %f Hz\n", SwFreq);
	printf("DSS Limit         : %f Vs\n", Perr_sum_limit);
	printf("DSS Average       : %f Vs\n", Perr_sum_ave);
	printf("CPU Comp. Time    : %f s\n", Timer_clock[1]);

	// ファイルを閉じる
    for(i = 0; i < FileNum; i++){
        fclose(fps[i]);
	}
}


// ********************************************************************** //
// ************************* DFSによるASPSの処理 ************************* //
// ********************************************************************** //
void DFSbasedASPS(){
	uint8_t VV_node;  // DFSレイヤのイテレータ

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
		Step_Current = dfs.Layer+Step_Total;
		if(AdjoinVector(Step_Current, VV_node, omega_ref, Ts)){
			goto NonRecursive;
		}

		VV_node_tmp[dfs.Layer] = VV_node;  // レイヤごとのノード探索状況を保存
		dfs.VV_Num[2] = dfs.VV_Num[1];     // 前々回のベクトル番号として保存
		dfs.VV_Num[1] = dfs.VV_Num[0];     // 前回のベクトル番号として保存
		dfs.VV_Num[0] = VV_node;           // 現在のVV番号を更新
		dfs.Layer++;   				       // レイヤを更新

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

		VV_node = VV_node_tmp[dfs.Layer]; // バックトラック後のレイヤのノードを代入

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
		dfs.VV_Num[1] = VV_seq_tmp;
	}else{
		dfs.VV_Num[1] = dfs.VV_seq[dfs.Layer-2];
	}
	dfs.Layer--;  // レイヤを戻す

	return;	
}




