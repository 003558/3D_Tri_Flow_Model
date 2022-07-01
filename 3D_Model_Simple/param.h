// -*-c-*-
/**************************

 計算用パラメータ定数の設定

 Last Update: 2017.5.16

***************************/
#define STEP_MAX 3600*240 //計算終了時のstep数
//#define STEP_MAX 7200*72//計算終了時のstep数
#define DAT_OUT 3600//計算結果を出力するstep間隔
//#define DAT_OUT 7200//計算結果を出力するstep間隔
#define DT ((DATA_TYPE)1.0) //1stepの計算時間幅(sec)
//#define DT ((DATA_TYPE)0.5) //1stepの計算時間幅(sec)
#define DDT ((DATA_TYPE)3600.0)
//#define DDT ((DATA_TYPE)7200.0)

#define NOZ 50 //z方向の分割数
#define NOB 240//97 //境界条件読込み数

#define Zmax ((DATA_TYPE)0.0)
#define Zmin ((DATA_TYPE)-10.0)
#define ZZ ((DATA_TYPE)10.0)

#define GV ((DATA_TYPE)9.81)
#define Kappa ((DATA_TYPE)0.4)

#define N_mol ((DATA_TYPE)0.000001)
#define C_1 ((DATA_TYPE)1.44)
#define C_2 ((DATA_TYPE)1.92)
#define C_3 ((DATA_TYPE)1.0)
#define C_mu ((DATA_TYPE)0.09)
#define C_d ((DATA_TYPE)0.0005)
//#define C_d ((DATA_TYPE)0.005)

#define S_r ((DATA_TYPE)0.8)
#define S_k ((DATA_TYPE)1.0)
#define S_e ((DATA_TYPE)1.3)

#define R_a ((DATA_TYPE)1.366)
#define R_w ((DATA_TYPE)1000.0)
#define R_s ((DATA_TYPE)1000.0)

#define Smin ((DATA_TYPE)0.0)
#define Smax ((DATA_TYPE)33.0)

#define D_h ((DATA_TYPE)500.0)

#define M_PI 3.14159265358979323846 /* pi 2018.04.17 harada Visula Studioでコンパイルする際に必要 */