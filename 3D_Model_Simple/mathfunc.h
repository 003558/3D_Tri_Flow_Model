// -*-c-*-
/***********************************

 浮動小数点型と数学関数の定義

 Last Update: 2017.5.16

***********************************/

#define TYPE 1 //浮動小数点型[0:float,1:double]

#if TYPE == 0
#define DATA_TYPE float
#define POW powf //指数で累乗した値を返す(float型)
#define FABS fabsf //絶対値を返す(float型)
#elif TYPE == 1
#define DATA_TYPE double
#define POW pow //指数で累乗した値を返す(double型)
#define FABS fabs //絶対値を返す(double型)
#endif

#define MAX(a,b) (((a)>(b)) ? (a) : (b)) //最大値を返す
#define MIN(a,b) (((a)<(b)) ? (a) : (b)) //最小値を返す
#define SGN(a) ((fabs(a)>0) ? (fabs(a)/(a)) : 0) //符号を返す
