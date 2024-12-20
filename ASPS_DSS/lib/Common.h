#ifndef COMMOM_H
#define COMMOM_H

// *** 変数 *** //
// 円周率
#define pi   (3.141592653589793)
// 1/2*pi
#define PI12 (1.570796326794894)
// 2/3*pi
#define PI23 (2.094395102393195)
// 1/6*pi
#define PI16 (0.523598775598298)
// 5/6*pi
#define PI56 (2.617993877991494)
// sqrt(2/3)
#define SQ23 (0.816496580927726)
// sqrt(3/2)
#define SQ32 (1.224744871391589)
//

// *** マクロ関数 *** //
#define PI(n)   ((float)M_PI*(n))
// 直交座標(x, y)からノルムを求める
#define	SQRT2(x,y)  (sqrt((x)*(x)+(y)*(y)))
// 符号関数
#define SGN(x)  (((x)>(0.0))?(1.0):(((x)<(0.0))?(-1.0):(0.0)))
// 絶対値
#define ABS(x)  (((x)>=(0)?(x):(-(x))))
// 四捨五入
#define RND(x)  ((int)((x)+0.5))

#endif