#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<Windows.h>

#define N 4
#define PI 3.1415926

typedef struct cart_pos {
    long x; //单位0.001mm
    long y;
    long z;
    long Rx; //单位0.0001deg
    long Ry;
    long Rz;
}CART_POS;

//取平方
double getSquart(double a);

//四舍五入
long round_fun(double a);

//两点间距离
long getdis(CART_POS a,CART_POS b);

//求4*4逆矩阵
void getAniMat(double A[][N], double B[][N], int n);

//矩阵相乘
void multi(double *A, double *B, double *C,int num);

//获取平面法向量
void getvector(CART_POS a, CART_POS b, CART_POS c);

//3*3旋转矩阵
void RotationMatrix3X3(double a, double b, double c, double R[3][3]);

//用户坐标转世界坐标
CART_POS UserToBt(CART_POS a, CART_POS b);

//世界坐标转用户坐标
CART_POS BtToUser(CART_POS a, CART_POS b);

//获得圆心坐标
CART_POS getcenter(CART_POS a, CART_POS b, CART_POS c);

//求姿态和圆心坐标
CART_POS ObtWorCtr(CART_POS a, CART_POS b, CART_POS c);

//获得工件顶点
CART_POS ObtWorVer(CART_POS a, CART_POS b, CART_POS c, double age);
//求偏移量
CART_POS Theoffset(CART_POS a, CART_POS b);
