#include <DIDiDa.h>


//取平方
double getSquart(double a){
    return a*a;
}
//四舍五入
long round_fun(double a){
    long b;
    if(a>0){
        b = a+0.5;
    }
    else {
        b = a-0.5;
    }
    return b;
}

//两点间距离
long getdis(CART_POS a,CART_POS b){
    double distance = (double)(b.x-a.x)*(b.x-a.x)+(double)(b.y-a.y)*(b.y-a.y)+(double)(b.z-a.z)*(b.z-a.z);
    return (long)sqrt(distance);
}

//求4*4逆矩阵
void getAniMat(double A[][N], double B[][N], int n)
{
    int i, j, k;
    double max, temp;
    double t[N][N];                //临时矩阵
    //将A矩阵存放在临时矩阵t[n][n]中
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            t[i][j] = A[i][j];
        }
    }
    //初始化B矩阵为单位阵
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = (i == j) ? (double)1 : 0;
        }
    }
    for (i = 0; i < n; i++)
    {
        //寻找主元
        max = t[i][i];
        k = i;
        for (j = i + 1; j < n; j++)
        {
            if (fabs(t[j][i]) > fabs(max))
            {
                max = t[j][i];
                k = j;
            }
        }
        //如果主元所在行不是第i行，进行行交换
        if (k != i)
        {
            for (j = 0; j < n; j++)
            {
                temp = t[i][j];
                t[i][j] = t[k][j];
                t[k][j] = temp;
                //B伴随交换
                temp = B[i][j];
                B[i][j] = B[k][j];
                B[k][j] = temp;
            }
        }
        //判断主元是否为0, 若是,则矩阵A不是满秩矩阵,不存在逆矩阵
        if (t[i][i] == 0)
        {
            printf("There is no inverse matrix!");
        }
        //消去A的第i列除去i行以外的各行元素
        temp = t[i][i];
        for (j = 0; j < n; j++)
        {
            t[i][j] = t[i][j] / temp;        //主对角线上的元素变为1
            B[i][j] = B[i][j] / temp;        //伴随计算
        }
        for (j = 0; j < n; j++)        //第0行->第n行
        {
            if (j != i)                //不是第i行
            {
                temp = t[j][i];
                for (k = 0; k < n; k++)        //第j行元素 - i行元素*j列i行元素
                {
                    t[j][k] = t[j][k] - t[i][k] * temp;
                    B[j][k] = B[j][k] - B[i][k] * temp;
                }
            }
        }
    }
}

//矩阵相乘
void multi(double *A, double *B, double *C,int num){
    int i, j, k;
    double temp = 0;
    for(i = 0; i < num; i++){
        for(k = 0; k< num;k++){
            for(j = 0; j < num; j++){//当前行的每个元素
                temp += (*(A+i*num+j))*(*(B+j*num+k));
          }
         (*(C+i*num+k))=temp;
         temp = 0;
        }
    }
}

void getvector(CART_POS a, CART_POS b, CART_POS c)//获取平面法向量
{
    CART_POS v1,v2,vec;
    v1.x = b.x - a.x;
    v1.y = b.y - a.y;
    v1.z = b.z - a.z;
    printf("v1.x, v1.y, v1.z = %ld %ld %ld\n", v1.x, v1.y, v1.z);

    v2.x = c.x - a.x;
    v2.y = c.y - a.y;
    v2.z = c.z - a.z;
    printf("v2.x, v2.y, v2.z = %ld %ld %ld\n", v2.x, v2.y, v2.z);

    vec.x = v1.y * (double)v2.z - v1.z * v2.y;
    vec.y = v1.z * (double)v2.x - v1.x * v2.z;
    vec.z = v1.x * (double)v2.y - v1.y * v2.x;
    printf("v2x, v2y, v2z = %ld %ld %ld\n", vec.x, vec.y, vec.z);
    printf("The normal vector: (%lf  %lf  %lf)\n", vec.x*(double)(b.x - a.x), vec.y*(double)(b.y - a.y), vec.z*(double)(b.z - a.z));
}
//3*3旋转矩阵
void RotationMatrix3X3(double a, double b, double c, double R[3][3])
{
 /* function:将角度值转换为3x3的旋转矩阵 */
 double alpha = a * PI / 180;
 double belta= b * PI / 180;
 double theta= c * PI / 180;

 //分别绕Z,Y,X旋转theta,belta,alpha
 double A[3][3] = { {cos(theta), -sin(theta), 0}, {sin(theta), cos(theta), 0}, {0, 0, 1 } };
 double B[3][3] = { {cos(belta), 0, sin(belta)}, {0, 1, 0}, {-sin(belta) , 0, cos(belta)} };
 double C[3][3] = { {1, 0, 0}, {0, cos(alpha), -sin(alpha)}, {0, sin(alpha), cos(alpha)} };

 double temp[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
 multi((double *)A, (double *)B, (double *)temp,N-1);
 multi((double *)temp, (double *)C, (double *)R,N-1);
}

//UserToBt
CART_POS UserToBt(CART_POS a, CART_POS b){
    CART_POS offset;

    double n[N-1][N-1];
    double b_a[N-1][N-1];
    double b_BT[N-1][N-1];
    double m[1][3] = {{b.x,b.y,b.z}};
    double xyz[3]={0};

    RotationMatrix3X3((a.Rx)/10000.0,a.Ry/10000.0,a.Rz/10000.0,n);
    for (int i=0;i<3;i++) {
        int k =0;
        for (int j=0;j<3;j++) {
            xyz[i]+=m[k][j]*n[i][j];
        }
    }
    RotationMatrix3X3((b.Rx)/10000.0,b.Ry/10000.0,b.Rz/10000.0,b_a);
    multi((double *)n, (double *)b_a, (double *)b_BT,N-1);

    double alp0 = atan2(-b_BT[2][0],sqrt(b_BT[0][0]*b_BT[0][0]+b_BT[1][0]*b_BT[1][0]))*180/3.1415926;
    double alp1 = atan2(b_BT[1][0]/cos(alp0*3.1415926/180),b_BT[0][0]/cos(alp0*3.1415926/180))*180/3.1415926;
    double alp2 = atan2(b_BT[2][1]/cos(alp0*3.1415926/180),b_BT[2][2]/cos(alp0*3.1415926/180))*180/3.1415926;

    offset.x = a.x+xyz[0];
    offset.y = a.y+xyz[1];
    offset.z = a.z+xyz[2];
    offset.Rx = alp2*10000;
    offset.Ry = alp0*10000;
    offset.Rz = alp1*10000;
    return offset;
}
//BtToUser
CART_POS BtToUser(CART_POS a, CART_POS b){
    CART_POS offset;

    double n[N-1][N-1];
    double _n[N][N];
    double b_BT[N-1][N-1];
    double b_a[N-1][N-1];
    double m[1][3] = {b.x-a.x,b.y-a.y,b.z-a.z};
    double xyz[3]={0};

    //a->BT
    RotationMatrix3X3(a.Rx/10000.0,a.Ry/10000.0,a.Rz/10000.0,n);
    double n1[N][N] = {
        {n[0][0], n[0][1], n[0][2], 1.0},
        {n[1][0], n[1][1], n[1][2], 1.0},
        {n[2][0], n[2][1], n[2][2], 1.0},
        {0.0, 0.0, 0.0, 1.0}
    };
    getAniMat(n1,_n,N);
    double n2[N-1][N-1] = {
        {_n[0][0], _n[0][1], _n[0][2]},
        {_n[1][0], _n[1][1], _n[1][2]},
        {_n[2][0], _n[2][1], _n[2][2]}
    };

    //计算b在a坐标系中XYZ
    for (int i=0;i<3;i++) {
        int k =0;
        for (int j=0;j<3;j++) {
            xyz[i]+=n[j][i]*m[k][j];
        }
    }

    //b->BT
    RotationMatrix3X3((b.Rx)/10000.0,b.Ry/10000.0,b.Rz/10000.0,b_BT);
    multi((double *)n2, (double *)b_BT, (double *)b_a,N-1);

    double alp0 = atan2(-b_a[2][0],sqrt(b_a[0][0]*b_a[0][0]+b_a[1][0]*b_a[1][0]))*180/3.1415926;
    double alp1 = atan2(b_a[1][0]/cos(alp0*3.1415926/180),b_a[0][0]/cos(alp0*3.1415926/180))*180/3.1415926;
    double alp2 = atan2(b_a[2][1]/cos(alp0*3.1415926/180),b_a[2][2]/cos(alp0*3.1415926/180))*180/3.1415926;

    offset.x = xyz[0];
    offset.y = xyz[1];
    offset.z = xyz[2];
    offset.Rx = alp2*10000;
    offset.Ry = alp0*10000;
    offset.Rz = alp1*10000;
    return offset;
}

//获得圆心坐标
CART_POS getcenter(CART_POS a, CART_POS b, CART_POS c)//获取圆心坐标及半径，参数1，2，3为三点坐标，4为圆心坐标，5为半径
{
    CART_POS center;
    float v1 = (b.x - a.x)* (c.x - a.x) + (b.y - a.y)* (c.y - a.y) + (b.z - a.z)* (c.z - a.z);
    float v2 = (b.x - a.x)* (b.x - a.x) + (b.y - a.y)* (b.y - a.y) + (b.z - a.z)* (b.z - a.z);
    float v3 = (c.x - a.x)* (c.x - a.x) + (c.y - a.y)* (c.y - a.y) + (c.z - a.z)* (c.z - a.z);

    //验证是否共线
    if ((v1*v1)/(v2*v3) == 1)
    {
        printf("error: collinear\n");
    }

    double x1 = a.x;
    double x2 = b.x;
    double x3 = c.x;
    double y1 = a.y;
    double y2 = b.y;
    double y3 = c.y;
    double z1 = a.z;
    double z2 = b.z;
    double z3 = c.z;
    double ABCD1[4];
    double ABCD2[4];
    double ABCD3[4];

    ABCD1[0] = y1 * z2 - y1 * z3 - z1 * y2 + z1 * y3 + y2 * z3 - y3 * z2;
    ABCD1[1] = -x1 * z2 + x1 * z3 + z1 * x2 - z1 * x3 - x2 * z3 + x3 * z2;
    ABCD1[2] = x1 * y2 - x1 * y3 - y1 * x2 + y1 * x3 + x2 * y3 - x3 * y2;
    ABCD1[3] = -x1 * y2 * z3 + x1 * y3 * z2 + x2 * y1 * z3 - x3 * y1 * z2 - x2 * y3 * z1 + x3 * y2 * z1;

    ABCD2[0] = 2 * (x2 - x1);
    ABCD2[1] = 2 * (y2 - y1);
    ABCD2[2] = 2 * (z2 - z1);
    ABCD2[3] = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

    ABCD3[0] = 2 * (x3 - x1);
    ABCD3[1] = 2 * (y3 - y1);
    ABCD3[2] = 2 * (z3 - z1);
    ABCD3[3] = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

    double mtx[4][4] = {
        {ABCD1[0],ABCD1[1],ABCD1[2],0.0},
        {ABCD2[0],ABCD2[1],ABCD2[2],0.0},
        {ABCD3[0],ABCD3[1],ABCD3[2],0.0},
        {0.0,0.0,0.0,0.0},
    };

    double mtxi[4][4];
    double mtxi1[3][3];
    getAniMat(mtx,mtxi,N-1);
    multi((double*)mtx,(double*)mtxi,(double*)mtxi1,N-1);

    center.x = -(mtxi[0][0] * ABCD1[3] + mtxi[0][1] * ABCD2[3] + mtxi[0][2] * ABCD3[3]);
    center.y = -(mtxi[1][0] * ABCD1[3] + mtxi[1][1] * ABCD2[3] + mtxi[1][2] * ABCD3[3]);
    center.z = -(mtxi[2][0] * ABCD1[3] + mtxi[2][1] * ABCD2[3] + mtxi[2][2] * ABCD3[3]);
    return center;
}

//求三点偏移量和圆心坐标
CART_POS ObtWorCtr(CART_POS a, CART_POS b, CART_POS c){
    //3点求圆心x,y,z
    CART_POS center = getcenter(a, b, c);
    //3点求Rx,Ry,Rz
    CART_POS pos;
        double i,j,k;
        i = b.x-a.x;
        j = b.y-a.y;
        k = b.z-a.z;
        double abc = (getSquart(i)+getSquart(j)+getSquart(k));
        if(abc == 0.0){
            printf("please input the corrent point!");
        }
        else {
            double ansx=(getSquart(i)*c.x+(getSquart(j)+getSquart(k))*a.x-i*j*(a.y-c.y)-i*k*(a.z-c.z))/abc;
            double ansy=(i*j*(c.x-a.x)-j*j*(a.y-c.y)-j*k*(a.z-c.z))/abc+a.y;
            double ansz=(i*k*(c.x-a.x)-k*j*(a.y-c.y)-k*k*(a.z-c.z))/abc+a.z;
            double abc1 = sqrt(getSquart(c.x-ansx)+getSquart(c.y-ansy)+getSquart(c.z-ansz));

            //法向量求RXRYRZ
            double vecx = (b.y - a.y) * (double)(c.z - a.z) - (b.z - a.z) * (double)(c.y - a.y);
            double vecy = (b.z - a.z) * (double)(c.x - a.x) - (b.x - a.x) * (double)(c.z - a.z);
            double vecz = (b.x - a.x) * (double)(c.y - a.y) - (b.y - a.y) * (double)(c.x - a.x);
            double abc2 = sqrt(getSquart(vecx)+getSquart(vecy)+getSquart(vecz));

            //新坐标系矩阵
            double Rot_metri[N-1][N-1] = {
                {(double)(b.x-a.x)/sqrt(abc),(double)(c.x-ansx)/abc1,vecx/abc2},
                {(double)(b.y-a.y)/sqrt(abc),(double)(c.y-ansy)/abc1,vecy/abc2},
                {(double)(b.z-a.z)/sqrt(abc),(double)(c.z-ansz)/abc1,vecz/abc2},
            };
            double angRy = atan2(-Rot_metri[2][0],sqrt(Rot_metri[0][0]*Rot_metri[0][0]+Rot_metri[1][0]*Rot_metri[1][0]));
            double angRz = atan2(Rot_metri[1][0]/cos(angRy),Rot_metri[0][0]/cos(angRy));
            double angRx = atan2(Rot_metri[2][1]/cos(angRy),Rot_metri[2][2]/cos(angRy));
            pos.x = round_fun(center.x);
            pos.y = round_fun(center.y);
            pos.z = round_fun(center.z);
            pos.Rz = round_fun(angRz*10000*180/PI);
            pos.Ry = round_fun(angRy*10000*180/PI);
            pos.Rx = round_fun(angRx*10000*180/PI);
            printf("The pos is:\nx = %ld\ny = %ld\nz = %ld\nRx = %ld\nRy = %ld\nRz = %ld\n",pos.x,pos.y,pos.z,pos.Rx,pos.Ry,pos.Rz);
            return pos;
  }
}

//获得工件顶点
CART_POS ObtWorVer(CART_POS a, CART_POS b, CART_POS c, double age){
    CART_POS pos;
        double i,j,k;
        i = b.x-a.x;
        j = b.y-a.y;
        k = b.z-a.z;
        double abc = (getSquart(i)+getSquart(j)+getSquart(k));
        if(abc == 0.0){
            printf("please input the corrent point!");
        }
        else {
            double ansx=(getSquart(i)*c.x+(getSquart(j)+getSquart(k))*a.x-i*j*(a.y-c.y)-i*k*(a.z-c.z))/abc;
            double ansy=(i*j*(c.x-a.x)-j*j*(a.y-c.y)-j*k*(a.z-c.z))/abc+a.y;
            double ansz=(i*k*(c.x-a.x)-k*j*(a.y-c.y)-k*k*(a.z-c.z))/abc+a.z;
            double abc1 = sqrt(getSquart(c.x-ansx)+getSquart(c.y-ansy)+getSquart(c.z-ansz));

            //法向量求RXRYRZ
            double vecx = (b.y - a.y) * (double)(c.z - a.z) - (b.z - a.z) * (double)(c.y - a.y);
            double vecy = (b.z - a.z) * (double)(c.x - a.x) - (b.x - a.x) * (double)(c.z - a.z);
            double vecz = (b.x - a.x) * (double)(c.y - a.y) - (b.y - a.y) * (double)(c.x - a.x);
            double abc2 = sqrt(getSquart(vecx)+getSquart(vecy)+getSquart(vecz));

            //新坐标系矩阵
            double Rot_metri[N-1][N-1] = {
                {(double)(b.x-a.x)/sqrt(abc),(double)(c.x-ansx)/abc1,vecx/abc2},
                {(double)(b.y-a.y)/sqrt(abc),(double)(c.y-ansy)/abc1,vecy/abc2},
                {(double)(b.z-a.z)/sqrt(abc),(double)(c.z-ansz)/abc1,vecz/abc2},
            };
            double angRy = atan2(-Rot_metri[2][0],sqrt(Rot_metri[0][0]*Rot_metri[0][0]+Rot_metri[1][0]*Rot_metri[1][0]));
            double angRz = atan2(Rot_metri[1][0]/cos(angRy),Rot_metri[0][0]/cos(angRy));
            double angRx = atan2(Rot_metri[2][1]/cos(angRy),Rot_metri[2][2]/cos(angRy));
            pos.x = round_fun(ansx);
            pos.y = round_fun(ansy);
            pos.z = round_fun(ansz);
            pos.Rz = round_fun(angRz*10000*180/PI);
            pos.Ry = round_fun(angRy*10000*180/PI);
            pos.Rx = round_fun(angRx*10000*180/PI);
    if(age == 90.0){
        printf("angle is %lf:\nx = %ld\ny = %ld\nz = %ld\nRx = %ld\nRy = %ld\nRz = %ld\n",age,pos.x,pos.y,pos.z,pos.Rx,pos.Ry,pos.Rz);
        return pos;
    }
    else if(age>90.0){
        long dis1 = getdis(pos,c);
        long dis = -dis1/tan(age*PI/180);
        pos.x = pos.x-dis*cos(angRy)*cos(angRz);
        pos.y = pos.y-dis*cos(angRy)*sin(angRz);
        pos.z = pos.z-dis*cos(angRy)*sin(angRy);
        printf("angle is %lf:\nx = %ld\ny = %ld\nz = %ld\nRx = %ld\nRy = %ld\nRz = %ld\n",age,pos.x,pos.y,pos.z,pos.Rx,pos.Ry,pos.Rz);
        return pos;
    }
    else {
        long dis1 = getdis(pos,c);
        long dis = dis1/tan(age*PI/180);
        pos.x = pos.x+dis*cos(angRy)*cos(angRz);
        pos.y = pos.y+ dis*cos(angRy)*sin(angRz);
        pos.z = pos.z+dis*sin(angRy);
        printf("angle is %lf:\nx = %ld\ny = %ld\nz = %ld\nRx = %ld\nRy = %ld\nRz = %ld\n",age,pos.x,pos.y,pos.z,pos.Rx,pos.Ry,pos.Rz);
        return pos;
    }
  }
}

CART_POS Theoffset(CART_POS a, CART_POS b){
    CART_POS offset;
    double n[N-1][N-1];
    double m[N-1][N-1];

    RotationMatrix3X3((a.Rx)/10000.0,a.Ry/10000.0,a.Rz/10000.0,n);
    RotationMatrix3X3(b.Rx/10000.0,b.Ry/10000.0,b.Rz/10000.0,m);
    double n1[N][N] = {
        {n[0][0], n[0][1], n[0][2], (double)a.x},
        {n[1][0], n[1][1], n[1][2], (double)a.y},
        {n[2][0], n[2][1], n[2][2], (double)a.z},
        {0.0, 0.0, 0.0, 1.0}
    };
    double m1[N][N] = {
        {m[0][0], m[0][1], m[0][2], (double)b.x},
        {m[1][0], m[1][1], m[1][2], (double)b.y},
        {m[2][0], m[2][1], m[2][2], (double)b.z},
        {0.0, 0.0, 0.0, 1.0}
    };

    double _n[N][N];
    getAniMat(n1,_n,N);
    double xyzOffset[N][N];
    multi((double*)m1,(double*)_n,(double*)xyzOffset,N);

    offset.x = round_fun(xyzOffset[0][3]);
    offset.y = round_fun(xyzOffset[1][3]);
    offset.z = round_fun(xyzOffset[2][3]);

    double alp0 = atan2(-xyzOffset[2][0],sqrt(xyzOffset[0][0]*xyzOffset[0][0]+xyzOffset[1][0]*xyzOffset[1][0]));
    double alp1 = atan2(xyzOffset[1][0]/cos(alp0),xyzOffset[0][0]/cos(alp0));
    double alp2 = atan2(xyzOffset[2][1]/cos(alp0),xyzOffset[2][2]/cos(alp0));
    offset.Rx = round_fun(alp2/PI*180*10000);
    offset.Ry = round_fun(alp0/PI*180*10000);
    offset.Rz = round_fun(alp1/PI*180*10000);
    return offset;
}
