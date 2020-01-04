#include<DIDiDa.h>
int main()
{
    CART_POS a,b,c,d,e,f,u,m,v;
    double angle;
    BOOL bool = 0;
    CART_POS pos0;
    CART_POS pos1;

if(bool == 0){
//    printf("circenter----------\na:\n");
//    scanf("%ld %ld %ld",&a.x,&a.y,&a.z);
//    printf("b:\n");
//    scanf("%ld %ld %ld",&b.x,&b.y,&b.z);
//    printf("c:\n");
//    scanf("%ld %ld %ld",&c.x,&c.y,&c.z);

//    printf("d:\n");
//    scanf("%ld %ld %ld",&d.x,&d.y,&d.z);
//    printf("e:\n");
//    scanf("%ld %ld %ld",&e.x,&e.y,&e.z);
//    printf("f:\n");
//    scanf("%ld %ld %ld",&f.x,&f.y,&f.z);
    a.x = 994665;
    a.y = -233402;
    a.z = -363114;

    b.x = 931479;
    b.y = -281748;
    b.z = -362134;

    c.x = 994544;
    c.y = -353073;
    c.z = -362951;

    d.x = 994517;
    d.y = -248477;
    d.z = -359280;

    e.x = 929340;
    e.y = -281675;
    e.z = -362154;

    f.x = 994427;
    f.y = -365966;
    f.z = -358534;
    pos0 = ObtWorCtr(a,b,c);
    pos1 = ObtWorCtr(d,e,f);
}
else {
    a.x = 1035199;
    a.y = -4049;
    a.z = -456173;

    b.x = 1094319;
    b.y = -3173;
    b.z = -454329;

    c.x = 1211561;
    c.y = -46205;
    c.z = -463763;

    d.x = 1035408;
    d.y = -5353;
    d.z = -454927;

    e.x = 1094964;
    e.y = 1722;
    e.z = -451149;

    f.x = 1232836;
    f.y = -46916;
    f.z = -459181;

    angle = 78;

    pos0 = ObtWorVer(a,b,c,angle);
    pos1 = ObtWorVer(d,e,f,angle);
}
    CART_POS offset = Theoffset(pos0,pos1);
    printf("The offset is:\nx = %ld\ny = %ld\nz = %ld\nRx = %ld\nRy = %ld\nRz = %ld\n",offset.x,offset.y,offset.z,offset.Rx,offset.Ry,offset.Rz);
    system("pause");
    return 0;
}


