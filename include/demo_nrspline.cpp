// demo for an nr3-independent interp_1d routine.
#include "interp_1d.h"
#include <cstdio>

// require:  interp_1d.h  interp_1d.cpp
// declaration and definitions separated by S.Z.
// http://www.nr.com/forum/showthread.php?t=1332

// compile with
//g++ -o demo_nrspline interp_1d.cpp demo_nrspline.cpp

int main()
{ 
   int n=5;
   std::vector<double> X(n), Y(n);
   X[0]=0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=2.0;
   Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;
   
    Spline_interp myfunc(X,Y);
    
    double x,y;
    x=1.5;
    y = myfunc.interp(x);
   
    printf("spline at %f is %f\n", x, y);

    return 0;
}
