#ifndef INTERP_1D_H
#define INTERP_1D_H

#include <vector>
#include <cmath>   //for pow() and abs()

// extracted from  interp_1d_old.h. get rid of nr3
// declaration and definitions separated by S.Z.
// http://www.nr.com/forum/showthread.php?t=1332

template<class T>
inline T SQR(const T a) {return a*a;}


template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}
        
        
struct Base_interp
{
	int n, mm, jsav, cor, dj;
	const double *xx, *yy;
	Base_interp(){}              //default ctor add by S.Z.
	Base_interp(const std::vector<double> &x, const double *y, int m)
		: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
		dj = MAX(1,(int)pow((double)n,0.25));
		                         //dj = MIN(1,(int)pow((double)n,0.25)); bug fix
	}

	double interp(double x) {
		int jlo = cor ? hunt(x) : locate(x);
		return rawinterp(jlo,x);
	}

	int locate(const double x);
	int hunt(const double x);
	
	double virtual rawinterp(int jlo, double x) = 0;
    
    virtual ~Base_interp(){}      
	// added S.Z. virtual destructor
    // This is Item 7 in Scott Meyers' Effective C++ . 
    // http://stackoverflow.com/questions/461203/when-to-use-virtual-destructors
};


struct Spline_interp : Base_interp
{
	std::vector<double> y2;
	Spline_interp(){}            //added by S.Z
	Spline_interp(const std::vector<double> &xv, const std::vector<double> &yv, double yp1=1.e99, double ypn=1.e99)
	: Base_interp(xv,&yv[0],2), y2(xv.size())
	{sety2(&xv[0],&yv[0],yp1,ypn);}

	Spline_interp(const std::vector<double> &xv, const double *yv, double yp1=1.e99, double ypn=1.e99)
	: Base_interp(xv,yv,2), y2(xv.size())
	{sety2(&xv[0],yv,yp1,ypn);}

	void sety2(const double *xv, const double *yv, double yp1, double ypn);
	double rawinterp(int jl, double xv);
};


#endif
