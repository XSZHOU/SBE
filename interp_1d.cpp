#include "interp_1d.h"
#include <cstdio>

// extracted from interp_1d_old.h, get rid of nr3
// declaration and definitions separated by S.Z.
// http://www.nr.com/forum/showthread.php?t=1332

int Base_interp::locate(const double x)
{
	int ju,jm,jl;
	if (n < 2 || mm < 2 || mm > n) { printf("t1\n");throw("locate size error");}
	bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if ((x >= xx[jm]) == ascnd)                        //added parenthese,S.Z.
			jl=jm;
		else
			ju=jm;
	}
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}


int Base_interp::hunt(const double x)
{
	int jl=jsav, jm, ju, inc=1;
	if (n < 2 || mm < 2 || mm > n) {printf("t2\n");throw("hunt size error");}
	bool ascnd=(xx[n-1] >= xx[0]);
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1;
	} else {
		if ((x >= xx[jl]) == ascnd) {                    //added parenthese,S.Z.
			for (;;) {
				ju = jl + inc;
				if (ju >= n-1) { ju = n-1; break;}
				else if ((x < xx[ju]) == ascnd) break;   //added parenthese,S.Z.
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if ((x >= xx[jl]) == ascnd) break;  //added parenthese,S.Z.
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if ((x >= xx[jm]) == ascnd)                      //added parenthese,S.Z.
			jl=jm;
		else
			ju=jm;
	}
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

void Spline_interp::sety2(const double *xv, const double *yv, double yp1, double ypn)
{
	int i,k;
	double p,qn,sig,un;
	int n=y2.size();                                     // y2(xv.size), bug fix
    std::vector<double> u(n-1);
	if (yp1 > 0.99e99)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
		u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e99)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}

double Spline_interp::rawinterp(int jl, double x)
{
	int klo=jl,khi=jl+1;
	double y,h,b,a;
	h=xx[khi]-xx[klo];
	if (h == 0.0) {printf("t5 %e %e\n",x,xx[khi]);throw("Bad input to routine splint");}
	a=(xx[khi]-x)/h;
	b=(x-xx[klo])/h;
	y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]
		+(b*b*b-b)*y2[khi])*(h*h)/6.0;
	return y;
}

