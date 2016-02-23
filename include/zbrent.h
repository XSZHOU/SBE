// extracted from roots.h, get rid of nr3
#include <limits>

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
	

template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol)
{
	const int ITMAX=200;//changed from 100
	const double EPS=std::numeric_limits<double>::epsilon();
	double a=x1,b=x2,c=x2,d=0.0,e=0.0,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
       {   //printf("kmin/fa kmax/fb | %e/%e : %e/%e\n",a,fa,b,fb);      // S.Z.
        throw("Root must be bracketed in zbrent");}
	fc=fb; // for the first iteration, then e will get initialized there // S.Z.
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (abs(fc) < abs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*abs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (abs(xm) <= tol1 || fb == 0.0) return b;
		if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=abs(p);
			double min1=3.0*xm*q-abs(tol1*q);
			double min2=abs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (abs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
}
