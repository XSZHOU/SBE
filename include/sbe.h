#ifndef SBE_H
#define SBE_H


#include <vector>
#include <complex>
#include "interp_1d.h"



/* c++ template functions */
template <typename T> 
inline T sgn(T val) {return ( (val > 0) - (val < 0));}

template<typename T> 
inline T sqr(T x) { return x*x; }


class CInfPass_SBEHF // for SBE @HF level
{
public:
    Spline_interp gPPEc;
    Spline_interp gPPEv;
    Spline_interp gPPVcc;
    Spline_interp gPPVvv;
    Spline_interp gPPVcv;
    double arr[3];
    
    CInfPass_SBEHF ( const std::vector<double> &kgrid, const std::vector<double> &Ec,  const std::vector<double> &Ev ,  
	                 const std::vector<double> &Vcc, const std::vector<double> &Vvv, const std::vector<double> &Vcv) : 
	                 gPPEc(kgrid, Ec), gPPEv(kgrid, Ev), gPPVcc(kgrid, Vcc), gPPVvv(kgrid,Vvv), gPPVcv(kgrid,Vcv) {}
	
};

class CInfPass_SBECC // for SBE @2nd Born level
{
public:
	//direct access vectors
    std::vector<double> EC;
    std::vector<double> EV;
    std::vector<double> KGRID;
    //objects HF
    Spline_interp gPPEc;
    Spline_interp gPPEv;
    Spline_interp gPPVcc;
    Spline_interp gPPVvv;
    Spline_interp gPPVcv;
    //objects CC
    Spline_interp gPPEc2k;
    Spline_interp gPPEv2k;
	Spline_interp gPPdkdEc;
    Spline_interp gPPdkdEv;
    Spline_interp gPPWcc;
    Spline_interp gPPWvv;
    Spline_interp gPPWcv;
    
    double arr[10];
                     
    CInfPass_SBECC(const std::vector<double> &kgrid, const std::vector<double> &Ec,  const std::vector<double> &Ev,
                   const std::vector<double> &ECm, const std::vector<double> &EVm,  const std::vector<double> &dkdEc,
                   const std::vector<double> &dkdEv, const std::vector<double> &Vcc,  const std::vector<double> &Vvv, 
				   const std::vector<double> &Vcv, const std::vector<double> &Wcc,  const std::vector<double> &Wvv, 
				   const std::vector<double> &Wcv) : EC(Ec), EV(Ev), KGRID(kgrid),
                   gPPEc(kgrid, Ec), gPPEv(kgrid, Ev),gPPVcc(kgrid, Vcc), gPPVvv(kgrid,Vvv), gPPVcv(kgrid,Vcv),
                   gPPEc2k(Ec,kgrid), gPPEv2k(Ev,kgrid),gPPdkdEc(ECm,dkdEc),gPPdkdEv(EVm,dkdEv),
                   gPPWcc(kgrid, Wcc), gPPWvv(kgrid,Wvv),gPPWcv(kgrid,Wcv) {}
                  

};

/* c++ generic functor invocation */

class InfPeak_diag
{  
public:
    double Ek1, Ek4;
    std::complex<double> qD; int Flag;
        
    std::vector<double> EC;
    std::vector<double> EV;

    Spline_interp gPPEc;
    Spline_interp gPPEv;
    Spline_interp gPPEc2k;
    Spline_interp gPPEv2k;
    
	InfPeak_diag(const double ek1, const double ek4, const std::complex<double>& qd, const int flag,   
	        const std::vector<double>& Ec, const std::vector<double>& Ev, const std::vector<double>& kgrid):
		    Ek1(ek1), Ek4(ek4), qD(qd), Flag(flag), EC{Ec}, EV{Ev}, gPPEc(kgrid, Ec), gPPEv(kgrid, Ev), 
			gPPEc2k(Ec, kgrid), gPPEv2k(Ev, kgrid) {}
		
	double operator()(const double x){ // Peakfinder<-operator overloading
		
		double k2  = x;
        double k3, Ek2, Ek3;
        std::vector<double> E3;
		
		if(Flag==1||Flag==4) // Ek2,Ek3-> c, SB1
		{
        E3=EC;
		Ek2= gPPEc.interp(k2);
		Ek3= Ek1+Ek2-Ek4;
	    }
		else                 // Ek2,Ek3-> v, SB2
		{
        E3=EV;
		Ek2= gPPEv.interp(k2);
		Ek3= Ek1+Ek2-Ek4;
	    }
		
		if (Ek3<=E3.front()||Ek3>=E3.back())
		{
		                     // reject, negative number sucks
		return 1.0e50;		 	
		}
        
        if(Flag==1||Flag==4)
		k3  = gPPEc2k.interp(Ek3);
		else
		k3  = gPPEv2k.interp(Ek3);
		
		double q_d = std::abs(qD);
		double b = k2*k2 + q_d*q_d - k3*k3;
		double ac= k2*k2 * q_d*q_d;
		return (b*b - 4.0*ac);
        
	}	

};


#endif /* SBE_H */
