#include <gsl/gsl_integration.h>
#include "scattering.h"
#include "interp_1d.h"
#include "zbrent.h"
#include "sbe.h"

//using namespace std;

int IntegrandHF2d_ndiag(unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval)
{  // HF non-diagonal part, use 2d cubature instead of discrete trapz 

	const std::complex<double> I(0,1);
    std::complex<double> K2;
    double k2,qD;
    
	CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(fdata);
    
    double K1=cubpass->arr[0];

    k2= x[0];
    K2= k2*exp(I*x[1]);
    qD= std::abs(K1-K2);

    fval[0]= k2* cubpass->gPPVcv.interp(qD);
    return 0;
}

int IntegrandHF2d_diag(unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval)
{   // HF diagonal part, corresponding to Screened exchange shift
    
	const std::complex<double> I(0,1);
    std::complex<double> K2;
    double k2,qD,Q=1.6021917e-19, kB = 1.38e-23, T=300;
    
	CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(fdata);
    
    double K1=cubpass ->arr[0];
    double EFc=cubpass->arr[1];
    double EFv=cubpass->arr[2];

    qD= x[0];
    K2= K1 + qD*exp(I*x[1]);
    k2= std::abs(K2);
    qD= std::abs(K1-K2);

    double Ek2c = cubpass->gPPEc.interp(k2);
    double Ek2v = cubpass->gPPEv.interp(k2);
    double rhok2c=1.0/(1.0 + std::exp( (Ek2c-EFc)/(kB*T/Q) ));
    double rhok2v=1.0/(1.0 + std::exp( (Ek2v-EFv)/(kB*T/Q) ));
    
    fval[0]= qD* (rhok2c*cubpass->gPPVcc.interp(qD)+rhok2v*cubpass->gPPVvv.interp(qD));
    return 0;
}

/* ------------------------------------------------------------------------------------------- */ 
                    
double IntegrandCC_ndiag  ( const double x,  void *userdata)
{  // for non-diagonal part, we eliminate one integration over k using delta(E)
    const std::complex<double> I(0,1);
	CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    cubpass->arr[9]=x;   //theta_ka
    int flag =(int)cubpass->arr[7];
    
    
    double kmin=0.0; double kmax=3e9;
    double res, abserr, fval=0.0; size_t neval;
    
    gsl_integration_cquad_workspace *ws = NULL ;
    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) 
	{
      printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
      abort();
    }
    
	gsl_function F;  
    F.params =   cubpass;
    
	switch (flag) 
	    {
	     case 1:
	       F.function = &Integrand1d_cquadn1;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     case 2:
	       F.function = &Integrand1d_cquadn2;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     case 3:
	       F.function = &Integrand1d_cquadn3;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     case 4:
	       F.function = &Integrand1d_cquadn4;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     case 5:
	       F.function = &Integrand1d_cquadn5;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     case 6:
	       F.function = &Integrand1d_cquadn6;
	       gsl_integration_cquad(&F, kmin,kmax,0, 5e-3, ws, &res, &abserr, &neval);
	       fval=res;
	     break;
	     
	     default:
	     	printf( "type illegal !" );
	     	
	   }
    
    gsl_integration_cquad_workspace_free( ws);
    return fval;
}

int IntegrandCC2d_diag(unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval)
{  // for diagonal part, we eliminate one integration over theta using delta(E)
	
	const std::complex<double> I(0,1);
    std::complex<double> K4, QD;
    double k4,qD,Ek1,Ek4;
    
	CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(fdata);

    double K1  =cubpass->arr[0];
	int flag =(int)cubpass->arr[7];
	
    cubpass->arr[8]=x[0];//qD
    cubpass->arr[9]=x[1];//theta
    
    qD= x[0];
    QD= qD*exp(I*x[1]);
    K4 = K1 - QD;
    k4= std::abs(K4);
    
    if(flag==1||flag==3)
    {   // Ek1,Ek4->c
    Ek1= cubpass->gPPEc.interp(K1);
    Ek4= cubpass->gPPEc.interp(k4);
    }
    else
    {   // Ek1,Ek4->v
    Ek1= cubpass->gPPEv.interp(K1);
    Ek4= cubpass->gPPEv.interp(k4);	
    }
    
    double k2min=0.0; double k2max=3e9;
    double res, abserr; size_t neval;
	    
    gsl_integration_cquad_workspace *ws = NULL ;
    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) 
	{
      printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
      abort();
    }
    
	gsl_function F;  
    F.params =   cubpass;
    
    switch (flag) 
    {
     case 1:
       F.function = &Integrand1d_cquad1;
       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
       fval[0]=res;
     break;
     
     case 2:
       F.function = &Integrand1d_cquad2;
       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
       fval[0]=res;
     break;
     
     case 3:
        F.function = &Integrand1d_cquad3;
        gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
        fval[0]=res;
     break;
     
     case 4:
       F.function = &Integrand1d_cquad4;
       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
       fval[0]=res;
     break;
     
     default: //debugging by plotting the 1d-integrand
     	{
     	  int nk2=10001; double xx, yy;
		  fval[0]=0.0;
          for (int ii=0; ii<nk2; ii++)
          {
            xx=3e9*ii/(nk2-1.0);
            yy=Integrand1d_cquad3(xx,cubpass);
            fval[0]+=yy*3e9/(nk2-1.0);
          }
            printf("->check ,%e \n",fval[0]);
     	}  	
     	
   }
                  
    
    if(std::isinf(fval[0])||std::isnan(fval[0])) // in case return NaN/Inf
    { 
    
      printf( "->divergence experienced,call InfPeak !\n" );
      
	    try
	      {   
	       InfPeak_diag Peak_ftor{Ek1,Ek4,QD,flag,cubpass->EC,cubpass->EV,cubpass->KGRID};
	       k2min=zbrent(Peak_ftor, k2min, k2max, 1e-1);
	      }
	      //catch(NRerror s) 
	      catch(const char* p)                   // by S.Z.
	      {
	       //NRcatch(s);
	       printf("%s\n", p);
	       fval[0]=0.0; return 0;                // no solution found by Peakfinder
	      }  
    
         gsl_function F;
         F.params =   cubpass;
       
	     switch (flag) {
	     case 1:
	       F.function = &Integrand1d_cquad1;
	       gsl_integration_cquad(&F, 0,k2min,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]=res;
	       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]+=res;
	     break;
	     
	     case 2:
	       F.function = &Integrand1d_cquad2;
	       gsl_integration_cquad(&F, 0,k2min,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]=res;
	       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]+=res;
	     break;
	     
	     case 3:
	       F.function = &Integrand1d_cquad3;
	       gsl_integration_cquad(&F, 0,k2min,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]=res;
	       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]+=res;
	     break;
	     
	     case 4:
	       F.function = &Integrand1d_cquad4;
	       gsl_integration_cquad(&F, 0,k2min,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]=res;
	       gsl_integration_cquad(&F, k2min,k2max,0, 1e-3, ws, &res, &abserr, &neval);
	       fval[0]+=res;
	     break;
	     
	     default:
	     	printf( "type illegal !" );
	                  }

    }
    
    gsl_integration_cquad_workspace_free( ws);
    
    return 0;	
	
}

/* Integrand1d_cquad1-4, called by IntegrandCC_diag, 1-cccc-2-vvvv-3-cvvc-4-vccv*/
double Integrand1d_cquad1 ( const double x , void *userdata )
{   // Vcccc,1234
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);
    
    std::complex<double> K4, QD;
    double k2,k4,Jacob,invJacob,fcc;
    double Ek2c,Ek3c,Ek4c,rhok2c,rhok3c,rhok4c;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    double K1    = cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double Ek1c  = cubpass->arr[3];
    double qD    = cubpass->arr[8];
    double theta = cubpass->arr[9];
    
    QD= qD*std::exp(I*theta);
    K4 = K1 - QD;
    k4= std::abs(K4);
    k2=x; double fval;
    
    Ek2c = cubpass->gPPEc.interp(k2);
    Ek4c = cubpass->gPPEc.interp(k4);
    rhok2c= 1.0/(1.0+std::exp( (Ek2c-EFc)/(kB*T/Q) ) );
    rhok4c= 1.0/(1.0+std::exp( (Ek4c-EFc)/(kB*T/Q) ) );
    Ek3c= Ek1c+ Ek2c -Ek4c;
    
    if (Ek3c<=cubpass->EC.front()||Ek3c>=cubpass->EC.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEc2k.interp(Ek3c);
    rhok3c= 1.0/(1.0+std::exp( (Ek3c-EFc)/(kB*T/Q) ) );
    double dkdE = cubpass->gPPdkdEc.interp(Ek3c)/Q;
    double WqDcc= cubpass->gPPWcc.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {
    fval =0; return fval;
    }
    
    Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = (1.0-rhok2c)*rhok3c*rhok4c+ rhok2c *(1.0-rhok3c)*(1.0-rhok4c);
    
    fval=(WqDcc*WqDcc) * fcc * k2 * qD * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquad2 ( const double x , void *userdata )
{   // Vvvvv,1234
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);
    
    std::complex<double> K4, QD;
    double k2,k4,Jacob,invJacob,fcc;
    double Ek2v,Ek3v,Ek4v,rhok2v,rhok3v,rhok4v;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    double K1    = cubpass->arr[0];
    double EFv   = cubpass->arr[2];
    double Ek1v  = cubpass->arr[4];
    double qD    = cubpass->arr[8];
    double theta = cubpass->arr[9];
    
    QD= qD*std::exp(I*theta);
    K4 = K1 - QD;
    k4= std::abs(K4);
    k2=x; double fval;
    
    Ek2v = cubpass->gPPEv.interp(k2);
    Ek4v = cubpass->gPPEv.interp(k4);
    rhok2v= 1.0/(1.0+std::exp( (Ek2v-EFv)/(kB*T/Q) ) );
    rhok4v= 1.0/(1.0+std::exp( (Ek4v-EFv)/(kB*T/Q) ) );
    Ek3v= Ek1v+ Ek2v -Ek4v; 
    
    if (Ek3v<=cubpass->EV.front()||Ek3v>=cubpass->EV.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEv2k.interp(Ek3v);
    rhok3v= 1.0/(1.0+std::exp( (Ek3v-EFv)/(kB*T/Q) ) );
    double dkdE = cubpass->gPPdkdEv.interp(Ek3v)/Q;
    double WqDvv=cubpass->gPPWvv.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {
    fval =0; return fval;
    }
    
    Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = (1.0-rhok2v)*rhok3v*rhok4v+ rhok2v *(1.0-rhok3v)*(1.0-rhok4v);
    fval=(WqDvv*WqDvv) * fcc * k2 * qD * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquad3 ( const double x , void *userdata )
{   // Vcvcv,1234
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);
    
    std::complex<double> K4, QD;
    double k2,k4,Jacob,invJacob,fcc;
    double Ek2v,Ek3v,Ek4c,rhok2v,rhok3v,rhok4c;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    double K1    = cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double EFv   = cubpass->arr[2];
    double Ek1c  = cubpass->arr[3];
    double qD    = cubpass->arr[8];
    double theta = cubpass->arr[9];
    
    QD= qD*std::exp(I*theta);
    K4 = K1 - QD;
    k4= std::abs(K4);
    k2=x; double fval;
    
    Ek2v = cubpass->gPPEv.interp(k2);
    Ek4c = cubpass->gPPEc.interp(k4);
    rhok2v= 1.0/(1.0+std::exp( (Ek2v-EFv)/(kB*T/Q) ) );
    rhok4c= 1.0/(1.0+std::exp( (Ek4c-EFc)/(kB*T/Q) ) );
    //Ek3v=  Ek1c+ Ek2v -Ek4c; electron->electron hole picture
	Ek3v = -(Ek1c -Ek2v - Ek4c); 
    
    if (Ek3v<=cubpass->EV.front()||Ek3v>=cubpass->EV.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEv2k.interp(Ek3v);
    rhok3v= 1.0/(1.0+std::exp( (Ek3v-EFv)/(kB*T/Q) ) );
    double dkdE = cubpass->gPPdkdEv.interp(Ek3v)/Q;  
    double WqDcv= cubpass->gPPWcv.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);

    if(std::abs(costt2)>1.0)
    {
    fval =0; return fval;
    }
     
    Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc =  rhok2v*(1.0-rhok3v)*rhok4c+ (1.0-rhok2v) *rhok3v*(1.0-rhok4c);
    
    fval=(WqDcv*WqDcv) * fcc * k2 * qD * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquad4 ( const double x , void *userdata )
{   // Vvcvc,1234 
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);
    
    std::complex<double> K4, QD;
    double k2,k4,Jacob,invJacob,fcc;
    double Ek2c,Ek3c,Ek4v,rhok2c,rhok3c,rhok4v;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    double K1    = cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double EFv   = cubpass->arr[2];
    double Ek1v  = cubpass->arr[4];
    double qD    = cubpass->arr[8];
    double theta = cubpass->arr[9];
    
    QD= qD*std::exp(I*theta);
    K4 = K1 - QD;
    k4= std::abs(K4);
    k2=x; double fval;
    
    Ek2c = cubpass->gPPEc.interp(k2);
    Ek4v = cubpass->gPPEv.interp(k4);
    rhok2c= 1.0/(1.0+std::exp( (Ek2c-EFc)/(kB*T/Q) ) );
    rhok4v= 1.0/(1.0+std::exp( (Ek4v-EFv)/(kB*T/Q) ) );
    //Ek3c= Ek1v+ Ek2c -Ek4v; electron->electron hole picture
    Ek3c = -Ek1v + Ek2c + Ek4v; 
    
    if (Ek3c<=cubpass->EC.front()||Ek3c>=cubpass->EC.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEc2k.interp(Ek3c);
    rhok3c= 1.0/(1.0+std::exp( (Ek3c-EFc)/(kB*T/Q) ) );
    double dkdE = cubpass->gPPdkdEc.interp(Ek3c)/Q;
    double WqDvc= cubpass->gPPWcv.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {
    fval =0; return fval;
    }
    
    Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = (1.0-rhok2c)*rhok3c*(1.0-rhok4v)+ rhok2c*(1.0-rhok3c)*rhok4v;
    fval=(WqDvc*WqDvc) * fcc * k2 * qD * 2.0 * invJacob;
    
	return fval;
}

/*Integrand1d_cquadn1-6, called by IntegrandCC_ndiag, k&ka known, kb&kc unknown*/
double Integrand1d_cquadn1( const double x , void *userdata )
{   //"case_11", corresponds to #2 in eq.4.66 , Wcv_qx*Wcv_qd
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka,PD,K3,K4;
    double k3=x; double fval,fcc;     //k3->kc
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1= cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double EFv   = cubpass->arr[2];
    double Ek1v  = cubpass->arr[4];
    double rhok1v = cubpass->arr[6];

    double ka= cubpass->arr[8];       //ka->k2
    Ka = ka*std::exp(I*cubpass->arr[9]);
    PD = K1+ Ka;
    double pd= std::abs(PD);         //PD=K1+K2=K3+K4
    
    double Ekac= cubpass->gPPEc.interp(ka);
    double Ek3c= cubpass->gPPEc.interp(k3);
    double Ek4v= -(-Ek1v+Ekac-Ek3c); //electron->electron hole picture
    
    if (Ek4v<=cubpass->EV.front()||Ek4v>=cubpass->EV.back())
    {fval=0.0; return fval;}
    
    double k4=cubpass->gPPEv2k.interp(Ek4v);
    
    double cospdk3= (k3*k3+pd*pd-k4*k4)/(2.0*k3*pd);
    double cospdk1= (K1*K1+pd*pd-ka*ka)/(2.0*K1*pd);
    
    if(std::abs(cospdk3)>1.0)      
    {fval=0.0; return fval;}
    
    // dk4/dE4* dtheta3/dk4, k1-k2-pd fixed-> dtheta3= dthetapdk3    
    double dkdE = cubpass->gPPdkdEv.interp(Ek4v)/Q;
    double Jacob= k3*pd*std::sqrt(1-cospdk3*cospdk3)/k4; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    // reconstruct K3,K4, acos->[0 pi]
    double thetapdk1 =std::acos ( cospdk1 ) * 180.0 / M_PI;
    double thetapdk3 =std::acos ( cospdk3 ) * 180.0 / M_PI;
    K3= k3*std::exp(I*(thetapdk1+thetapdk3));
    K4= PD-K3;
    
    // compute qD, qX
    double qD= std::abs(K1-K4);
    double qX= std::abs(K1-K3);
    if(qD>6e9||qX>6e9) {fval=0.0;return fval;}
    double WqDvc = cubpass->gPPWcv.interp(qD);
    double WqXvc = cubpass->gPPWcv.interp(qX);
    
    double rhok3c= 1.0/(1.0+std::exp( (Ek3c-EFc)/(kB*T/Q) ) );
    double rhok4v= 1.0/(1.0+std::exp( (Ek4v-EFv)/(kB*T/Q) ) );
    
    fcc = (1.0-rhok4v)*rhok3c*rhok1v+ rhok4v*(1.0-rhok3c)*(1.0-rhok1v);
    fval= (WqDvc*WqXvc) * fcc * ka * k3 * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquadn2( const double x , void *userdata )
{   //"case_21", corresponds to #2 in eq.4.66, Wcv_qx*Wcv_qd
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka,PD,K3,K4;
    double k3=x; double fval,fcc;   //k3v->kb
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1= cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double EFv   = cubpass->arr[2];
    double Ek1c  = cubpass->arr[3];
    double rhok1c= cubpass->arr[5];

    double ka= cubpass->arr[8];       //ka->k2
    Ka = ka*std::exp(I*cubpass->arr[9]);
    PD = K1+ Ka;
    double pd= std::abs(PD);         //PD=K1+K2=K3+K4
    
    double Ekav= cubpass->gPPEv.interp(ka);
    double Ek3v= cubpass->gPPEv.interp(k3);
    double Ek4c= Ek1c-Ekav+Ek3v;     //electron->electron hole picture
    
    if (Ek4c<=cubpass->EC.front()||Ek4c>=cubpass->EC.back())
    {fval=0.0; return fval;}
    
    double k4=cubpass->gPPEc2k.interp(Ek4c);
    
    double cospdk3= (k3*k3+pd*pd-k4*k4)/(2.0*k3*pd);
    double cospdk1= (K1*K1+pd*pd-ka*ka)/(2.0*K1*pd);
    
    if(std::abs(cospdk3)>1.0)      
    {fval=0.0; return fval;}
    
    // dk4/dE4* dtheta3/dk4, k1-k2-pd fixed-> dtheta3= dthetapdk3    
    double dkdE = cubpass->gPPdkdEc.interp(Ek4c)/Q;
    double Jacob= k3*pd*std::sqrt(1-cospdk3*cospdk3)/k4; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    // reconstruct K3,K4
    double thetapdk1 =std::acos ( cospdk1 ) * 180.0 / M_PI;
    double thetapdk3 =std::acos ( cospdk3 ) * 180.0 / M_PI;
    K3= k3*std::exp(I*(thetapdk1+thetapdk3));
    K4= PD-K3;
    
    // compute qD, qX
    double qD= std::abs(K1-K4);
    double qX= std::abs(K1-K3);
    if(qD>6e9||qX>6e9) {fval=0.0;return fval;}
    double WqDcv = cubpass->gPPWcv.interp(qD);
    double WqXcv = cubpass->gPPWcv.interp(qX);
    
    double rhok3v= 1.0/(1.0+std::exp( (Ek3v-EFv)/(kB*T/Q) ) );
    double rhok4c= 1.0/(1.0+std::exp( (Ek4c-EFc)/(kB*T/Q) ) );
    
    fcc = (1.0-rhok3v)*rhok4c*(1.0-rhok1c)+ rhok3v*(1.0-rhok4c)*rhok1c;
    fval= (WqDcv*WqXcv) * fcc * ka * k3 * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquadn3( const double x , void *userdata )
{   // -"case12"->3 + "case13"->5, Wcv_qd*(Wcc_qd-Wcc_qx)
    const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka, QD;
    double k2=x; double fval,fcc;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1= cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double Ek1v  = cubpass->arr[4];
    double rhok1v= cubpass->arr[6];

    double ka= cubpass->arr[8];       //ka->k4
    Ka = ka*std::exp(I*cubpass->arr[9]);
    QD = K1 - Ka;                     //q=k1-k4
    double qD= std::abs(QD);         
    
    double Ekav= cubpass->gPPEv.interp(ka);
    double Ek2c= cubpass->gPPEc.interp(k2);
    double Ek3c= -Ek1v+Ek2c+Ekav;     //electron->electron hole picture
    
    if (Ek3c<=cubpass->EC.front()||Ek3c>=cubpass->EC.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEc2k.interp(Ek3c);
    
    double rhok2c= 1.0/(1.0+std::exp( (Ek2c-EFc)/(kB*T/Q) ) );
    double rhok3c= 1.0/(1.0+std::exp( (Ek3c-EFc)/(kB*T/Q) ) );
    
	double dkdE = cubpass->gPPdkdEc.interp(Ek3c)/Q;
    double WqDvc= cubpass->gPPWcv.interp(qD);
    double WqDcc= cubpass->gPPWcc.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {fval =0; return fval;}
    
    double Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = rhok1v*(1.0-rhok2c)*rhok3c+ (1.0-rhok1v)*(1.0-rhok3c)*rhok2c;
    fval=(WqDcc*WqDvc) * fcc * k2 * ka * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquadn4( const double x , void *userdata )
{  //"case13"->6,7 , Wvc_qd*(Wvv_qd-Wvv_qx)
	const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka, QD;
    double k2=x; double fval,fcc;
   
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1    = cubpass->arr[0];
    double EFv   = cubpass->arr[2];
    double Ek1v  = cubpass->arr[4];
    double rhok1v= cubpass->arr[6];

    double ka= cubpass->arr[8];       //ka->k4
    Ka = ka*std::exp(I*cubpass->arr[9]);
    QD = K1 - Ka;                     //q=k1-k4
    double qD= std::abs(QD);         
    
    double Ekav= cubpass->gPPEv.interp(ka);
    double Ek2v= cubpass->gPPEv.interp(k2);
    double Ek3v= Ek1v+Ek2v-Ekav;     //electron->electron hole picture
    
    if (Ek3v<=cubpass->EV.front()||Ek3v>=cubpass->EV.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEv2k.interp(Ek3v);
    
    double rhok2v= 1.0/(1.0+std::exp( (Ek2v-EFv)/(kB*T/Q) ) );
    double rhok3v= 1.0/(1.0+std::exp( (Ek3v-EFv)/(kB*T/Q) ) );
    
	double dkdE = cubpass->gPPdkdEv.interp(Ek3v)/Q;
    double WqDcv= cubpass->gPPWcv.interp(qD);
    double WqDvv= cubpass->gPPWvv.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {fval =0; return fval;}
    
    double Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = rhok1v*rhok2v*(1.0-rhok3v)+ (1.0-rhok1v)*(1.0-rhok2v)*rhok3v;
    fval=(WqDcv*WqDvv) * fcc * k2 * ka * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquadn5( const double x , void *userdata )
{  //-"case22"->4 + "case23"->8 ,Wcv_qd*(Wvv_qd-Wvv_qx)
	const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka, QD;
    double k2=x; double fval,fcc;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1= cubpass->arr[0];
    double EFv   = cubpass->arr[2];
    double Ek1c  = cubpass->arr[3];
    double rhok1c= cubpass->arr[5];

    double ka= cubpass->arr[8];       //ka->k4
    Ka = ka*std::exp(I*cubpass->arr[9]);
    QD = K1 - Ka;                     //q=k1-k4
    double qD= std::abs(QD);         
    
    double Ekac= cubpass->gPPEc.interp(ka);
    double Ek2v= cubpass->gPPEv.interp(k2);
    double Ek3v=-(Ek1c-Ek2v-Ekac);     //electron->electron hole picture
    
    if (Ek3v<=cubpass->EV.front()||Ek3v>=cubpass->EV.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEv2k.interp(Ek3v);
    
    double rhok2v= 1.0/(1.0+std::exp( (Ek2v-EFv)/(kB*T/Q) ) );
    double rhok3v= 1.0/(1.0+std::exp( (Ek3v-EFv)/(kB*T/Q) ) );
    
	double dkdE = cubpass->gPPdkdEv.interp(Ek3v)/Q;
    double WqDcv= cubpass->gPPWcv.interp(qD);
    double WqDvv= cubpass->gPPWvv.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {fval =0; return fval;}
    
    double Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = (1.0-rhok1c)*rhok2v*(1.0-rhok3v)+ rhok1c*(1.0-rhok2v)*rhok3v;
    fval=(WqDvv*WqDcv) * fcc * k2 * ka * 2.0 * invJacob;
    
	return fval;
}

double Integrand1d_cquadn6( const double x , void *userdata )
{  //"case23"->9,10 Wvc_qd*(Wcc_qd-Wcc_qx)
	const double kB = 1.38e-23, T=300, Q = 1.6021917e-19;
    const std::complex<double> I(0,1);std::complex<double> Ka, QD;
    double k2=x; double fval,fcc;
    
    CInfPass_SBECC *cubpass=static_cast<CInfPass_SBECC*>(userdata);
    
    double K1    = cubpass->arr[0];
    double EFc   = cubpass->arr[1];
    double Ek1c  = cubpass->arr[3];
    double rhok1c= cubpass->arr[5];

    double ka= cubpass->arr[8];       //ka->k4
    Ka = ka*std::exp(I*cubpass->arr[9]);
    QD = K1 - Ka;                     //q=k1-k4
    double qD= std::abs(QD);         
    
    double Ekac= cubpass->gPPEc.interp(ka);
    double Ek2c= cubpass->gPPEc.interp(k2);
    double Ek3c= Ek1c+Ek2c-Ekac;     //electron->electron hole picture
    
    if (Ek3c<=cubpass->EC.front()||Ek3c>=cubpass->EC.back())
    {fval=0.0; return fval;}
    
    double k3=cubpass->gPPEc2k.interp(Ek3c);
    
    double rhok2c= 1.0/(1.0+std::exp( (Ek2c-EFc)/(kB*T/Q) ) );
    double rhok3c= 1.0/(1.0+std::exp( (Ek3c-EFc)/(kB*T/Q) ) );
    
	double dkdE = cubpass->gPPdkdEc.interp(Ek3c)/Q;
    double WqDvc= cubpass->gPPWcv.interp(qD);
    double WqDcc= cubpass->gPPWcc.interp(qD);
    
    double costt2= -(k2*k2+qD*qD-k3*k3)/(2.0*k2*qD);
    
    if(std::abs(costt2)>1.0)
    {fval =0; return fval;}
    
    double Jacob=k2*qD*std::sqrt(1-costt2*costt2)/k3; 
    double invJacob=std::abs(dkdE)/std::abs(Jacob);
    
    fcc = rhok1c*rhok3c*(1.0-rhok2c)+ (1.0-rhok1c)*(1.0-rhok3c)*rhok2c;
    fval=(WqDvc*WqDcc) * fcc * k2 * ka * 2.0 * invJacob;
    
	return fval;
}
