#include <flens/flens.cxx>
#include <gsl/gsl_integration.h>
#include <algorithm>
#include <functional>
#include <numeric>  
#include <mpi.h>
#include <omp.h>
#include "sbe.h"
#include "cubature.h"
#include "scattering.h"
#include "matio.h"

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

#define WORKTAG 1
#define ExitTAG 2

using namespace flens;

/* Local functions */
static void master(const int, const int);
static void slave (const int, char *argv[]);
static void slave_compute(char *argv[], std::vector<double>& , const int);

/* constructor preparing*/
void Prep_CC_ctor( std::vector<double>& EC,   std::vector<double>& EV, 
                   std::vector<double>& kgrid,std::vector<double>& eps,
                   std::vector<double>& Vcc,  std::vector<double>& Vvv, std::vector<double>& Vcv);


/* main fun */
int main(int argc, char * argv[]) 
{
	
  int myrank; int msize=300;
  double t1, t2, diff;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  /* MPI_Barrier(MPI_COMM_WORLD);  */

  if (myrank == 0) 
  {
    t1=MPI_Wtime();
    master(myrank,msize);
    t2=MPI_Wtime();
    diff=t2-t1;
    //printf("Toltal execution time: %3.f\n", diff);
    std::cout<<"Total execution time: "<< diff << std::endl;
  }
  else 
    slave(myrank, argv);

  /* MPI_Barrier(MPI_COMM_WORLD); */
  MPI_Finalize();
  
  std::cout <<"process #"<<myrank<<" finished !"<<std::endl;
  
  return 0;
  
}


static void master(const int myrank, const int msize)
{	
	
  int ntasks, rank; double work;
  
  // FLENs
  typedef FullStorage<double, ColMajor>              FS;
  typedef GeMatrix<FS>                               GEMatrix;
  typedef DenseVector<Array<double> >                DenseVec;
  typedef DenseVector<ArrayView<double> >            DenseVecView;
  typedef DenseVec::IndexType                        IndexType;
  const Underscore<IndexType> _;
  
  //allocate
  GEMatrix XHF(msize,msize);
  GEMatrix XCC(msize,msize);
  
  //initialize working queue  
  std::vector<int> mywork(msize+1);
  for( int i = 0; i <= msize; i++ )         //  msize
  mywork.push_back( i );
     
  //communication vector 
  std::vector<double> result(2*msize+1);    //  2
  std::vector<double> workqueue(mywork.begin(), mywork.end());
  
  //MPI
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  for (rank = 1; rank < ntasks; ++rank) {
    
    work = workqueue.back();
    workqueue.pop_back();

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* data item is a double */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
	         
  }
  
  while ( workqueue.back()> 1e-3 ) {
    
    MPI_Recv(&result.front(),   /* message buffer */
             result.size(),     /* vector item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */
             
             std::cout<<"job id is "<< result.back()<<std::endl;
             
			 work = workqueue.back();
             workqueue.pop_back();
             
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* data item is a double */
             status.MPI_SOURCE, /* to who we just received from */  /*<---*/
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
             
              /***  write to Matrix  ***/
			 int ik1= static_cast<int>( result.back() ) ; 
			 DenseVecView veck = XHF(ik1, _);
			 for(int j=0; j< msize; ++j)
             veck(j+1)=result[j];
			 DenseVecView vecm = XCC(ik1, _);
			 for(int j=0; j< msize; ++j)
			 vecm(j+1)=result[j+msize];		 		 
			 /***  write to Matrix  ***/
             
  }


  for (rank = 1; rank < ntasks; ++rank) {
  	
    MPI_Recv(&result.front(),   /* message buffer */
             result.size(),     /* vector item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */
             
             std::cout<<"job id is "<< result.back()<<std::endl;
    
          	 /***  write to Matrix  ***/
			 int ik1= static_cast<int>( result.back() ) ; 
			 DenseVecView veck = XHF(ik1, _);
			 for(int j=0; j< msize; ++j)
             veck(j+1)=result[j];		 
             DenseVecView vecm = XCC(ik1, _);
			 for(int j=0; j< msize; ++j)
			 vecm(j+1)=result[j+msize];		 
			 /***  write to Matrix  ***/
  }


  for (rank = 1; rank < ntasks; ++rank)
    MPI_Send(0, 0, MPI_INT, rank, ExitTAG, MPI_COMM_WORLD);
  
  std::cout<<"master number# "<<myrank<<" gonna write output"<<std::endl;
  
  h5writeflens(XHF,"Xhf.h5","/XHF");
  
  h5writeflens(XCC,"Xcc.h5","/XCC");
  
}


static void slave(const int myrank, char *argv[])
{
     std::cout<<"slave number # "<<myrank<<" is on-line"<<std::endl;
   
  double work;
  MPI_Status status;
  
  while (1) 
  {
    
    MPI_Recv(&work, 1, MPI_DOUBLE, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == ExitTAG) {
        
      std::cout<<"slave number # "<<myrank<<" gonna be off-line"<<std::endl;
    
      return;
    }
    
    /*communication vector*/
    std::vector<double> result;
    
    /*do the actual job*/
    slave_compute(argv, result, static_cast<int>(work) );
    
    /*add the job id */
	result.push_back(work);
	
    /* Send the result back */
    MPI_Send(&result.front(), result.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
  }
    
}

/***********************************main workload*************************************************/

static void slave_compute(char *argv[], std::vector<double>& result, const int ik1)
{
	// read input  
    std::vector<double> F14, F23, SB1, SB2, XYv, EPS, DZ, EFC, EFV, kgrid;
    h5readflens(argv[1], F14, F23, SB1, SB2, XYv, EPS, DZ, EFC, EFV, kgrid);
    
    double dz= DZ[0]; double EFc= EFC[0];  double EFv= EFV[0];
    
	int nk= EPS.size(); int nkk= EPS.size()/2; int np= F14.size(); 
	double dk= kgrid[3]-kgrid[2]; result.resize(nkk*2);  // 2
    std::vector<double> stdEc(SB1);
    std::vector<double> stdEv(SB2);
    
    // FLENs
    typedef FullStorage<double, ColMajor>            FS;
    typedef GeMatrix<FS>                             GEMatrix;
    typedef FullStorageView<double, ColMajor>        FSView;
    typedef GeMatrix<FSView>                         GEMatrixView;
    typedef DenseVector<Array<double>>               DenseVec;
    typedef DenseVector<ArrayView<double> >          DenseVecView;
    //typedef Array<double> 				           DArray;
    typedef ArrayView<double> 			             DArrayView;
    const   Underscore<GEMatrix::IndexType> _;
    
    const double E0=8.854e-12, Q=1.6021917e-19, kB = 1.38e-23, T=300;
    const double ER_STA = 10.38;//HBAR=6.63e-34/(2.0*M_PI);

    DenseVecView f14 =  DArrayView(F14.size(),F14.data()); 
    DenseVecView f23 =  DArrayView(F23.size(),F23.data());    
    GEMatrixView  XY =  FSView(np, np, XYv.data(), np);

	std::vector<double> Vcc,Vvv,Vcv;
	
      for(int iq=0; iq<nk; ++iq)
       {

        double  qD, Vps;
        GEMatrix::IndexVariable ii, jj;
        DenseVec Vtmp(np); GEMatrix XYe(np, np);
       
        qD =kgrid[iq] ;
        XYe(ii,jj)= Exp(-(qD)*XY(ii,jj));
     
        Vtmp = XYe * f14; Vps = f14 * Vtmp;
        Vcc.push_back( sqr(Q)/(2*E0*ER_STA*qD)*sqr(dz)*Vps );
     
        Vtmp = XYe * f23; Vps = f23 * Vtmp;
        Vvv.push_back( sqr(Q)/(2*E0*ER_STA*qD)*sqr(dz)*Vps );
        
        Vtmp = XYe * f23; Vps = f14 * Vtmp;
        Vcv.push_back( sqr(Q)/(2*E0*ER_STA*qD)*sqr(dz)*Vps );
     
	   }
	    
	   std::vector<double> ECm(stdEc);std::vector<double> EVm(stdEv);
	   std::vector<double> dkdEc(kgrid);std::vector<double> dkdEv(EPS);
	   std::vector<double> Wcc(Vcc);std::vector<double> Wvv(Vvv);std::vector<double> Wcv(Vcv);
	    
	   Prep_CC_ctor(ECm,EVm,dkdEc,dkdEv,Wcc,Wvv,Wcv);//prepare for CC constructor
	   CInfPass_SBECC cubpass{kgrid,stdEc,stdEv,ECm,EVm,dkdEc,dkdEv,Vcc,Vvv,Vcv,Wcc,Wvv,Wcv};
	   cubpass.arr[1]=EFc; cubpass.arr[2]=EFv;
	   
       double k1norm= kgrid[ik1-1];
       double Ek1c =  stdEc[ik1-1];
       double rhok1c= 1.0/(1.0+std::exp( (Ek1c-EFc)/(kB*T/Q) ) );
       double Ek1v =  stdEv[ik1-1];
       double rhok1v= 1.0/(1.0+std::exp( (Ek1v-EFv)/(kB*T/Q) ) );
     
       cubpass.arr[0]=k1norm;
       cubpass.arr[3]=Ek1c;
       cubpass.arr[4]=Ek1v;
       cubpass.arr[5]=rhok1c;
       cubpass.arr[6]=rhok1v;
                    
   
        std::cout<<"->computing HF ..."<<std::endl;
    /******************************************************************************************HF*/
        #pragma omp parallel for firstprivate(cubpass)
		for (int ik2=0; ik2<nkk; ++ik2)
          {
           	
           	//std::cout<<omp_get_num_threads()<<std::endl;
           	//std::cout<<omp_get_thread_num()<<std::endl;
           	
            double k2norm= kgrid[ik2]; 

            double valnHF=0.0, err; 
            
            double xmin[2] = {std::max(k2norm-dk/2.0,0.0), 0}; 
			double xmax[2] = {k2norm+dk/2.0, 2.0*M_PI};
			cubature(1, IntegrandHF2d_ndiag, &cubpass, 2, xmin, xmax, \
			         1e4, 0.0, 1e-4, ERROR_INDIVIDUAL, &valnHF, &err);
			 
            /*** HF non-diag ***/
            result[ik2] = -(rhok1c+rhok1v-1)*valnHF;
            
          }
            /*** HF diagonal ***/
            double valdHF=0.0, err;
            double xmin[2] = {1e6, 0.0}; double xmax[2] = {3.0e9, 2*M_PI};
			cubature(1, IntegrandHF2d_diag, &cubpass, 2, xmin, xmax, \
			         1e4, 0.0, 1e-4, ERROR_INDIVIDUAL, &valdHF, &err);
            result[ik1-1] = valdHF;
    /******************************************************************************************HF*/
    
    
            
        std::cout<<"->computing CC ..."<<std::endl; omp_set_num_threads(8);
    /******************************************************************************************CC*/
        #pragma omp parallel for firstprivate(cubpass)
	for (int ik2=0; ik2<nkk; ++ik2)
          {
           	
            double k2norm= kgrid[ik2];  
            
                 if( (ik2+1) != ik1)
                 {
                    double valnCC=0.0;
                    for (int i=1 ;i<=6;i++)
                    { 
                    double res, abserr; size_t neval;
		    cubpass.arr[8]= k2norm; cubpass.arr[7]= i;	
					    
                    gsl_integration_cquad_workspace *ws = NULL ;
                    
		    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) 
	            {printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
		             abort(); }
	                
	            gsl_function F;  
                    F.params =  &cubpass; 
                    F.function =&IntegrandCC_ndiag;
                    
                    gsl_integration_cquad(&F, 0, 2*M_PI,0, 5e-3, ws, &res, &abserr, &neval);
	            valnCC += 0.5*dk*res;         /*PI/HBAR->0.5, trapz->dk*/
                    
                    gsl_integration_cquad_workspace_free( ws);
                    }
                 
		    result[ik2+nkk] = valnCC;     /*** CC non-diag ***/
		    std::cout<<omp_get_thread_num()<<" ik2: "<<ik2<<"ik1: "<<ik1<<std::endl;
                 }
                 
          } 
           
                       double valdCC=0.0;
	               for (int i=1;i<=4;i++)
			{
			double val=0.0,err; cubpass.arr[7]= i;
			cubature(1, IntegrandCC2d_diag, &cubpass, 2, xmin, xmax, \
			         1e4, 0.0, 3e-3, ERROR_INDIVIDUAL, &val, &err);
			/* Vcv/Vvc spin for secondary processes, without exchange cancelation */ 
                        if(i==3||i==4) val=val*2.0;  
			valdCC+=val;        
	         	}
			
                   result[ik1-1+nkk] = -valdCC*0.5;     /*** CC diagonal ***/
    /******************************************************************************************CC*/
			      
	
}


void Prep_CC_ctor( std::vector<double>& EC, std::vector<double>& EV, 
                   std::vector<double>& kgrid,std::vector<double>&eps,
                   std::vector<double>& Vcc, std::vector<double>& Vvv, std::vector<double>& Vcv)
  {				                 
   //shift, add, average, 0.5* ( SB1(ic,1:end-1)+SB2(iv,2:end) )
   std::vector<double> ECm, EVm;
   ECm.resize(EC.size()-1);
   std::transform(EC.begin(), EC.end()-1, EC.begin()+1, ECm.begin(), std::plus<double>());
   std::transform(ECm.begin(), ECm.end(), ECm.begin(),std::bind1st(std::multiplies<double>(),0.5));
   EVm.resize(EV.size()-1);
   std::transform(EV.begin(), EV.end()-1, EV.begin()+1, EVm.begin(), std::plus<double>());
   std::transform(EVm.begin(), EVm.end(), EVm.begin(),std::bind1st(std::multiplies<double>(),0.5));
   
   //dkdE diff(kgrid)/diff( SB(i,:) )
   std::vector<double>   dEc(EC.size());
   std::vector<double>   dEv(EV.size());
   std::vector<double> dk(kgrid.size());
   std::adjacent_difference (EC.begin(), EC.end(), dEc.begin());           /*[Ec(1) diff(Ec)]*/
   std::adjacent_difference (EV.begin(), EV.end(), dEv.begin());
   std::adjacent_difference (kgrid.begin(), kgrid.end(), dk.begin());
   dEc.erase (dEc.begin()); dEv.erase (dEv.begin()); dk.erase (dk.begin());/*delete first ele*/
   std::vector<double> dkdEc, dkdEv;
   dkdEc.resize(dk.size());
   dkdEv.resize(dk.size());
   std::transform(dk.begin(), dk.end(), dEc.begin(), dkdEc.begin(), std::divides<double>());
   std::transform(dk.begin(), dk.end(), dEv.begin(), dkdEv.begin(), std::divides<double>());
   
   // W=V/eps
   std::vector<double> Wcc, Wvv, Wcv;
   Wcc.resize(Vcc.size());
   Wvv.resize(Vvv.size());
   Wcv.resize(Vcv.size());
   std::transform(Vcc.begin(), Vcc.end(), eps.begin(), Wcc.begin(), std::divides<double>());
   std::transform(Vvv.begin(), Vvv.end(), eps.begin(), Wvv.begin(), std::divides<double>());
   std::transform(Vcv.begin(), Vcv.end(), eps.begin(), Wcv.begin(), std::divides<double>());
   
   // return by reference
   EC=ECm;EV=EVm;kgrid=dkdEc;eps=dkdEv;Vcc=Wcc;Vvv=Wvv;Vcv=Wcv;
                  
  }  
