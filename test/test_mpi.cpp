//test the master-slave working paradigm
#pragma GCC diagnostic ignored "-Wwrite-strings"

//#include "mat.h"
#include <vector>
#include <iostream>
#include <mpi.h>
#include "matread.h"

#define WORKTAG 1
#define ExitTAG 2

/* Local functions */
static void master(int);
static void slave(int, char *argv[]);

/* modulized into matread.c
int matread( const char *file,       std::vector<double>& p1, std::vector<double>& p2,
            std::vector<double>& p3, std::vector<double>& p4, std::vector<double>& p5, std::vector<double>& p6,
            std::vector<double>& p7, std::vector<double>& p8, std::vector<double>& p9, std::vector<double>& p10);
*/
       
int main(int argc, char * argv[]) {
	
  int myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  /* MPI_Barrier(MPI_COMM_WORLD);  */

  if (myrank == 0) 
    master(myrank);
  else 
    slave(myrank, argv);

  /* MPI_Barrier(MPI_COMM_WORLD); */
  MPI_Finalize();
  
  std::cout <<"process #"<<myrank<<" finished !"<<std::endl;
  
  return 0;
}


static void master(int myrank)
{
  int ntasks, rank, msize=300;
  double work;
  
  std::vector<int> mywork;
  for( int i = 0; i <= msize; i++ )
  mywork.push_back( i );
    
  std::vector<double> result(msize+1);  
  std::vector<double> workqueue(mywork.begin(), mywork.end());
  
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  for (rank = 1; rank < ntasks; ++rank) {

    work = workqueue.back();
    workqueue.pop_back();

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
  }

  
  while ( workqueue.back()> 1e-6 ) {  //changed by S.Z. here

    work = workqueue.back();
    workqueue.pop_back();
    
    MPI_Recv(&result.front(),   /* message buffer */
             result.size(),     /* vector item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */
             
             std::cout<<"job id is "<< result.back()<<std::endl;

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */  /*<---*/
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
             
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
  }


  for (rank = 1; rank < ntasks; ++rank)
    MPI_Send(0, 0, MPI_INT, rank, ExitTAG, MPI_COMM_WORLD);
  
  std::cout<<"master number# "<<myrank<<" gonna be off line"<<std::endl;
  
}


static void slave(int myrank, char *argv[])
{
    
  double work;
  MPI_Status status;
  
  while (1) 
  {
    
    MPI_Recv(&work, 1, MPI_DOUBLE, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == ExitTAG) {
        
      std::cout<<"slave number # "<<myrank<<" gonna be off line"<<std::endl;
    
      return;
    }
    
    /* Do the work */  
    std::vector<double> SB1, SB2, epsfun, XY, dz, F14, F23, efc, efv, kgrid;
    matread(argv[1],SB1, SB2, F14, F23, XY, dz,efc, efv, epsfun, kgrid);
    std::vector<double> result(kgrid.begin(), kgrid.begin()+300); 
	result.push_back(work);//std::cout <<result.size()<<std::endl;
    
    /* Send the result back */
    MPI_Send(&result.front(), result.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
  }
    
}
