//plain-mpi test for multi-processes reading mat-file
#pragma GCC diagnostic ignored "-Wwrite-strings"

#include "mat.h"
#include <vector>
#include <iostream>
#include <mpi.h>

using namespace std;
                   
int matread(const char *file, std::vector<double>& p1, std::vector<double>& p2,
                              std::vector<double>& p3, std::vector<double>& p4,
			                  std::vector<double>& p5, std::vector<double>& p6,
                              std::vector<double>& p7, std::vector<double>& p8,
                              std::vector<double>& p9, std::vector<double>& p10);

void matwrite(const std::vector<double>& p1);

//main

int main(int argc, char * argv[]) {
    
   int rank, size;
   MPI_Init (&argc, &argv);      
   MPI_Comm_size (MPI_COMM_WORLD, &size);        
   MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
  
   std::vector<double> SB1, SB2, epsfun, XY, dz, F14, F23, efc, efv, kgrid;
    
   if (argc == 2)
    {
      if (std::string(argv[1]) == "sbe_input.mat")
      {
      matread(argv[1],SB1, SB2, F14, F23, XY, dz,efc, efv, epsfun, kgrid);
      cout<< "->read data finished"<<std::endl;
      }
      else
      {
      cout<< std::string(argv[1])<< " <- illegal file !"<<std::endl;
      return 0;
      }
    }
	else
	{
    printf("require input .mat file");
    return 0;
	}

    cout<< efc[0]<<endl;
    cout<< efv[0]<<endl;
    cout<< dz[0]<<endl;
    cout<< XY[0]<<"  "<<XY[1]<<"  "<<XY[2]<<endl;

    MPI_Finalize();
    
return 0;

}

/***********************************IO-MAT-API********************************************/
void matwrite(const std::vector<double>& p1) {
    
    MATFile *pmat;
    mxArray *pma1;
    const char *file = "Xcc.mat";
    
    pmat = matOpen(file, "w");
    if (pmat == NULL) {
        printf("Error creating file %s\n", file);}
    
    pma1 = mxCreateDoubleMatrix(p1.size(),1,mxREAL);
    
    //memcpy(mxGetPr(pma1), &p1[0], p1.size() * sizeof(double));
    std::copy(p1.begin(),p1.end(),mxGetPr(pma1));
    
    matPutVariable(pmat, "Xcc", pma1);
    
    matClose(pmat);
    mxDestroyArray(pma1);
    
    //printf("*---------output written to mat--------*\n");
    
}

int matread(const char *file, std::vector<double>& p1, std::vector<double>& p2,
     std::vector<double>& p3, std::vector<double>& p4, std::vector<double>& p5,
     std::vector<double>& p6, std::vector<double>& p7, std::vector<double>& p8,
     std::vector<double>& p9, std::vector<double>& p10) {
    
    //printf("*----------input read from mat---------*\n");
    
    MATFile *pmat;
    const char **dir;
    const char *name;
    int	  ndir;
    mxArray *pma1,*pma2,*pma3,*pma4,*pma5,*pma6,*pma7,*pma8,*pma9,*pma10;
    double  *pda1,*pda2,*pda3,*pda4,*pda5,*pda6,*pda7,*pda8,*pda9,*pda10;
    
    pmat=matOpen(file,"r");
    if (pmat == NULL) {
        printf("Error opening file %s\n", file);
        return(1);
    }
    
    /* checking fields */
    /*
    dir = (const char **)matGetDir(pmat, &ndir);
    if (dir == NULL) {
        printf("Error reading directory of file %s\n", file);
        return(1);
    } else {
        printf("Directory of %s:\n", file);
        for (int i=0; i < ndir; i++)
            printf("%s\n",dir[i]);
    }
    mxFree(dir);
    */
    
    /*close and reopen to read*/
    if(matClose(pmat)!=0) {
        printf("error closing\n ");return 0;}
    pmat=matOpen(file,"r");
    
    pma1 = matGetNextVariable(pmat, &name);
    pma2 = matGetNextVariable(pmat, &name);
    pma3 = matGetNextVariable(pmat, &name);
    pma4 = matGetNextVariable(pmat, &name);
    pma5 = matGetNextVariable(pmat, &name);
    pma6 = matGetNextVariable(pmat, &name);
    pma7 = matGetNextVariable(pmat, &name);
    pma8 = matGetNextVariable(pmat, &name);
    pma9 = matGetNextVariable(pmat, &name);
    pma10 = matGetNextVariable(pmat, &name);
    
    if(pma10==NULL){
        printf("error reading info\n ");return 0;}
    
    int rows1 = mxGetM(pma1); int cols1 = mxGetN(pma1);
    int rows2 = mxGetM(pma2); int cols2 = mxGetN(pma2);
    int rows3 = mxGetM(pma3); int cols3 = mxGetN(pma3);
    int rows4 = mxGetM(pma4); int cols4 = mxGetN(pma4);
    int rows5=  mxGetM(pma5); int cols5 = mxGetM(pma5);
    int rows6=  mxGetM(pma6); int cols6 = mxGetM(pma6);
    int rows7=  mxGetM(pma7); int cols7 = mxGetM(pma7);
    int rows8=  mxGetM(pma8); int cols8 = mxGetM(pma8);
    int rows9=  mxGetM(pma9); int cols9 = mxGetM(pma9);
    int rows10=  mxGetM(pma10); int cols10 = mxGetM(pma10);
    
    pda1 = mxGetPr(pma1);
    pda2 = mxGetPr(pma2);
    pda3 = mxGetPr(pma3);
    pda4 = mxGetPr(pma4);
    pda5 = mxGetPr(pma5);
    pda6 = mxGetPr(pma6);
    pda7 = mxGetPr(pma7);
    pda8 = mxGetPr(pma8);
    pda9 = mxGetPr(pma9);
    pda10 = mxGetPr(pma10);
        
    p1.resize(rows1*cols1);
    p2.resize(rows2*cols2);
    p3.resize(rows3*cols3);
    p4.resize(rows4*cols4);
    p5.resize(rows5*cols5);
    p6.resize(rows6*cols6);
    p7.resize(rows7*cols7);
    p8.resize(rows8*cols8);
    p9.resize(rows9*cols9);
    p10.resize(rows10*cols10);
    p1.assign(pda1, pda1 + rows1*cols1);
    p2.assign(pda2, pda2 + rows2*cols2);
    p3.assign(pda3, pda3 + rows3*cols3);
    p4.assign(pda4, pda4 + rows4*cols4);
    p5.assign(pda5, pda5 + rows5*cols5);
    p6.assign(pda6, pda6 + rows6*cols6);
    p7.assign(pda7, pda7 + rows7*cols7);
    p8.assign(pda8, pda8 + rows8*cols8);
    p9.assign(pda9, pda9 + rows9*cols9);
    p10.assign(pda10, pda10 + rows10*cols10);
    
  /*memcpy(&p1[0],pda1,rows*cols*sizeof(double));
    memcpy(&p2[0],pda2,rows*cols*sizeof(double));
    memcpy(&p3[0],pda3,rows*cols*sizeof(double));
    memcpy(&p4[0],pda4,rows*cols*sizeof(double));*/
  
    if(matClose(pmat)!=0) {
        printf("error closing\n ");return 0;}
    
    mxDestroyArray(pma1);
    mxDestroyArray(pma2);
    mxDestroyArray(pma3);
    mxDestroyArray(pma4);
    mxDestroyArray(pma5);
    mxDestroyArray(pma6);
    mxDestroyArray(pma7);
    mxDestroyArray(pma8);
    mxDestroyArray(pma9);
    mxDestroyArray(pma10);
    
    return 0;
    
}
/***********************************IO-MAT-API*******************************************/
