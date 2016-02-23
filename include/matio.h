#ifndef matio_H
#define matio_H

#include <flens/flens.cxx>
#include <vector>
#include <string>
#include "H5Cpp.h"

template <typename MA>
int h5writeflens(flens::GeMatrix<MA> &A,const std::string& filename, const std::string& datasetname)
{
    
    hsize_t nx = A.numCols();
    hsize_t ny = A.numRows();
    
    H5::H5File fp(filename.c_str(), H5F_ACC_TRUNC);
    
    hsize_t dims[2] = {nx, ny};
    H5::DataSpace dataspace(2, dims);  //2 is the rank
    
    H5::DataSet dataset = fp.createDataSet(datasetname.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
    
    dataset.write(A.data(), H5::PredType::NATIVE_DOUBLE);
    
    fp.close();
    
    return 0;
	
}

int h5readflens( const char *file,   std::vector<double>& p1, std::vector<double>& p2,
            std::vector<double>& p3, std::vector<double>& p4, std::vector<double>& p5, std::vector<double>& p6,
            std::vector<double>& p7, std::vector<double>& p8, std::vector<double>& p9, std::vector<double>& p10);            
            
            
#endif 
