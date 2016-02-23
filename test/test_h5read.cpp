#include <iostream>
#include <vector>
#include <string>
#include "H5Cpp.h"

//g++ test_matio.cpp  -std=c++11 -fdiagnostics-color=always -lhdf5_cpp -lhdf5
//test HDF5-read API

int main()
{
    std::string ifn = "SBE_input.h5";
    H5::H5File fp(ifn.c_str(),H5F_ACC_RDONLY);
    
	std::string datasetPath = "/F14";
    H5::DataSet dset1 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace1 = dset1.getSpace();
    //H5T_class_t type_class = dset.getTypeClass();
    hsize_t rank; hsize_t dims[2];
    rank = dspace1.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace1(rank,dims);
    std::vector<double> p1;
    p1.resize(dims[0]*dims[1]);
    dset1.read(p1.data(), H5::PredType::NATIVE_DOUBLE, memspace1, dspace1);
    
    datasetPath = "/F23";
    H5::DataSet dset2 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace2 = dset2.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace2.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace2 (rank,dims);
    std::vector<double> p2;
    p2.resize(dims[0]*dims[1]);
    dset2.read(p2.data(), H5::PredType::NATIVE_DOUBLE, memspace2, dspace2);
    
    datasetPath = "/SB1";
    H5::DataSet dset3 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace3 = dset3.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace3.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace3 (rank,dims);
    std::vector<double> p3;
    p3.resize(dims[0]*dims[1]);
    dset3.read(p3.data(), H5::PredType::NATIVE_DOUBLE, memspace3, dspace3);
    
    datasetPath = "/SB2";
    H5::DataSet dset4 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace4 = dset4.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace4.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace4 (rank,dims);
    std::vector<double> p4;
    p4.resize(dims[0]*dims[1]);
    dset4.read(p4.data(), H5::PredType::NATIVE_DOUBLE, memspace4, dspace4);
    
    datasetPath = "/XY";
    H5::DataSet dset5 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace5 = dset5.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace5.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace5 (rank,dims);
    std::vector<double> p5;
    p5.resize(dims[0]*dims[1]);
    dset5.read(p5.data(), H5::PredType::NATIVE_DOUBLE, memspace5, dspace5);
    
    datasetPath = "/epsfun";
    H5::DataSet dset6 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace6 = dset6.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace6.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace6 (rank,dims);
    std::vector<double> p6;
    p6.resize(dims[0]*dims[1]);
    dset6.read(p6.data(), H5::PredType::NATIVE_DOUBLE, memspace6, dspace6);
    
    datasetPath = "/dz";
    H5::DataSet dset7 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace7 = dset7.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace7.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace7 (rank,dims);
    std::vector<double> p7;
    p7.resize(dims[0]*dims[1]);
    dset7.read(p7.data(), H5::PredType::NATIVE_DOUBLE, memspace7, dspace7);
    
    datasetPath = "/efc";
    H5::DataSet dset8 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace8 = dset8.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace8.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace8 (rank,dims);
    std::vector<double> p8;
    p8.resize(dims[0]*dims[1]);
    dset8.read(p8.data(), H5::PredType::NATIVE_DOUBLE, memspace8, dspace8);
    
    datasetPath = "/efv";
    H5::DataSet dset9 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace9 = dset9.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace9.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace9 (rank,dims);
    std::vector<double> p9;
    p9.resize(dims[0]*dims[1]);
    dset9.read(p9.data(), H5::PredType::NATIVE_DOUBLE, memspace9, dspace9);
    
    datasetPath = "/kgrid";
    H5::DataSet dset10 = fp.openDataSet(datasetPath.c_str());
    H5::DataSpace dspace10 = dset10.getSpace();
    //H5T_class_t type_class = dset1.getTypeClass();
    rank = dspace10.getSimpleExtentDims(dims, NULL);
    H5::DataSpace memspace10 (rank,dims);
    std::vector<double> p10;
    p10.resize(dims[0]*dims[1]);
    dset10.read(p10.data(), H5::PredType::NATIVE_DOUBLE, memspace10, dspace10);

    std::cout<< p7[0] << " " <<p8[0]<< " "<<p9[0] <<std::endl;
    fp.close();
    
    return 0;
    
}
