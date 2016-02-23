CC =gcc
CFLAG = -O3
CXX = mpiCC
CXXFLAG = -DWITH_MKLBLAS -DMKL_ILP64
CXXFLAG+= -std=c++11 -O3 -fopenmp -Wall -pedantic -fexceptions -fdiagnostics-color=always
CXXINC  = -I${GSL_INC} -I${HDF5_INC} -I${MKL_INC} -I./FLENS -I./include 
LDFLAGS = -lgomp  -L${GSL_LIB} -lgslcblas -lgsl -lpthread -lm
LDFLAGS+= -L${MKL_LIB} -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -L${HDF5_LIB} -lhdf5_cpp -lhdf5

all: Xcc.x 

Xcc.x: test_mpi_hdf5.o matio.o scattering.o interp_1d.o hcubature.o
	        $(CXX)  -o Xcc.x test_mpi_hdf5.o matio.o scattering.o interp_1d.o hcubature.o $(LDFLAGS)

test_mpi_hdf5.o: test_mpi_hdf5.cpp
	        $(CXX) $(CXXFLAG) $(CXXINC) -c test_mpi_hdf5.cpp

interp_1d.o: interp_1d.cpp
		$(CXX) $(CXXFLAG) $(CXXINC) -c interp_1d.cpp

matio.o: matio.cpp
	        $(CXX) $(CXXFLAG) $(CXXINC) -c matio.cpp

scattering.o: scattering.cpp
	        $(CXX) $(CXXFLAG) $(CXXINC) -c scattering.cpp

hcubature.o: hcubature.c
	        $(CXX) $(CXXFLAG) $(CXXINC) -c hcubature.c 

clean:
	        rm -f *.o *.h5 


