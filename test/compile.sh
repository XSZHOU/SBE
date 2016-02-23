#!/bin/bash
## export for run time library linking
export MATLAB="/Applications/MATLAB_R2014b.app"
export DYLD_LIBRARY_PATH=$MATLAB/bin/maci64/:$MATLAB/sys/os/maci64/:DYLD_LIBRARY_PATH
echo `clear`
set -eu # makes program exit on error or unbound variable
## the same path specified /bin/maci64/ for compiling
mpic++ -o Xcc.x -Wall -O3 -fopenmp -std=c++11 -fdiagnostics-color=auto -Wextra -pedantic -fno-common -fexceptions\
       matio.cpp hcubature.c test_mpi_mat.cpp \
       -I/Users/seanzhou/Documents/FLENS \
       -I/Applications/MATLAB_R2014b.app/extern/include/ \
       -L/Applications/MATLAB_R2014b.app/bin/maci64/ \
       -lgsl -lmx -lmat -lgomp
mpirun -n 2 ./Xcc.x sbe_input.mat
echo -end '\n'

