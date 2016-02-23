#!/bin/bash

#-------------------- hw section -------------------------------

#### Ask for a number of chunks (select=) each with cpus "ncpus"
#PBS -l select=32:ncpus=8:mpiprocs=1:mem=1gb

#-------------------- queue system -------------------------------

#### Queue of submission
#PBS -q route

#### Maximum length of the job (hh:mm:ss)
#PBS -l walltime=00:30:00

# --------------- accounting/budget  -----------------------------

#### account number (type saldo -b)
#PBS -A try15_xianzhou 

# ---------------- other info ------------------------------------

#### File for standard output and error
#PBS -o job.out
#PBS -e job.err

#### Job name
#PBS -N Xcc

#### send email to the following address
##PBS -M 

#### send emai after abort or end
##PBS -m 

# -----------end of PBS keywords section -------------------------

cd $PBS_O_WORKDIR
#module load profile/advanced
module load autoload gsl/1.16--gnu--4.9.2
module load mkl/11.2--binary 
module load openmpi/1.8.4--gnu--4.9.2
module load autoload hdf5/1.8.14_ser--gnu--4.9.2
cat $PBS_NODEFILE
##export OMP_PROC_BIND= true
##mpirun -n 3 --report-bindings ./Xcc.x sbe_input.mat 
mpirun -n 32 --bind-to none ./Xcc.x sbe_input.mat 
