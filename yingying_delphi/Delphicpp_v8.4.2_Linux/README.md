To all users:
==============

Current DelPhi C++ (Version 8.0) allows users to compile the DelPhi program 
into different versions by turning on and off appropriate flags from a single 
distribution. After downloading the DelPhi source code from DelPhi website, 
the following software must be pre-installed before the DelPhi program
is compiled:

1. appropriate C++ compiler such as GCC 5.4.0 and above;
2. boost library installed in /usr/include and its path is recognized in the user 
   environment;
3. latest version of OpenMPI or MPICH if MPI version of DelPhi executable is desired.

To compile the DelPhi C++ into:

a. a regular executable (simple single thread/CPU version):

1. go to the folder Release and type "make";
2. an executable named delphicpp_release will be generated when the compilation process 
   is finished.

b. a multi-threading OpenMP executable (able to utilize the maximum computing power of 
   a multi-core CPU to accelerate the calculations):
   
1. go to the source code folder src and uncomment the line "//#define PARALLEL_OMP" 
   in the file of src/interface/environment.h;
2. go to the compilation folder Release_omp and type "make";
3. an executable named delphicpp_omp_release will be generated when the compilation 
   process is finished.

c. a multi-CPU MPI executable (able to utilize the computing power of CPUs across 
   multiple computing nodes on one HPC cluster):
   
1. Go to the folder Release_mpi and type "make".
2. an executable named delphicpp_mpi_release will be generated when the compilation 
   process is finished.

To developers:
==============

Compilation folders (Debug, Debug_omp, and Debug_mpi) are also provided to developers so
that DelPhi can be compiled into corresponding "Debug" versions, which can be loaded into 
a debugger for finding the bugs, and studying how the program is executed, etc.. 

Regular users are NOT suggested to compile the DelPhi program in these folders!


Running the MPI version of the DelPhi program on a HPC cluster:
===============================================================

Users who are interested in running the MPI version of DelPhi program are advised to contact 
your administrator first before compiling and running the MPI version of the DelPhi program
on your HPC cluster. A sample PBS script to submit computing job on the Palmetto cluster 
(www.palmetto.clemson.edu) is provided in below for your easy reference but it is subject to 
necessary changes depending on your own computing environment.

#PBS -q workq
#PBS -l select=3:ncpus=1:mpiprocs=1:mem=120gb:interconnect=fdr
#PBS -l walltime=72:00:00

### ---------------------------------------
### BEGINNING OF EXECUTION
### ---------------------------------------

cd $PBS_O_WORKDIR
module purge
module add gcc/5.4.0 mpich/3.1.4

ln -sf ../../delphicpp_mpi_v4/Release_mpi/delphicpp_release ./delphicpp
ln -sf ../../test_cases/giant/4udf.pqr ./4udf.pqr
ln -sf ../../test_cases/work/a.pdb.charmm ./a.pdb.charmm
ln -sf ../../test_cases/work/charmm.crg ./charmm.crg
ln -sf ../../test_cases/work/charmm.siz ./charmm.siz

/bin/time -v mpiexec -n 3 ./delphicpp ./giant.para > 4udf_cpu3_linit_run3.txt

sleep 5
