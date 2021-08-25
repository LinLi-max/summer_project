#!/bin/bash
#SBATCH -n 16         
#SBATCH -t 02:00:00   
#SBATCH -p compute      
#SBATCH -J blocking_shock_wave 

# load the correct modules
module load gcc openmpi scorep papi

# set the SCOREP experiment variables
export SCOREP_EXPERIMENT_DIRECTORY=$HOME/mpi_blocking/mpi_blocking_profiling
export SCOREP_ENABLE_PROFILING=true
export SCOREP_OVERWRITE_EXPERIMENT_DIRECTORY=true

# for tracing
export SCOREP_ENABLE_TRACING=true
export SCOREP_METRIC_PAPI=PAPI_L2_DCM

# launch the code
make
mpirun -n 16 ./shock_wave
make clean
