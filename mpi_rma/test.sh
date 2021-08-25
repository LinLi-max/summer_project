#!/bin/sh
#SBATCH -n 64         
#SBATCH -t 12:00:00   
#SBATCH -p compute      
#SBATCH -J rma_shock_wave  

module load gcc openmpi

make

np=(2 4 8 16 32 64)

for(( i = 0; i < ${#np[@]} ; i ++ ))
do
    mpirun -n ${np[i]} ./shock_wave >> rma.out
done

make clean