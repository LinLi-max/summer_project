#!/bin/bash

module load gcc openmpi

make

np=(1 2 4 6 8 10 12)

for(( i = 0; i < ${#np[@]} ; i ++ ))
do
    mpirun -n ${np[i]} ./shock_wave >> result.out
done

make clean