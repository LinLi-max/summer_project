#!/bin/bash

module load cuda
make

TPB=(4 8 16 32 64 128 512 1024)

for(( k = 0; k < ${#TPB[@]} ; k ++ ))
do
    ./shock_wave ${TPB[k]} >> result.out
done

make clean
