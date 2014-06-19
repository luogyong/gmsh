#!/bin/bash

for i in {1..5}
do
    OMP_NUM_THREADS=6 \
        numactl --cpunodebind=0 --membind=0 \
        swave -msh circle3D.msh -k 16 -type vector -nopos -o $i
done
