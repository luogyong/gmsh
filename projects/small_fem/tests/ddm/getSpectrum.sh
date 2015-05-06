#!/bin/bash

#MSH='./circle/circle_concentric.msh'
MSH='./cylinder/cylinder_concentric.msh'
TYPE='vector'
DDM='osrc'
ROOT='spectrumVectCir3D'
K=5
NDOM=2
OMP_NUM_THREADS=2

mkdir $ROOT

for V in $(seq 4 4)
do
    for B in $(seq 0 $V)
    do
        echo 'Volume: '$V' - Border: '$B
        mpirun -x OMP_NUM_THREADS -np $NDOM cirsp -msh $MSH -type $TYPE -ddm $DDM -k $K -ov $V -ob $B -max 1000 -ck 0 -pade 4 -name $ROOT'/'$V$B -solver -eps_monitor
    done
done
echo 'Done!'
