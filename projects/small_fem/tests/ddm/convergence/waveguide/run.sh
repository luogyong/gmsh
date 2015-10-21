#!/bin/bash
MSH='../tests/ddm/waveguide/guide3dStruct.msh' # 20 mesh elements per lambda
O='1'
TYPE='scalar'
MODE='te'
DDM='emda'

OMP_NUM_THREADS=2 mpirun -np 2 ./ddm -msh $MSH -k 25 -max 1000 -chi 6.25 -lc 0.0125664 -ck 0 -pade 8 -hist -ov $O -ob $O -type $TYPE -mode $MODE -ddm $DDM -solver -ksp_rtol 1e-6 -ksp_monitor
