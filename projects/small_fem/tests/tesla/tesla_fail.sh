#!/bin/bash
NPROCS=4

OMP_NUM_THREADS=1 mpirun -bind-to none -np $NPROCS emode -msh ./fail/tesla_1cell_UQparam.msh -type vector -o 2 -n 10 -shift 0.00074 -tol 1e-6 -solver -eps_monitor
