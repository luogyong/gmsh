MESH="haroche_geoorder_2_air_5_pml_5_mir_5_part_4.msh"
PROCS="4"
ORDER="2"
NEIG="4"
SYM="0"
TOL="1e-15"
export OMP_NUM_THREADS=1

mpirun -np $PROCS har -msh $MESH -o $ORDER -n $NEIG -shift 1.03082e+23 -sym $SYM -tol $TOL -maxit 100 -pml pml.dat -solver -eps_monitor #-eps_view
