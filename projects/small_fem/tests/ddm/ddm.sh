#!/bin/bash

## Binaries
BINREF='waveg'
BINDDM='ddm'

## Data
THREADS=2
TYPE='vector'
NDOM=2
DDM='osrc'
NMSH=3
MSH='./guide3d_'
REF='./guide3d_16.msh'
OR=3
OV=3
OB=3
K=10

## Useful
NDOMMinus=$(bc <<< $NDOM'-1')

## Reference Solution
echo '#### Reference ####'
MYMESH=$MSH$(bc <<< 2^$NMSH)'.msh'
$BINREF -msh $MYMESH -k $K -n $NDOM -type $TYPE -o $OR -sigma 0 \
        -interp $REF -name 'ref'

## DDM
for m in $(seq 1 $NMSH)
do
    MYMESH=$MSH$(bc <<< 2^$m)'.msh'
    for v in $(seq 1 $OV)
    do
        for b in $(seq 1 $OB)
        do
            NAME=$TYPE'_'$DDM'_'$m'_'$v'_'$b
            echo '#### '$NAME' ####'

            OMP_NUM_THREADS=$THREADS
            mpirun -np $NDOM $BINDDM -msh $MYMESH -k $K \
                   -max 250 -ddm $DDM \
                   -pade 4 -ck 0 -chi 0 -lc 0.06 -type $TYPE -ov $v -ob $b \
                   -name $NAME -interp $REF -hist $NAME'.hist' \
                   -solver -ksp_rtol 1e-9

            ## L2 Error
            for d in $(seq 0 $NDOMMinus)
            do
                octave -q l2.m $NAME$d'.dat' 'ref'$d'.dat' > $NAME'_'$d'.tmp'
            done

            ## Clear
            for d in $(seq 0 $NDOMMinus)
            do
                rm $NAME$d'.dat'
            done
        done
    done
done

## Clear Reference
for d in $(seq 0 $NDOMMinus);
do
    rm 'ref'$d'.dat'
done

## Grab every residual in a single file
FILE=$TYPE'_'$DDM'_'$NMSH'_'$OV'_'$OB'.l2.dat'
OCTAVE='grabL2(''"'$TYPE'"'', ''"'$DDM'"'', '
OCTAVE=$OCTAVE$NDOM', '$NMSH', '$OV', '$OB', '
OCTAVE=$OCTAVE'"'$FILE'"'')'
octave -q <<< "$OCTAVE"

## Clear temp files
for m in $(seq 1 $NMSH)
do
    for v in $(seq 1 $OV)
    do
        for b in $(seq 1 $OB)
        do
            for d in $(seq 0 $NDOMMinus)
            do
                NAME=$TYPE'_'$DDM'_'$m'_'$v'_'$b'_'$d'.tmp'
                rm $NAME
            done
        done
    done
done
