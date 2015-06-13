#!/bin/bash

## Binaries
ROOT=''
BINANA=$ROOT'wavea'
BINDDM=$ROOT'ddm'

## Data
THREADS=2
TYPE='vector'
MODE='tm'
DDM='emda'
NDOM=2
NMSH=6
MSH='./guide2d_'
INTERP='./guide2d_32.msh'
OV=1
OB=1
K=20

## Useful
NDOMMinus=$(bc <<< $NDOM'-1')

## Analytical Solution
echo '#### Reference ####'
$BINANA -msh $INTERP -k $K -n $NDOM -type $TYPE -mode $MODE -name 'ref'

## DDM
for m in $(seq 1 $NMSH)
do
    MYMESH=$MSH$(bc <<< 2^'('$m-1')')'.msh'
    for v in $(seq 1 $OV)
    do
        for b in $(seq 1 $v)
        do
            NAME=$TYPE'_'$DDM'_'$m'_'$v'_'$b
            echo '#### '$NAME' ####'

            OMP_NUM_THREADS=$THREADS
            mpirun -bind-to none \
                   -np $NDOM \
                   $BINDDM -msh $MYMESH -k $K \
                   -max 250 -ddm $DDM -mode $MODE \
                   -pade 4 -ck 0 -chi 0 -lc 0.06 -type $TYPE -ov $v -ob $b \
                   -name $NAME -interp $INTERP -hist $NAME'.hist' \
                   -solver -ksp_rtol 1e-9 -ksp_monitor

            # -machinefile machine
            ## Synching nodes
            # rsync -avPH ace43:ddm/ .

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

## Grab every history in a single file
FILE=$TYPE'_'$DDM'_'$NMSH'_'$OV'_'$OB'.hist.dat'
OCTAVE='grabHist(''"'$TYPE'"'', ''"'$DDM'"'', '
OCTAVE=$OCTAVE$NMSH', '$OV', '$OB', '
OCTAVE=$OCTAVE'"'$FILE'"'')'
octave -q <<< "$OCTAVE"

## Clear temp files
for m in $(seq 1 $NMSH)
do
    for v in $(seq 1 $OV)
    do
        for b in $(seq 1 $v)
        do
            rm $TYPE'_'$DDM'_'$m'_'$v'_'$b'.hist'

            for d in $(seq 0 $NDOMMinus)
            do
                rm $TYPE'_'$DDM'_'$m'_'$v'_'$b'_'$d'.tmp'
            done
        done
    done
done
