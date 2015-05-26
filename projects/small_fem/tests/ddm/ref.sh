#!/bin/bash

## Binaries
ROOT=''
BINANA=$ROOT'wavea'
BINFEM=$ROOT'waveg'

## Data
THREADS=4
TYPE='vector'
NDOM=2
NMSH=5
MSH='./guide2d_'
INTERP='./guide2d_64.msh'
nO=4
K=100

## Useful
NDOMMinus=$(bc <<< $NDOM'-1')

## Analytical Solution
echo '#### Reference ####'
$BINANA -msh $INTERP -k $K -n $NDOM -type $TYPE -name 'ref'

## Iterate on mesh
for m in $(seq 1 $NMSH)
do
    MYMESH=$MSH$(bc <<< 2^$m)'.msh'

    ## Iterate on orders
    for o in $(seq 1 $nO)
    do
        NAME=$TYPE'_'$m'_'$o
        echo '#### '$NAME' ####'

        ## FEM solution
        OMP_NUM_THREADS=$THREADS
        $BINFEM -msh $MYMESH -k $K -o $o -n $NDOM \
                -type $TYPE -interp $INTERP -name $NAME

        ## L2 Error
        for d in $(seq 0 $NDOMMinus)
        do
            octave -q l2.m $NAME$d'.dat' 'ref'$d'.dat' > $NAME'_'$d'.tmp'
        done

        ## Clear
        for d in $(seq 0 $NDOMMinus);
        do
            rm $NAME$d'.dat'
        done
    done
done

## Clear Reference
for d in $(seq 0 $NDOMMinus)
do
    rm 'ref'$d'.dat'
done

## Grab every residual in a single file
FILE=$TYPE'_'$NMSH'_'$nO'.l2.dat'
OCTAVE='grabL2Lite(''"'$TYPE'"'', '
OCTAVE=$OCTAVE$NDOM', '$NMSH', '$nO', '
OCTAVE=$OCTAVE'"'$FILE'"'')'
octave -q <<< "$OCTAVE"

## Clear temp files
for m in $(seq 1 $NMSH)
do
    for o in $(seq 1 $nO)
    do
        for d in $(seq 0 $NDOMMinus)
        do
            NAME=$TYPE'_'$m'_'$o'_'$d'.tmp'
            rm $NAME
        done
    done
done
