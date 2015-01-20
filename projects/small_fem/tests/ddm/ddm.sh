#!/bin/bash

## Binaries
BINREF='disc'
BINDDM='cir'

## Data
THREADS=1
TYPE='vector'
NDOM=4
DOM=0
DDM='osrc'
MSH='./mesh_all.msh'
REF='./fine.msh'
OR=3
OV=3
OB=3
K=5

## Useful
NDOMMinus=$(bc <<< $NDOM'-1')

## Reference Solution
echo '#### Reference ####'
$BINREF -msh $MSH -k $K -n $NDOM -type $TYPE -interp $REF -o $OR -name 'ref'

## DDM
for v in $(seq 1 $OV);
do
    for b in $(seq 1 $OB);
    do
        NAME=$TYPE'_'$DDM'_'$v'_'$b
        echo '#### '$NAME' ####'

        OMP_NUM_THREADS=$THREADS mpirun -np $NDOM $BINDDM -msh $MSH -k $K \
            -max 250 -ddm $DDM \
            -pade 4 -ck 0 -chi 0 -lc 0.06 -type $TYPE -ov $v -ob $b \
            -name $NAME -interp $REF -hist $NAME'.hist' -solver -ksp_rtol 1e-9

        ## L2 Error
        octave -q l2.m $NAME$DOM'.dat' 'ref'$DOM'.dat' > $NAME'.tmp'

        ## Clear
        for i in $(seq 0 $NDOMMinus);
        do
            rm $NAME$i'.dat'
        done
    done
done

## Clear Reference
for i in $(seq 0 $NDOMMinus);
do
    rm 'ref'$i'.dat'
done

## Grab every residual in a single file
FILE=$TYPE'_'$DDM'_'$OV'_'$OB'.l2'
echo $FILE                > $FILE
echo 'Row: Volume order' >> $FILE
echo 'Col: Border order' >> $FILE

for v in $(seq 1 $OV);
do
    TMP=''
    for b in $(seq 1 $OB);
    do
        NAME=$TYPE'_'$DDM'_'$v'_'$b'.tmp'
        TMP=$TMP' '$(echo $(cat $NAME))
    done
    echo $TMP >> $FILE
done

## Clear temp files
for v in $(seq 1 $OV);
do
    for b in $(seq 1 $OB);
    do
        NAME=$TYPE'_'$DDM'_'$v'_'$b'.tmp'
        rm $NAME
    done
done
