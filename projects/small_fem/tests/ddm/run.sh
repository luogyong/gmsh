#!/bin/bash

## DATA
TYPE='scalar'
DDM='emda'
DOM=4
MSH='./cylinder/out/mesh_all.msh'
REF='./cylinder/ref.msh'
OR=4
OV=2
OB=2
K=5

## Useful
DOMMinus=$(bc <<< $DOM'-1')

## Reference Solution
echo '#### Reference ####'
disc -msh $MSH -k $K -n $DOM -type $TYPE -interp $REF -o $OR -name 'ref'

## DDM
for v in $(seq 1 $OV);
do
    for b in $(seq 1 $OB);
    do
        NAME=$TYPE'_'$DDM'_'$v'_'$b
        echo '#### '$NAME' ####'

        OMP_NUM_THREADS=1 mpirun -np $DOM cir -msh $MSH -k $K \
            -max 1000 -ddm $DDM \
            -pade 4 -ck 0 -chi 0 -lc 0.06 -type $TYPE -ov $v -ob $b \
            -name $NAME -interp $REF -hist $NAME'.hist' -solver -ksp_rtol 1e-9

        ## L2 Error
        octave -q l2.m $NAME'0.dat' 'ref0.dat' > $NAME'.tmp'

        ## Clear
        for i in $(seq 0 $DOMMinus);
        do
            rm $NAME$i'.dat'
        done
    done
done

## Clear Reference
for i in $(seq 0 $DOMMinus);
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
