#!/bin/bash

## Binaries
BINREF='disc'

## Data
TYPE='vector'
NDOM=4
DOM=0
MSH='./mesh_all.msh'
REF='./fine.msh'
OR=3
K=5

## Useful
NDOMMinus=$(bc <<< $NDOM'-1')

## Reference Solution
echo '#### Reference ####'
$BINREF -msh $MSH -k $K -n $NDOM -type $TYPE -interp $REF -o $OR -name 'ref'

## Lower order solution
for o in $(seq 1 $OR);
do
    NAME='order_'$o
    echo '#### '$NAME' ####'

    $BINREF -msh $MSH -k $K -n $NDOM -type $TYPE -interp $REF -o $o -name $NAME
    octave -q l2.m $NAME$DOM'.dat' 'ref'$DOM'.dat' > $NAME'.tmp'

    ## Clear
    for i in $(seq 0 $NDOMMinus);
    do
        rm $NAME$i'.dat'
    done
done

## Clear Reference
for i in $(seq 0 $NDOMMinus);
do
    rm 'ref'$i'.dat'
done

## Grab every residual in a single file
FILE=$'ref_'$OR'.l2'
echo $FILE             > $FILE
echo 'Row: FEM order' >> $FILE

for o in $(seq 1 $OR);
do
    NAME='order_'$o'.tmp'
    echo $(cat $NAME) >> $FILE
done

## Clear temp files
for o in $(seq 1 $OR);
do
    NAME='order_'$o'.tmp'
    rm $NAME
done
