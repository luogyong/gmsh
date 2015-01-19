#!/bin/bash

## DATA
TYPE='vector'
DOM=2
MSH='./mesh_all.msh'
REF='./fine.msh'
OR=4
K=5

## Useful
DOMMinus=$(bc <<< $DOM'-1')

## Reference Solution
echo '#### Reference ####'
waveg -msh $MSH -k $K -n $DOM -type $TYPE -interp $REF -o $OR -name 'ref'

## Lower order solution
for o in $(seq 1 $OR);
do
    NAME='order_'$o
    echo '#### '$NAME' ####'

    waveg -msh $MSH -k $K -n $DOM -type $TYPE -interp $REF -o $o -name $NAME
    octave -q l2.m $NAME'0.dat' 'ref0.dat' > $NAME'.tmp'

    ## Clear
    for i in $(seq 0 $DOMMinus);
    do
        rm $NAME$i'.dat'
    done
done

## Clear Reference
for i in $(seq 0 $DOMMinus);
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
