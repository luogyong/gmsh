#!/usr/bin/env python
from gmshpy import *
import sys
import os

## Get Data ##
##############
if(len(sys.argv) != 3):
    raise ValueError('Bad argument: '
                     'partition filename number_of_partitions')

name = str(sys.argv[1])
part = int(sys.argv[2])

## Read Mesh ##
###############
model = GModel()
model.readMSH(name)

## Partition ##
###############
partitionOpt = meshPartitionOptions()
partitionOpt.setNumOfPartitions(part)
PartitionMesh(model, partitionOpt)

## Save ##
##########
model.save(os.path.splitext(name)[0] + '_part_' + str(part) + '.msh')
