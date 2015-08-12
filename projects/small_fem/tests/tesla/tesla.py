#!/bin/python

## Imports
import onelab
import numpy
import os

## ONELAB
ol = onelab.client()

## Gmsh geometry
name = 'tesla'
ol.openProject(name + '.geo')

## Get SmallFEM data
mshOrder = ol.defineNumber(name = 'Input/00Mesh/00Order',       value = 2)
femOrder = ol.defineNumber(name = 'Input/01FEM/00Order',        value = 2)
target   = ol.defineNumber(name = 'Input/02Eigen/00Target',     value = 7e-4)
nEig     = ol.defineNumber(name = 'Input/02Eigen/01Number',     value = 10)
postpro  = ol.defineNumber(name = 'Input/03Post-Pro/00Draw',    choices = {0,1},
                                                                  value = 1)
nProc    = ol.defineNumber(name = 'Input/04Solver/00Process',   value = 1)
tol      = ol.defineNumber(name = 'Input/04Solver/01Tolerance', value = 1e-6)
#maxit    = ol.defineNumber(name = 'Input/04Solver/02Iteration', value = 10000)

mshOrder = int(mshOrder)
femOrder = int(femOrder)
target   = float(target)
nEig     = int(nEig)
postpro  = int(postpro)
nProc    = int(nProc)
tol      = float(tol)
#maxit    = int(maxit)

## Check is done
if ol.action == 'check':
    exit(1)

## Mesh
ol.run('mesh',
       'gmsh'         + ' ' + \
       '-3 '          + ' ' + \
       '-order'       + ' ' + str(mshOrder) + ' ' + \
       '-part'        + ' ' + str(nProc)    + ' ' + \
       '-optimize_ho' + ' ' + str(name)     + '.geo')
ol.mergeFile(name + '.msh')

## Simulate
cmd  = 'mpirun --bind-to none -np '   + str(nProc) + ' '
cmd += 'emode'  + ' ' + \
       '-msh'   + ' ' + name + '.msh' + ' ' + \
       '-type'  + ' ' + 'vector'      + ' ' + \
       '-o'     + ' ' + str(femOrder) + ' ' + \
       '-n'     + ' ' + str(nEig)     + ' ' + \
       '-shift' + ' ' + str(target)   + ' ' + \
       '-tol'   + ' ' + str(tol)      + ' '
#       '-maxit' + ' ' + str(maxit)    + ' '

if postpro == 0:
    cmd += '-nopos' + ' '

cmd += '-solver -eps_monitor' # -eps_view'

os.environ['OMP_NUM_THREADS'] = str(4)
ol.call('sf', cmd)

## Draw
if postpro == 1:
    if(nProc == 1):
        ol.mergeFile('eigenModes.msh')

    else:
        for i in range(1, nProc + 1):
            ol.mergeFile('eigenModes_part' + str(i) + '.msh')
