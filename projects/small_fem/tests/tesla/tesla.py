#!/bin/python

## Imports
import onelab
import numpy
import os

## ONELAB
ol = onelab.client()

## Gmsh geometry
name = 'tesla_1cell_UQparam'
ol.openProject(name + '.geo')

## Get SmallFEM data
mshOrder = ol.defineNumber(name = 'Input/00Mesh/00Order',       value = 2)
femOrder = ol.defineNumber(name = 'Input/01FEM/00Order',        value = 2)
target   = ol.defineNumber(name = 'Input/02Eigen/00Target',     value = 0.00074)
nEig     = ol.defineNumber(name = 'Input/02Eigen/01Number',     value = 10)
postpro  = ol.defineNumber(name = 'Input/03Post-Pro/00Draw',    choices = {0,1},
                                                                  value = 1)
nProc    = ol.defineNumber(name = 'Input/04Solver/00Process',   value = 4)
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
cmd  = 'mpirun -host localhost -np ' + str(nProc) + ' '
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

os.environ['OMP_NUM_THREADS'] = str(1)
ol.call('sf', cmd)

## Draw
if postpro == 1:
    for i in range(0, nProc):
        ol.mergeFile('eigen_mode_proc' + str(i) + '.msh')
