#!/bin/python

## Imports
import onelab
import numpy
import os

## ONELAB
ol = onelab.client()

## Gmsh geometry
ol.openProject('cavity_haroche_2D.geo')
freq = float(ol.getNumber(name = 'Input/00Haroche/00Frequency'))

## Get SmallFEM data
femOrder = ol.defineNumber(name  = 'Input/03FEM/00Order',
                           value = 2)
nEig     = ol.defineNumber(name  = 'Input/04Eigenproblem/00Eigenvalue',
                           value = 4)
target   = ol.defineNumber(name  = 'Input/04Eigenproblem/01Target',
                           value = (freq * 2 * numpy.pi) ** 2, readOnly = 1)
nProc    = ol.defineNumber(name  = 'Input/05Solver/02Process',
                           value = 1)
tol      = ol.defineNumber(name  = 'Input/05Solver/03Tolerance',
                           value = 1e-12)
maxit    = ol.defineNumber(name  = 'Input/05Solver/04Iteration',
                           value = 100)
postpro  = ol.defineNumber(name  = 'Input/06Post-Pro/00Draw',
                           value = 1, choices = {0,1})

femOrder = int(femOrder)
nEig     = int(nEig)
target   = float(target)
nProc    = int(nProc)
tol      = float(tol)
maxit    = int(maxit)
postpro  = int(postpro)

## Check is done
if ol.action == 'check':
    exit(1)

## Mesh
ol.run('mesh', 'gmsh -2 -part ' + str(nProc) + ' cavity_haroche_2D.geo')
ol.mergeFile('cavity_haroche_2D.msh')

## Simulate
cmd  = 'mpirun -host localhost -np ' + str(nProc) + ' '
cmd += 'har2d'  + ' ' + \
       '-msh'   + ' ' + 'cavity_haroche_2D.msh' + ' ' + \
       '-o'     + ' ' + str(femOrder)           + ' ' + \
       '-n'     + ' ' + str(nEig)               + ' ' + \
       '-shift' + ' ' + str(target)             + ' ' + \
       '-tol'   + ' ' + str(tol)                + ' ' + \
       '-maxit' + ' ' + str(maxit)              + ' '

if postpro == 0:
    cmd += '-nopos' + ' '

cmd += '-solver -eps_monitor' # -eps_view'

os.environ['OMP_NUM_THREADS'] = str(1)
ol.call('sf', cmd)

## Draw
if postpro == 1:
    if(nProc == 1):
        ol.mergeFile('harocheModes.msh')

    else:
        for i in range(1, nProc + 1):
            ol.mergeFile('harocheModes_part' + str(i) + '.msh')

## Post pro
freq = numpy.loadtxt('harocheValues.txt', usecols = [3], skiprows = 1)
time = numpy.loadtxt('harocheValues.txt', usecols = [4], skiprows = 1)
maxi = numpy.argmax(time)

ol.sendInfo('')
ol.sendInfo('Found eigen frequencies [MHz] and life time [ms]:')
for i in range(0, len(freq)):
    ol.sendInfo(str(i) + ': ' + \
                str(freq[i])  + ' [MHZ]\t' + str(time[i]) + ' [ms]')


ol.sendInfo('')
ol.sendInfo('Best candidate is:')
ol.sendInfo(str(maxi) + ': ' + \
            str(freq[maxi])  + ' [MHZ]\t' + str(time[maxi]) + ' [ms]')

run = int(ol.defineNumber(name = 'Output/00Run', value = 0, visible = 0))

ol.defineNumber(name = 'Output/Freq ' + str(run), value=freq[maxi], readOnly=1)
ol.defineNumber(name = 'Output/Time ' + str(run), value=time[maxi], readOnly=1)

ol.setNumber(name = 'Output/00Run', value = run + 1)
