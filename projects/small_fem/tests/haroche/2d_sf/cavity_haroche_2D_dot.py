#!/bin/python

## Imports
import onelab
import numpy

## ONELAB
ol = onelab.client()

## Gmsh geometry
ol.openProject('cavity_haroche_2D.geo')
freq   = float(ol.getNumber(name = 'Input/00Haroche/00Frequency'))
target = (freq * 2 * numpy.pi) ** 2

## For remote mesh
S_PML = int(ol.getNumber(name = 'Input/01Geometry/00PML size'))
D_PML = int(ol.getNumber(name = 'Input/01Geometry/01PML distance'))
MSH_A = int(ol.getNumber(name = 'Input/02Mesh/00Size Air'))
MSH_P = int(ol.getNumber(name = 'Input/02Mesh/00Size PML'))
MSH_M = int(ol.getNumber(name = 'Input/02Mesh/00Size Mirror'))
ORDER = int(ol.getNumber(name = 'Input/02Mesh/01Order'))
TYPE  = int(ol.getNumber(name = 'Input/02Mesh/02Type'))

## Get SmallFEM data
femOrder = ol.defineNumber(name = 'Input/03FEM/00Order',        value = 1)
nEig     = ol.defineNumber(name = 'Input/03FEM/01Eigenvalue',   value = 4)
sym      = ol.defineNumber(name = 'Input/03FEM/02Symmetry',     choices = {0,1})
nProc    = ol.defineNumber(name = 'Input/04Solver/00Process',   value = 12)
tol      = ol.defineNumber(name = 'Input/04Solver/01Tolerance', value = 1e-6)
maxit    = ol.defineNumber(name = 'Input/04Solver/02Iteration', value = 10000)
postpro  = ol.defineNumber(name = 'Input/05Post-Pro/00Draw',    choices = {0,1})

femOrder = int(femOrder)
nEig     = int(nEig)
sym      = int(sym)
nProc    = int(nProc)
tol      = float(tol)
maxit    = int(maxit)
postpro  = int(postpro)

## Check is done
if ol.action == 'check':
    exit(1)

## Mesh
mesh = '/home/nicolas/src/gmsh/build/'             + \
       'gmsh'                                + ' ' + \
       '-2'                                  + ' ' + \
       '-part'            + ' ' + str(nProc) + ' ' + \
       '-bin'                                + ' ' + \
       '-setnumber F_HAR' + ' ' + str(freq)  + ' ' + \
       '-setnumber S_PML' + ' ' + str(S_PML) + ' ' + \
       '-setnumber D_PML' + ' ' + str(D_PML) + ' ' + \
       '-setnumber MSH_A' + ' ' + str(MSH_A) + ' ' + \
       '-setnumber MSH_P' + ' ' + str(MSH_P) + ' ' + \
       '-setnumber MSH_M' + ' ' + str(MSH_M) + ' ' + \
       '-setnumber ORDER' + ' ' + str(ORDER) + ' ' + \
       '-setnumber TYPE'  + ' ' + str(TYPE)  + ' ' + \
       'cavity_haroche_2D.geo'

ol.upload('cavity_haroche_2D.geo', 'onelab/', 'ace42')
ol.call('mesh', mesh, 'ace42', 'onelab/')
#ol.run('mesh', 'gmsh -2 -part ' + str(nProc) + ' cavity_haroche_2D.geo')
ol.download('onelab/cavity_haroche_2D.msh', '.', 'ace42')
ol.mergeFile('cavity_haroche_2D.msh')

## Simulate
cmd  = 'mpirun -host localhost -np ' + str(nProc) + ' ' + \
       '-x OMP_NUM_THREADS=1'                     + ' '

cmd += '/home/nicolas/bin/har2d'                  + ' ' + \
       '-msh'   + ' ' + 'cavity_haroche_2D.msh'   + ' ' + \
       '-o'     + ' ' + str(femOrder)             + ' ' + \
       '-n'     + ' ' + str(nEig)                 + ' ' + \
       '-shift' + ' ' + str(target)               + ' ' + \
       '-sym'   + ' ' + str(sym)                  + ' ' + \
       '-tol'   + ' ' + str(tol)                  + ' ' + \
       '-maxit' + ' ' + str(maxit)                + ' '

if postpro == 0:
    cmd += '-nopos' + ' '

cmd += '-solver -eps_monitor' # -eps_view'

ol.call('sf', cmd, 'ace42', 'onelab/')

## Draw
if postpro == 1:
    for i in range(0, nProc):
        ol.download('onelab/harocheModes_proc' + str(i) + '.msh', '.', 'ace42')

    for i in range(0, nProc):
        ol.mergeFile('harocheModes_proc' + str(i) + '.msh')

## Post pro
ol.download('onelab/harocheValues.txt', '.', 'ace42')

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
