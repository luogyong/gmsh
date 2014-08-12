import os
import sys
import time
import numpy

## Physical / Math constants ##
pi       = numpy.pi
nm       = 1.e-9
epsilon0 = 8.854187817e-3 * nm
mu0      = 400. * pi * nm
cel      = 1.0 / (numpy.sqrt(epsilon0 * mu0))

## User data ##
if(len(sys.argv) != 6):
    raise ValueError('Bad argument: '
                     'cavity_haroche_sf meshName symmetry femOrder nEig nProc')

meshName = str(sys.argv[1])
symmetry = int(sys.argv[2])
femOrder = int(sys.argv[3])
nEig     = int(sys.argv[4])
nProc    = int(sys.argv[5])

## My Data ##
tol   = 1e-6
maxIt = 100

## Eigen Problem Shift ##
freqTarget = 51.099e9
lambdaVp   = cel / freqTarget
eigTarget  = (2. * pi * cel / lambdaVp)**2

## Start ##
start = time.time()

print '## Haroche'
print ' -- Process       : %d' %nProc
print ' -- Mesh          : %s' %meshName
print ' -- FEM order     : %d' %femOrder
print ' -- N. Eigen      : %d' %nEig
print ' -- Shift         : %e' %eigTarget
print ' -- Symmetry      : %e' %symmetry
print ' -- Tolerance     : %e' %tol
print ' -- Max Iteration : %d' %maxIt
print '                      '

## Calling small_fem ##
print '## Simulating'
os.system('mpirun -np %d' %nProc     + \
          ' har'                     + \
          ' -msh %s'      %meshName  + \
          ' -o %d'        %femOrder  + \
          ' -n %d'        %nEig      + \
          ' -shift %e'    %eigTarget + \
          ' -sym %d'      %symmetry  + \
          ' -tol %e'      %tol       + \
          ' -maxit %d'    %maxIt     + \
          ' -solver -eps_monitor -eps_view')

## Renaming Output ##
for i in range(nProc):
    os.rename('harocheModes_proc%d.msh' %i         + \
              'post_proc%d_'            %i         + \
              'femorder_%d_'            %femOrder  + \
              'sym_%d_'                 %symmetry  + \
              meshName)

os.rename('harocheValues.txt', \
          'eig_femorder_%d_' %femOrder  + \
          'sym_%d_'          %symmetry  + \
          meshName                      + '.txt')

## Done ##
stop = time.time()
print '## Done in %e' %(stop - start) + ' [s]'
