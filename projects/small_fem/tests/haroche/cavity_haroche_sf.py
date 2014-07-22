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
meshName  = str(sys.argv[1])
symmetry  = int(sys.argv[2])
fem_order = int(sys.argv[3])
neig      = int(sys.argv[4])

## My Data ##
tol   = 1e-6
maxIt = 100

## Eigen Problem Shift ##
freq_target = 51.099e9
lambda_vp   = cel / freq_target
eig_target  = (2. * pi * cel / lambda_vp)**2

## Start ##
start = time.time()

print '## Haroche'
print ' -- Mesh          : %s' %meshName
print ' -- FEM order     : %d' %fem_order
print ' -- N. Eigen      : %d' %neig
print ' -- Shift         : %e' %eig_target
print ' -- Symmetry      : %e' %symmetry
print ' -- Tolerance     : %e' %tol
print ' -- Max Iteration : %d' %maxIt
print '                      '

## Calling small_fem ##
print '## Simulating'
os.system('har'    + \
          ' -msh ' +       meshName   + \
          ' -o %d'        %fem_order  + \
          ' -n %d'        %neig       + \
          ' -shift %e'    %eig_target + \
          ' -sym %d'      %symmetry   + \
          ' -tol %e'      %tol        + \
          ' -maxit %d'    %maxIt      + \
          ' -solver -eps_monitor -eps_view')

## Renaming Output ##
os.rename('harocheModes.msh', \
          'post_femorder_%d_' %fem_order + \
          'sym_%d_'           %symmetry  + \
          meshName)

os.rename('harocheValues.txt', \
          'eig_femorder_%d_' %fem_order + \
          'sym_%d_'          %symmetry  + \
          meshName                      + '.txt')

## Done ##
stop = time.time()
print '## Done in %e' %(stop - start) + ' [s]'
