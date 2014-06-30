import os
import sys
import time
import numpy

## Physical / Math constant ##
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

## Eigen Problem Shift ##
freq_target = 51.099e9
lambda_vp   = cel / freq_target
eig_target  = (2. * pi * cel / lambda_vp)**2

## Start ##
start = time.time()

print '## Haroche'
print ' -- Mesh      : %s' %meshName
print ' -- FEM order : %d' %fem_order
print ' -- N. Eigen  : %d' %neig
print ' -- Shift     : %e' %eig_target
print ' -- symmetry  : %e' %symmetry
print '                  '

## Calling small_fem ##
print '## Simulating'
os.system('har'    + \
          ' -msh ' + meshName         + \
          ' -o %d'        %fem_order  + \
          ' -n %d'        %neig       + \
          ' -shift %3.5e' %eig_target + \
          ' -sym %d'      %symmetry   + \
          ' -solver -eps_monitor')

## Renaming Output ##
os.rename('harocheModes.msh', \
          'post_femorder_%d_' %fem_order + \
          'sym_%d_'           %symmetry  + \
          meshName)

          # 'haroche_femorder_%d' %fem_order  + \
          # '_geoorder_%d'        %geo_order  + \
          # '_paramaille_%d'      %refinement + \
          # '_sym_%d'             %symmetry   + '.msh')

os.rename('harocheValues.txt', \
          'eig_femorder_%d_' %fem_order + \
          'sym_%d_'          %symmetry  + \
          meshName                      + '.txt')

          # 'eigenValues_femorder_%d' %fem_order  + \
          # '_geoorder_%d'            %geo_order  + \
          # '_paramaille_%d'          %refinement + \
          # '_sym_%d'                 %symmetry   + '.txt')

## Done ##
stop = time.time()
print '## Done in %e' %(stop - start) + ' [s]'

## Post Pro ##
#vp_real = np.loadtxt('./EigenValuesReal.txt',usecols=[5])
#vp_imag = np.loadtxt('./EigenValuesImag.txt',usecols=[5])
#print 'EigenFreq (GHz)', 2e9*np.pi*vp_real.transpose()
#print 'vp imag', vp_imag.transpose()

#os.system(str_gmsh_path+'gmsh haroche_postplot.geo &')
# os.system(str_gmsh_path+'gmsh eigenVectors_CompX_faces.pos geometry_haroche_realistic.geo &')

# omega2=(vp_real+1j*vp_imag)**2
# pl.plot(np.sqrt(-omega2.real).sort())
# print np.sort(2*pi*cel/vp_real*1e9)
# pl.figure()
# pl.plot(vp_real,vp_imag,'o')
# pl.figure()
# pl.plot(-vp_real,-vp_imag,'o')
# pl.show()

# omega_real_sort = np.sort(-omega2.real)
# pl.figure()
# pl.plot(omega_real_sort,'o')
# pl.figure()
# pl.plot(omega2.imag,'-')
# pl.show()
# vp_real_adj = np.loadtxt('./EigenValuesReal_adj.txt',usecols=[5])
# vp_imag_adj = np.loadtxt('./EigenValuesImag_adj.txt',usecols=[5])
