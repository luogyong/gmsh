import os
import time
import numpy as np

## Physical / Math constant ##
pi       = np.pi
nm       = 1.e-9
epsilon0 = 8.854187817e-3 * nm
mu0      = 400. * pi * nm
cel      = 1.0 / (np.sqrt(epsilon0 * mu0))

## Mesh parameters : MeshElementSize = lambda / (paramaille * n) ##
mesh_filename   = 'haroche_mesh.msh'
refinement      = 3
lambda_m        = 5.8e-3

paramaille_air  = lambda_m / refinement
paramaille_pml  = lambda_m / refinement
paramaille_mir  = lambda_m / refinement

## Solver Parameters ##
neig      = 10
geo_order = 2
fem_order = 2
symmetry  = 1

## Eigen Problem Shift ##
freq_target = 51.099e9
lambda_vp   = cel / freq_target
eig_target  = (2. * pi * cel / lambda_vp)**2

## Start ##
start = time.time()

print '## Haroche'
print ' -- Geo order : %d' %geo_order
print ' -- Refinement: %d' %refinement
print ' -- FEM order : %d' %fem_order
print ' -- N. Eigen  : %d' %neig
print ' -- Shift     : %e' %eig_target
print '                  '

## Write Gmsh and GetDP data file ##
par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'w')

par_gmsh_getdp.write('nm             = %3.15e;\n' %(nm))
par_gmsh_getdp.write('epsilon0       = %3.15e;\n' %(epsilon0))
par_gmsh_getdp.write('mu0            = %3.15e;\n' %(mu0))
par_gmsh_getdp.write('cel            = %3.15e;\n' %(cel))
par_gmsh_getdp.write('lambda_vp      = %3.15e;\n' %(lambda_vp))
par_gmsh_getdp.write('lambda_m       = %3.15e;\n' %(lambda_m))
par_gmsh_getdp.write('paramaille_air = %3.15e;\n' %(paramaille_air))
par_gmsh_getdp.write('paramaille_pml = %3.15e;\n' %(paramaille_pml))
par_gmsh_getdp.write('paramaille_mir = %3.15e;\n' %(paramaille_mir))
par_gmsh_getdp.write('eig_target     = %3.15e;\n' %(eig_target))
par_gmsh_getdp.write('neig           = %d    ;\n' %(neig))

par_gmsh_getdp.close()

## Calling Gmsh ##
print '## Meshing'
os.system('gmsh geometry_haroche_realistic.geo -3 -order %d ' %geo_order + \
          '-o ' + mesh_filename + ' -v 3')

## Calling small_fem ##
print '## Simulating'
os.system('har'  + \
          ' -msh ' + mesh_filename    + \
          ' -o %d'        %fem_order  + \
          ' -n %d'        %neig       + \
          ' -shift %3.5e' %eig_target + \
          ' -sym %d'      %symmetry   + \
          ' -solver -eps_monitor')

## Renaming Output ##
os.rename('harocheModes.msh', \
          'haroche_femorder_%d' %fem_order  + \
          '_geoorder_%d'        %geo_order  + \
          '_paramaille_%d'      %refinement + \
          '_sym_%d'             %symmetry   + '.msh')

os.rename('harocheValues.txt', \
          'eigenValues_femorder_%d' %fem_order  + \
          '_geoorder_%d'            %geo_order  + \
          '_paramaille_%d'          %refinement + \
          '_sym_%d'                 %symmetry   + '.txt')

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
