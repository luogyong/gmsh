import os 
import scipy as sc
import scipy.io as sio
import numpy as np
import pylab as pl
import sys
pi    = np.pi
# os.system('rm '+str_path+'/*.pos')
# os.system('rm '+str_path+'/*.dat')
# os.system('rm '+str_path+'/*.out')
# os.system('rm '+str_path+'/*.txt')
# os.system('rm '+str_path+'/*.pre')
# os.system('rm '+str_path+'/*.res')


#### your configuration (nproc is your number of available processors)
str_gmsh_path    = '~/programs/gmsh/Gmsh-2.8.4-stable.app/Contents/MacOS/'
str_getdp_path   = '~/programs/getdp/getdp_build_svn/'
str_mesh_filename = 'haroche.msh'
str_pro_filename  = 'cavity_haroche.pro'

#### scaling the problem (*10^12)
nm	     = 1.e-9
epsilon0 = 8.854187817e-3*nm
mu0      = 400.*pi*nm
cel      = 1.0/(np.sqrt(epsilon0 * mu0))

#### mesh parameters : MeshElementSize = lambda/(paramaille*n)
lambda_m        = 5.8e-3
paramaille_air  = lambda_m/3.
paramaille_pml  = lambda_m/3.
paramaille_mir  = lambda_m/3.

#### opto-geometric parameters
#lambda_vp = 5.8e-3
#eig_target= (2.*np.pi*cel/lambda_vp)**2
freq_target = 51.099e9
lambda_vp = cel/freq_target
eig_target= (2.*np.pi*cel/lambda_vp)**2

neig = 70

#### write all above variables in a text file, used as an input to both gmsh and getdp
par_gmsh_getdp =open('parameters_gmsh_getdp.dat', 'w')
par_gmsh_getdp.write('nm      		 = %3.15e;\n'%(nm))        
par_gmsh_getdp.write('epsilon0  	 = %3.15e;\n'%(epsilon0))
par_gmsh_getdp.write('mu0  		     = %3.15e;\n'%(mu0))
par_gmsh_getdp.write('cel    		 = %3.15e;\n'%(cel))
par_gmsh_getdp.write('lambda_vp      = %3.15e;\n'%(lambda_vp))
par_gmsh_getdp.write('lambda_m       = %3.15e;\n'%(lambda_m))
par_gmsh_getdp.write('paramaille_air = %3.15e;\n'%(paramaille_air))
par_gmsh_getdp.write('paramaille_pml = %3.15e;\n'%(paramaille_pml))
par_gmsh_getdp.write('paramaille_mir = %3.15e;\n'%(paramaille_mir))
par_gmsh_getdp.write('eig_target     = %3.15e;\n'%(eig_target))
par_gmsh_getdp.write('neig       	 = %d    ;\n'%(neig))
par_gmsh_getdp.close()

### define strings to be "oscommanded"
str_getdp1        = str_getdp_path+'getdp '+str_pro_filename+' -pre all -msh  '+str_mesh_filename+' -cal -pos postop_eigenvectors_full -v 50 -bin'
str_slepc_options =  ' -eps_monitor -eps_view \
                -eps_type krylovschur \
               -st_ksp_type preonly \
               -st_pc_type lu \
               -st_pc_factor_mat_solver_package mumps \
               -eps_max_it 100 \
               -eps_mpd 400 \
               -eps_target '+ "%3.5e" % eig_target +' \
               -eps_target_real \
               -eps_nev '+"%d" % (neig)
               # \
               #-v 0'
# ### calling gmsh : mesh & save mesh!
os.system(str_gmsh_path+'gmsh geometry_haroche_realistic.geo -3 -o '+str_mesh_filename)
# ### calling getdp : preprocess and solve, using mpi if required
os.system(str_getdp1+str_slepc_options)

vp_real = np.loadtxt('./EigenValuesReal.txt',usecols=[5])
vp_imag = np.loadtxt('./EigenValuesImag.txt',usecols=[5])
print 'EigenFreq (GHz)', 2e9*np.pi*vp_real.transpose()
print 'vp imag', vp_imag.transpose()

os.system(str_gmsh_path+'gmsh haroche_postplot.geo &')
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
