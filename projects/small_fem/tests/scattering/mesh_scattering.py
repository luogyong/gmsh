import os

### adapt path/names to ur config
str_gmsh_path     = ''#'/Users/demesy/programs/gmsh/Gmsh-2.8.4-stable.app/Contents/MacOS/'
str_geo_filename  = 'scattering'
str_mesh_filename = 'scattering'

### incident plane wave parameters
nm      = 1.e-9
lambda0 = 1000. * nm
theta0  = 0.0
phi0    = 0.0
psi0    = 0.0

### mesh size = lambda_medium/(paramaille)
paramaille_list = [5.0, 10.0]

### geom : sphere of radius ro of permittivity eps_In
###        inside medium of permittivity eps_Out in a cartesian PML box
dom       = 2000 * nm
dom_x     = dom
dom_y     = dom
dom_z     = dom
PML_top   = 600 * nm
PML_bot   = 600 * nm
PML_lat   = 600 * nm
ro        = 500 * nm
eps_In    = 9. - 1. * 1j
eps_Out   = 1.00

### loop on paramaille_list
while len(paramaille_list) != 0:
    paramaille = paramaille_list.pop()

    ## write common parameters file for gmsh&getdp
    par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'w')
    par_gmsh_getdp.write('paramaille = %3.15e;\n' %(paramaille))
    par_gmsh_getdp.write('nm         = %3.15e;\n' %(nm))
    par_gmsh_getdp.write('lambda0    = %3.15e;\n' %(lambda0))
    par_gmsh_getdp.write('theta0     = %3.15e;\n' %(theta0))
    par_gmsh_getdp.write('phi0       = %3.15e;\n' %(phi0))
    par_gmsh_getdp.write('psi0       = %3.15e;\n' %(psi0))
    par_gmsh_getdp.write('dom        = %3.15e;\n' %(dom))
    par_gmsh_getdp.write('dom_x      = %3.15e;\n' %(dom_x))
    par_gmsh_getdp.write('dom_y      = %3.15e;\n' %(dom_y))
    par_gmsh_getdp.write('dom_z      = %3.15e;\n' %(dom_z))
    par_gmsh_getdp.write('PML_top    = %3.15e;\n' %(PML_top))
    par_gmsh_getdp.write('PML_bot    = %3.15e;\n' %(PML_bot))
    par_gmsh_getdp.write('PML_lat    = %3.15e;\n' %(PML_lat))
    par_gmsh_getdp.write('ro         = %3.15e;\n' %(ro))
    par_gmsh_getdp.write('eps_re_In  = %3.15e;\n' %(eps_In.real))
    par_gmsh_getdp.write('eps_im_In  = %3.15e;\n' %(eps_In.imag))
    par_gmsh_getdp.write('eps_re_Out = %3.15e;\n' %(eps_Out.real))
    par_gmsh_getdp.write('eps_im_Out = %3.15e;\n' %(eps_Out.imag))
    par_gmsh_getdp.close()

    print 'Paramaille: %d' %(paramaille)
    os.system(str_gmsh_path+'gmsh '+str_geo_filename+'.geo'+\
              ' -3 -v 1 -o '+str_mesh_filename+'_%d' %(paramaille)+'.msh')
