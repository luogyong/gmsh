#coding-Utf-8-*-
from gmshpy import *
from dG3Dpy import*

#script to launch beam problem with a python script

# material law
lawnum   = 1 # unique number of law
rho      = 7850.
young    = 1.e10
nu       = 0. #

# geometry
meshfile="cube.msh" # name of mesh file

# solver
sol = 2  # Gmm=0 (default) Taucs=1 PETsc=2
soltype = 1 # StaticLinear=0 (default) StaticNonLinear=1
nstep = 1   # number of step (used only if soltype=1)
ftime =1.   # Final time (used only if soltype=1)
tol=1.e-6   # relative tolerance for NR scheme (used only if soltype=1)
nstepArch=1 # Number of step between 2 archiving (used only if soltype=1)
fullDg = 1 #O = CG, 1 = DG
space1 = 0 # function space (Lagrange=0)


#  compute solution and BC (given directly to the solver
# creation of law
law1 = J2LinearDG3DMaterialLaw(lawnum,rho,young,nu,1.e100,0.)

# creation of ElasticField
nfield = 101 # number of the field (physical number of surface)
myfield1 = dG3DDomain(1000,nfield,space1,lawnum,fullDg)
#myfield1.stabilityParameters(30.)
# creation of Solver
mysolver = nonLinearMechSolver(1000)
mysolver.loadModel(meshfile)
mysolver.addDomain(myfield1)
mysolver.addMaterialLaw(law1)
mysolver.Scheme(soltype)
mysolver.Solver(sol)
mysolver.snlData(nstep,ftime,tol)
mysolver.stepBetweenArchiving(nstepArch)

mysolver.displacementBC("Face",2376,0,0.0)
mysolver.displacementBC("Face",2376,1,0.0)
mysolver.displacementBC("Face",2376,2,0.0)
mysolver.pressureOnPhysicalGroupBC(1584,1.e6,0.)

mysolver.internalPointBuildView("svm",IPField.SVM, 1, 1)

mysolver.solve()

