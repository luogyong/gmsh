--[[
    Script to launch beam problem with a lua script
  ]]
--[[
    data 
  ]]
-- material law
lawnum = 1 -- unique number of law
E = 1.e6 -- Young's modulus
nu = 0.3   -- Poisson's ratio

-- geometry
h = 0.1  -- thickness
meshfile="plate.msh" -- name of mesh file
-- integration
nsimp = 3 -- number of Simpson's points (odd)

-- solver
sol = 1 --Gmm=0 (default) Taucs=1 PETsc=2
beta1 = 100. -- value of stabilization parameter
beta2 = 100.
beta3 = 100.
soltype = 0 -- StaticLinear=0 (default) StaticNonLinear=1
nstep = 60   -- number of step (used only if soltype=1)
ftime =1.   -- Final time (used only if soltype=1)
tol=1.e-6   -- relative tolerance for NR scheme (used only if soltype=1)
nstepArch=6 -- Number of step between 2 archiving (used only if soltype=1)

--[[
    compute solution and BC (given directly to the solver
  ]]
-- creation of law
law1 = linearElasticLawPlaneStress(lawnum,E,nu)

-- creation of ElasticField
nfield =99 -- number of the field (physical number of surface)
fullDg = 1 --  formulation CgDg=0 fullDg =1
myfield1 = dgLinearShellDomain()
myfield1:tag(1000)
myfield1:thickness(h)
myfield1:simpsonPoints(nsimp)
myfield1:formulation(fullDg)
myfield1:lawnumber(lawnum)

-- creation of Solver
mysolver = DgC0PlateSolver(1000)
mysolver:readmsh(meshfile)
mysolver:addDgLinearElasticShellDomain(myfield1,nfield,2)
mysolver:AddLinearElasticLawPlaneStress(law1)
mysolver:setScheme(soltype)
mysolver:whichSolver(sol)
mysolver:SNLData(nstep,ftime,tol)
mysolver:stepBetweenArchiving(nstepArch)
mysolver:stabilityParameters(beta1,beta2,beta3)
-- BC
mysolver:prescribedDisplacement("Edge",11,0,0.)
mysolver:prescribedDisplacement("Edge",11,1,0.)
mysolver:prescribedDisplacement("Edge",21,2,0.)
mysolver:prescribedDisplacement("Edge",31,2,0.)
mysolver:prescribedDisplacement("Edge",41,0,0.)
mysolver:prescribedDisplacement("Edge",41,1,0.)
mysolver:prescribedForce("Node",111,0.,0.,-50.)
mysolver:AddThetaConstraint(11)
mysolver:AddThetaConstraint(41)
mysolver:solve()
