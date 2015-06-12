#include<cstring>
#include "dG3DMaterialLaw.h"
#include "nonLinearMechSolver.h"
#include "dG3DDomain.h"
#include "GmshMessage.h"
#include "Context.h"
int main(int argc, char ** argv)
{
  Msg::Init(argc,argv);
  CTX::instance()->terminal = 1; // How to avoid this?

 	# material law
	lawnumEP = 2 
	rhoEP   = 1000.
	youngEP = 3.2e9
	nuEP    = 0.3
	sy0EP   = 25.e6
	hEP     = 71.e9


	lawcnumEP = 4
	GcEP   =  83.
	sigmacEP = 78.e6
	deltacEP = 2*GcEP/sigmacEP
	KcEP = sigmacEP/(0.001*deltacEP)



	beta =0.87 # ratio KII/KI
	mu = 0. # friction coefficient ??
	fsmin = 0.9
	fsmax = 1.1

	# geometry
	meshfile="twoHole.msh" # name of mesh file

	# solver
	sol = 2  # Gmm=0 (default) Taucs=1 PETsc=2
	soltype = 1 # StaticLinear=0 (default) StaticNonLinear=1
	nstep = 300   # number of step (used only if soltype=1)
	ftime =1 # Final time (used only if soltype=1)
	tol=1.e-4   # relative tolerance for NR scheme (used only if soltype=1)
	nstepArch=10 # Number of step between 2 archiving (used only if soltype=1)
	fullDg =1 #O = CG, 1 = DG
	space1 = 0 # function space (Lagrange=0)
	beta1  = 100

	#  compute solution and BC (given directly to the solver
	# creation of law
	law1EP = J2LinearDG3DMaterialLaw(lawnumEP,rhoEP,youngEP,nuEP,sy0EP,hEP)
	law2EP = LinearCohesive3DLaw(lawcnumEP,GcEP,sigmacEP,beta,mu,fsmin,fsmax,KcEP)
	law3EP = FractureByCohesive3DLaw(6,lawnumEP,lawcnumEP)


	# creation of ElasticField

	myfieldEP = dG3DDomain(1000,54,space1,6,fullDg)
	#myfieldEP.matrixByPerturbation(1,1,1,1e-8)
	myfieldEP.stabilityParameters(beta1)


	# creation of Solver
	mysolver = nonLinearMechSolver(1000)
	mysolver.loadModel(meshfile)
	mysolver.addDomain(myfieldEP)
	mysolver.addMaterialLaw(law1EP)
	mysolver.addMaterialLaw(law2EP)
	mysolver.addMaterialLaw(law3EP)
	mysolver.Scheme(soltype)
	mysolver.Solver(sol)
	mysolver.snlData(nstep,ftime,tol)
	mysolver.snlManageTimeStep(10,5,20,10)
	mysolver.stepBetweenArchiving(nstepArch)


	#mysolver.explicitSpectralRadius(ftime,0.5,0.)
	#mysolver.dynamicRelaxation(0.1, ftime, 1.e-3,2)
	#mysolver.explicitTimeStepEvaluation(nstep)


	#mysolver.pathFollowing(0)
	#mysolver.setPathFollowingControlType(0)

	# BC

	mysolver.displacementBC("Volume",54,2,0.)

	mysolver.displacementBC("Face",55,0,0.)
	mysolver.displacementBC("Face",55,1,0.)

	mysolver.displacementBC("Face",56,0,0.)
	mysolver.displacementBC("Face",56,1,3e-6)



	mysolver.internalPointBuildView("svm",IPField.SVM, 1, 1)
	mysolver.internalPointBuildView("sig_xx",IPField.SIG_XX, 1, 1)
	mysolver.internalPointBuildView("sig_yy",IPField.SIG_YY, 1, 1)
	mysolver.internalPointBuildView("sig_zz",IPField.SIG_ZZ, 1, 1)
	mysolver.internalPointBuildView("sig_xy",IPField.SIG_XY, 1, 1)
	mysolver.internalPointBuildView("sig_yz",IPField.SIG_YZ, 1, 1)
	mysolver.internalPointBuildView("sig_xz",IPField.SIG_XZ, 1, 1)
	mysolver.internalPointBuildView("epl",IPField.PLASTICSTRAIN, 1, 1)

	mysolver.archivingForceOnPhysicalGroup("Face", 55, 1, 1)
	mysolver.archivingNodeDisplacement(57,1,1)

  Msg::Info("Init done");
  mysolver.solve();

  Msg::Exit(0);

  return 0;
}
