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

  // material law
  int lawnum   = 1;
  double rho      = 7850.;
  double young    = 1.e10;
  double nu       = 0.;

  //geometry
  std::string meshfile("cube.msh");

  //solver
  int sol = 2;
  int soltype = 1;
  int nstep = 1;
  double ftime =1.;
  double tol=1.e-6;
  int nstepArch=1;
  int fullDg = 1;
  int space1 = 0;

  //compute solution and BC (given directly to the solver
  // creation of law
  J2LinearDG3DMaterialLaw law1(lawnum,rho,young,nu,1.e100,0.);

  // creation of ElasticField
  int nfield = 101;
  dG3DDomain myfield1(1000,nfield,space1,lawnum,fullDg);
  nonLinearMechSolver mysolver(1000);

  mysolver.loadModel(meshfile);
  mysolver.addDomain(&myfield1);
  mysolver.addMaterialLaw(&law1);
  mysolver.Scheme(soltype);
  mysolver.Solver(sol);
  mysolver.snlData(nstep,ftime,tol);
  mysolver.stepBetweenArchiving(nstepArch);

  mysolver.displacementBC("Face",2376,0,0.0);
  mysolver.displacementBC("Face",2376,1,0.0);
  mysolver.displacementBC("Face",2376,2,0.0);
  mysolver.pressureOnPhysicalGroupBC(1584,1.e6,0.);

  mysolver.internalPointBuildView("svm",IPField::SVM, 1, 1);

  Msg::Info("Init done");
  mysolver.solve();

  Msg::Exit(0);

  return 0;
}
