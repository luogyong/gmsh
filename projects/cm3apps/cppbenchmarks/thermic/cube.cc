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
  int lawnum=1;
  double rho   = 7850.;
  double young = 2.1e11;
  double nu = 0.3;
  double MUxy = young/(2.*(1.+nu));
  double cv=1.;
  double t0=0.;
  double Kx=51.9;
  double alphax=12.e-6;


  // geometry
  std::string meshfile("cube.msh");

  // solver
  int sol = 2;
  int soltype = 1;
  int nstep = 100;
  double ftime =1.;
  double tol=1.e-6;
  int nstepArch=1;
  int fullDg = 1;
  int space1 = 0;
  double beta1  = 10;

  double alpha =0.;
  double cp = 0.;

  //  compute solution and BC (given directly to the solver
  // creation of law

  LinearThermoMechanicsDG3DMaterialLaw law1(lawnum,rho,young,young,young,nu,nu,nu,MUxy,MUxy,MUxy,alpha,alpha,alpha,cv,t0,Kx,Kx,Kx,alphax,alphax,alphax,cp);

  int nfield = 10; // number of the field (physical number of surface)
  ThermoMechanicsDG3DDomain myfield1(1000,nfield,space1,lawnum,fullDg,1.e6);
  myfield1.matrixByPerturbation(0,0,0,1e-8);
  myfield1.stabilityParameters(beta1);
  myfield1.ThermoMechanicsStabilityParameters(beta1,1);
  // creation of Solver
  nonLinearMechSolver mysolver(1000);
  mysolver.loadModel(meshfile);
  mysolver.addDomain(&myfield1);
  mysolver.addMaterialLaw(&law1);
  mysolver.Scheme(soltype);
  mysolver.Solver(sol);
  mysolver.snlData(nstep,ftime,tol);
  mysolver.stepBetweenArchiving(nstepArch);
  // BC

  //mechanical BC
  mysolver.displacementBC("Face",1234,2,0.);
  mysolver.displacementBC("Face",2376,0,0.);
  mysolver.displacementBC("Face",1265,1,0.);

  //thermal BC
  //mysolver.displacementBC("Face",1234,3,200.);
  //mysolver.displacementBC("Face",5678,3,100.);
  mysolver.forceBC("Flux",1234,3,5.19e6);

  mysolver.internalPointBuildView("svm",IPField::SVM, 1, 1);
  mysolver.internalPointBuildView("sig_xx",IPField::SIG_XX, 1, 1);
  mysolver.internalPointBuildView("sig_yy",IPField::SIG_YY, 1, 1);
  mysolver.internalPointBuildView("sig_zz",IPField::SIG_ZZ, 1, 1);
  mysolver.internalPointBuildView("sig_xy",IPField::SIG_XY, 1, 1);
  mysolver.internalPointBuildView("sig_yz",IPField::SIG_YZ, 1, 1);
  mysolver.internalPointBuildView("sig_xz",IPField::SIG_XZ, 1, 1);
  mysolver.internalPointBuildView("temperature",IPField::TEMPERATURE, 1, 1);
  mysolver.internalPointBuildView("qx",IPField::THERMALFLUX_X, 1, 1);
  mysolver.internalPointBuildView("qy",IPField::THERMALFLUX_Y, 1, 1);
  mysolver.internalPointBuildView("qz",IPField::THERMALFLUX_Z, 1, 1);
  mysolver.archivingForceOnPhysicalGroup("Face", 1234, 2);
  mysolver.archivingForceOnPhysicalGroup("Face", 5678, 2);
  mysolver.archivingForceOnPhysicalGroup("Face", 1234, 3);
  mysolver.archivingForceOnPhysicalGroup("Face", 5678, 3);

  mysolver.solve();

  Msg::Exit(0);
  return 0;
}
