#include "SmallFem.h"
#include "MPIOStream.h"
#include "BoubouchonsHelper.h"

#include "Mesh.h"
#include "GroupOfElement.h"

#include "FunctionSpaceVector.h"
#include "FunctionSpaceScalar.h"

#include "FormulationHelper.h"
#include "FormulationDummy.h"
#include "FormulationContainer.h"
#include "FormulationSteadyWave.h"

#include "FormulationEMDA.h"
#include "FormulationOSRCVector.h"
#include "FormulationUpdateEMDA.h"
#include "FormulationUpdateOSRCVector.h"

#include "SolverDDM.h"
#include "System.h"
#include "SystemHelper.h"


using namespace std;

typedef FormulationSteadyWave<Complex> Wave;

const double Pi        = 4 * atan(1);
const double C0        = 299792458;
const double Mu0       = 4 * Pi * 1e-7;

const double epsRRodRe = 6;
const double epsRRodIm = 0;

double Omega0;

// Rods //
void epsRRod(fullVector<double>& xyz, fullMatrix<Complex>& epsR);
void  nuRRod(fullVector<double>& xyz, fullMatrix<Complex>& nuR);

// Source //
fullVector<Complex> fSrc(fullVector<double>& xyz);

// Dummy volume source term //
fullVector<Complex> sVol(fullVector<double>& xyz);

void compute(const Options& option){
  // MPI //
  int nProcs;
  int myProc;
  MPIOStream cout(0, std::cout);

  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myProc);

  // Get parameters //
  // Wavenumber
  double f  = atof(option.getValue("-f")[1].c_str());
  Omega0    = 2 * Pi * f;
  double  k = Omega0 / C0;

  // FEM orders
  const int orderVol = atoi(option.getValue("-ov")[1].c_str());
  const int orderSur = atoi(option.getValue("-ob")[1].c_str());

  // EMDA
  //double chi = 0;

  // OSRC
  double    ck = 0;
  int    NPade = 4;
  Complex keps = k + Complex(0, k * ck);

  // Get Domains //
  cout << "Reading domain... " << endl << flush;
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement         air(msh.getFromPhysical(1007,   myProc + 1));
  GroupOfElement         rod(msh.getFromPhysical(1008,   myProc + 1));
  GroupOfElement         src(msh.getFromPhysical(1009,   myProc + 1));
  GroupOfElement      pmlXYZ(msh.getFromPhysical(1000,   myProc + 1));
  GroupOfElement       pmlXZ(msh.getFromPhysical(1001,   myProc + 1));
  GroupOfElement       pmlYZ(msh.getFromPhysical(1002,   myProc + 1));
  GroupOfElement       pmlXY(msh.getFromPhysical(1003,   myProc + 1));
  GroupOfElement        pmlZ(msh.getFromPhysical(1004,   myProc + 1));
  GroupOfElement        pmlY(msh.getFromPhysical(1005,   myProc + 1));
  GroupOfElement        pmlX(msh.getFromPhysical(1006,   myProc + 1));
  GroupOfElement ddmBoundary(msh.getFromPhysical(40000 + myProc + 1));

  // Full domain
  vector<const GroupOfElement*> domain(11);
  domain[0]  = &air;
  domain[1]  = &rod;
  domain[2]  = &src;
  domain[3]  = &pmlXYZ;
  domain[4]  = &pmlXZ;
  domain[5]  = &pmlYZ;
  domain[6]  = &pmlXY;
  domain[7]  = &pmlZ;
  domain[8]  = &pmlY;
  domain[9]  = &pmlX;
  domain[10] = &ddmBoundary;

  // Dirichlet boundary container
  vector<const GroupOfElement*> dirichlet(1);
  dirichlet[0] = &src;

  // DDM boundary container
  vector<const GroupOfElement*> ddmBoundaryDomain(1);
  ddmBoundaryDomain[0] = &ddmBoundary;

  // Function space //
  FunctionSpaceVector fs(           domain, orderVol);
  FunctionSpaceVector fG(ddmBoundaryDomain, orderSur);

  vector<const FunctionSpaceVector*> OsrcPhi(NPade, NULL);
  vector<const FunctionSpaceScalar*> OsrcRho(NPade, NULL);
  FunctionSpaceVector                OsrcR(ddmBoundaryDomain, orderSur);

  for(int j = 0; j < NPade; j++){
    OsrcPhi[j] = new FunctionSpaceVector(ddmBoundaryDomain, orderSur);

    if(orderSur == 0)
      OsrcRho[j] = new FunctionSpaceScalar(ddmBoundaryDomain, 1);
    else
      OsrcRho[j] = new FunctionSpaceScalar(ddmBoundaryDomain, orderSur);
  }

  // Formulations //
  cout << "FEM terms... " << endl << flush;

  // Waves
  Wave    waveAir(air,    fs, k);
  Wave    waveRod(rod,    fs, k, nuRRod, epsRRod, sVol);
  Wave wavePmlXYZ(pmlXYZ, fs, k, Material::XYZ::Nu, Material::XYZ::Eps, sVol);
  Wave  wavePmlXY(pmlXY,  fs, k,  Material::XY::Nu,  Material::XY::Eps, sVol);
  Wave  wavePmlYZ(pmlYZ,  fs, k,  Material::YZ::Nu,  Material::YZ::Eps, sVol);
  Wave  wavePmlXZ(pmlXZ,  fs, k,  Material::XZ::Nu,  Material::XZ::Eps, sVol);
  Wave   wavePmlX(pmlX,   fs, k,   Material::X::Nu,   Material::X::Eps, sVol);
  Wave   wavePmlY(pmlY,   fs, k,   Material::Y::Nu,   Material::Y::Eps, sVol);
  Wave   wavePmlZ(pmlZ,   fs, k,   Material::Z::Nu,   Material::Z::Eps, sVol);

  // Ddm context
  map<Dof, Complex> ddmG;
  map<Dof, Complex> rhsG;
  FormulationHelper::initDofMap(fG, ddmBoundary, ddmG);
  FormulationHelper::initDofMap(fG, ddmBoundary, rhsG);


  DDMContextOSRCVector context(ddmBoundary, dirichlet,
                               fs, fG, OsrcPhi, OsrcRho, OsrcR,
                               k, keps, NPade, Pi / 2);
  //DDMContextEMDA context(ddmBoundary, dirichlet, fs, fG, k, chi);
  context.setDDMDofs(ddmG);

  // Ddm
  FormulationOSRCVector         ddm(context);
  FormulationUpdateOSRCVector upDdm(context);
  //FormulationEMDA         ddm(context);
  //FormulationUpdateEMDA upDdm(context);

  // Dummy
  FormulationDummy<Complex> dummy;

  // Container
  FormulationContainer<Complex> allFem;
  allFem.addFormulation(waveAir);
  allFem.addFormulation(waveRod);
  allFem.addFormulation(wavePmlXYZ);
  allFem.addFormulation(wavePmlXY);
  allFem.addFormulation(wavePmlYZ);
  allFem.addFormulation(wavePmlXZ);
  allFem.addFormulation(wavePmlX);
  allFem.addFormulation(wavePmlY);
  allFem.addFormulation(wavePmlZ);

  // Solve non-homogenous problem //
  cout << "Solving non homogenous problem" << endl << flush;

  System<Complex>* nonHomogenous = new System<Complex>;
  nonHomogenous->addFormulation(allFem);
  nonHomogenous->addFormulation(ddm);

  SystemHelper<Complex>::dirichlet(*nonHomogenous, fs, src, fSrc);

  nonHomogenous->assemble();
  nonHomogenous->solve();

  // Solve non-homogenous DDM problem //
  cout << "Computing right hand side" << endl << flush;

  context.setSystem(*nonHomogenous);
  upDdm.update(); // update volume solution (at DDM border)

  System<Complex>* nonHomogenousDDM = new System<Complex>;
  nonHomogenousDDM->addFormulation(upDdm);

  nonHomogenousDDM->assemble();
  nonHomogenousDDM->solve();
  nonHomogenousDDM->getSolution(rhsG, 0);

  // Clear Systems //
  delete nonHomogenous;
  delete nonHomogenousDDM;

  // DDM Solver //
  cout << "Solving DDM problem" << endl << flush;
  SolverDDM* solver = new SolverDDM(allFem, dummy, context, ddm, upDdm, rhsG);

  // Solve
  int maxIt = 1000;
  solver->setMaximumIteration(maxIt);
  solver->setRestart(maxIt); // No restart!
  cout << " ! Warning: no restart ! " << endl;
  solver->solve();

  // Get Solution
  solver->getSolution(ddmG);
  context.setDDMDofs(ddmG);

  // Clear DDM
  delete solver;

  // Full Problem //
  cout << "Solving full problem" << endl << flush;
  ddm.update();

  System<Complex> full;
  full.addFormulation(allFem);
  full.addFormulation(ddm);

  SystemHelper<Complex>::dirichlet(full, fs, src, fSrc);

  full.assemble();
  full.solve();

  // Draw Solution //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    cout << "Writing full problem" << endl << flush;

    FEMSolution<Complex> feSol;
    full.getSolution(feSol, fs, domain);

    feSol.setSaveMesh(true);
    feSol.setBinaryFormat(true);
    feSol.setParition(myProc + 1);
    feSol.write("boubouchon");
  }

  // Clear //
  for(int j = 0; j < NPade; j++){
    delete OsrcPhi[j];
    delete OsrcRho[j];
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-ov,-ob,-f,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}

void epsRRod(fullVector<double>& xyz, fullMatrix<Complex>& epsR){
  epsR.scale(0);

  epsR(0, 0) = Complex(epsRRodRe, epsRRodIm);
  epsR(1, 1) = Complex(epsRRodRe, epsRRodIm);
  epsR(2, 2) = Complex(epsRRodRe, epsRRodIm);
}

void nuRRod(fullVector<double>& xyz, fullMatrix<Complex>& nuR){
  nuR.scale(0);

  nuR(0, 0) = 1;
  nuR(1, 1) = 1;
  nuR(2, 2) = 1;
}

fullVector<Complex> fSrc(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Complex(0, 0);
  ret(1) = Complex(0, 0);
  ret(2) = Complex(0, 1) * Omega0 * Mu0;

  return ret;
}

fullVector<Complex> sVol(fullVector<double>& xyz){
  fullVector<Complex> ret(3);
  ret(0) = 0; ret(1) = 0; ret(2) = 0;

  return ret;
}
