#include "SmallFem.h"
#include "MPIOStream.h"

#include "Mesh.h"
#include "GroupOfElement.h"

#include "FunctionSpaceVector.h"
#include "FunctionSpaceScalar.h"

#include "FormulationHelper.h"
#include "FormulationDummy.h"
#include "FormulationContainer.h"
#include "FormulationSteadyWave.h"
#include "FormulationSilverMuller.h"

#include "FormulationOSRCVector.h"
#include "FormulationUpdateOSRCVector.h"

#include "SolverDDM.h"
#include "System.h"
#include "SystemHelper.h"

using namespace std;

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

// Peak Memory //
double getPeakMemory(void);

// Compute //
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
  GroupOfElement    infinity(msh.getFromPhysical(1010,   myProc + 1));
  GroupOfElement ddmBoundary(msh.getFromPhysical(40000 + myProc + 1));

  // Full domain
  vector<const GroupOfElement*> domain(5);
  domain[0] = &air;
  domain[1] = &rod;
  domain[2] = &src;
  domain[3] = &infinity;
  domain[4] = &ddmBoundary;

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
  FunctionSpaceVector                OsrcR(ddmBoundaryDomain, orderVol);

  for(int j = 0; j < NPade; j++){
    OsrcPhi[j] = new FunctionSpaceVector(ddmBoundaryDomain, orderVol);

    if(orderVol == 0)
      OsrcRho[j] = new FunctionSpaceScalar(ddmBoundaryDomain, 1);
    else
      OsrcRho[j] = new FunctionSpaceScalar(ddmBoundaryDomain, orderVol);
  }

  // Formulations //
  cout << "FEM terms... " << endl << flush;

  // Waves
  FormulationSteadyWave<Complex>   waveAir(air, fs, k);
  FormulationSteadyWave<Complex>   waveRod(rod, fs, k, nuRRod, epsRRod, sVol);
  FormulationSilverMuller        radiation(infinity, fs, k);

  // Ddm context
  map<Dof, Complex> ddmG;
  map<Dof, Complex> rhsG;
  FormulationHelper::initDofMap(fG, ddmBoundary, ddmG);
  FormulationHelper::initDofMap(fG, ddmBoundary, rhsG);


  DDMContextOSRCVector context(ddmBoundary, dirichlet,
                               fs, fG, OsrcPhi, OsrcRho, OsrcR,
                               k, keps, NPade, Pi / 2);
  context.setDDMDofs(ddmG);

  // Ddm
  FormulationOSRCVector         ddm(context);
  FormulationUpdateOSRCVector upDdm(context);

  // Dummy
  FormulationDummy<Complex> dummy;

  // Container
  FormulationContainer<Complex> allFem;
  allFem.addFormulation(waveAir);
  allFem.addFormulation(waveRod);
  allFem.addFormulation(radiation);

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

    feSol.setSaveMesh(false);
    feSol.setBinaryFormat(true);
    feSol.setParition(myProc + 1);
    feSol.write("boubouchon");
  }

  // Give peak virtual memory //
  double  myVmPeak = getPeakMemory();
  double* alVmPeak = new double[nProcs];

  MPI_Allgather(&myVmPeak,1,MPI_DOUBLE, alVmPeak,1,MPI_DOUBLE, MPI_COMM_WORLD);

  cout << "Peak VM:" << endl << flush;
  for(int i = 0; i < nProcs; i++)
    cout << " ** Process " << i << ": " << alVmPeak[i] << " MB"
         << endl << flush;

  // Give max peak VMem accross all processes //
  int maxIdx = 0;
  for(int i = 1; i < nProcs; i++)
    if(alVmPeak[maxIdx] < alVmPeak[i])
      maxIdx = i;

  cout << "Maximum Peak VM accross MPI processes:" << endl
       << " ** Process " << maxIdx << ": " << alVmPeak[maxIdx] << " MB"
       << endl << flush;

  // Give systems sizes //
  int  mySize = full.getSize();
  int* alSize = new int[nProcs];

  MPI_Allgather(&mySize, 1, MPI_INT, alSize, 1, MPI_INT, MPI_COMM_WORLD);

  cout << "Volume system size:" << endl << flush;
  for(int i = 0; i < nProcs; i++)
    cout << " ** Process " << i << ": " << alSize[i] << " unknowns"
         << endl << flush;

  // Clear //
  for(int j = 0; j < NPade; j++){
    delete OsrcPhi[j];
    delete OsrcRho[j];
  }

  delete[] alVmPeak;
  delete[] alSize;
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

double getPeakMemory(void){
  // Stream //
  ifstream stream("/proc/self/status", ifstream::in);
  char        tmp[1048576];
  double   vmPeak;

  // Is open ? //
  if(!stream.is_open())
    throw Exception("Boubouchons: cannot open /proc/self/status for VmPeak");

  // Look for "VmPeak:" //
  stream >> tmp;
  while(strncmp(tmp, "VmPeak:", 1048576) != 0){
    stream.getline(tmp, 1048576);
    stream >> tmp;
  }

  // Read VmPeak //
  stream >> vmPeak;

  // Close & Return //
  stream.close();
  return vmPeak / 1024;
}
