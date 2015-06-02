#include "SmallFem.h"

#include "Mesh.h"
#include "GroupOfElement.h"

#include "FunctionSpaceVector.h"
#include "FormulationSteadyWave.h"
#include "FormulationSilverMuller.h"

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

void compute(const Options& option){
  // Get Domains //
  cout << "Reading domain... " << endl << flush;
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement      air(msh.getFromPhysical(1007));
  GroupOfElement      rod(msh.getFromPhysical(1008));
  GroupOfElement      src(msh.getFromPhysical(1009));
  GroupOfElement infinity(msh.getFromPhysical(1010));

  // Full Domain
  vector<const GroupOfElement*> domain(4);
  domain[0] = &air;
  domain[1] = &rod;
  domain[2] = &src;
  domain[3] = &infinity;

  // Wavenumber //
  double f  = atof(option.getValue("-f")[1].c_str());
  Omega0    = 2 * Pi * f;
  double  k = Omega0 / C0;

  // FunctionSpace //
  cout << "Functionspace... " << endl << flush;
  size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(domain, order);

  // Formulations for wave //
  cout << "FEM terms... " << endl << flush;
  FormulationSteadyWave<Complex>   waveAir(air, fs, k);
  FormulationSteadyWave<Complex>   waveRod(rod, fs, k, nuRRod, epsRRod, sVol);
  FormulationSilverMuller        radiation(infinity, fs, k);

  // System //
  System<Complex> system;
  system.addFormulation(waveAir);
  system.addFormulation(waveRod);
  system.addFormulation(radiation);

  // Dirichlet //
  cout << "Dirichlet conditions... " << endl << flush;
  SystemHelper<Complex>::dirichlet(system, fs, src, fSrc);

  // Assemble and solve //
  cout << "Boubouchons (Order: " << order
       << " --- Frequency: "     << f << ")" << endl;

  system.assemble();
  cout << "Assembled! " << system.getSize() << endl << flush;

  system.solve();
  cout << "Solved! " << endl << flush;

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<Complex> feSol;
    system.getSolution(feSol, fs, domain);
    feSol.write("boubouchons");
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-f,-nopos");
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
