#include "SmallFem.h"
#include "BoubouchonsHelper.h"

#include "Mesh.h"
#include "GroupOfElement.h"

#include "FunctionSpaceVector.h"
#include "FormulationSteadyWave.h"
#include "System.h"
#include "SystemHelper.h"

using namespace std;

typedef FormulationSteadyWave<Complex> Wave;

const double Pi        = 4 * atan(1);
const double C0        = 299792458;
const double Mu0       = 4 * Pi * 1e-7;

const double epsRRodRe = 6;
const double epsRRodIm = 0;
const double srcL      = 0.005;

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

  GroupOfElement     air(msh.getFromPhysical(1007));
  GroupOfElement     rod(msh.getFromPhysical(1008));
  GroupOfElement  srcVol(msh.getFromPhysical(1009));
  GroupOfElement srcLine(msh.getFromPhysical(1011));

  GroupOfElement pmlXYZ(msh.getFromPhysical(1000));
  GroupOfElement pmlXZ(msh.getFromPhysical(1001));
  GroupOfElement pmlYZ(msh.getFromPhysical(1002));
  GroupOfElement pmlXY(msh.getFromPhysical(1003));
  GroupOfElement pmlZ(msh.getFromPhysical(1004));
  GroupOfElement pmlY(msh.getFromPhysical(1005));
  GroupOfElement pmlX(msh.getFromPhysical(1006));

  // Full Domain
  vector<const GroupOfElement*> dom(11);
  dom[0]  = &air;
  dom[1]  = &rod;
  dom[2]  = &srcVol;
  dom[3]  = &srcLine;
  dom[4]  = &pmlXYZ;
  dom[5]  = &pmlXZ;
  dom[6]  = &pmlYZ;
  dom[7]  = &pmlXY;
  dom[8]  = &pmlZ;
  dom[9]  = &pmlY;
  dom[10] = &pmlX;

  // Wavenumber //
  double f  = atof(option.getValue("-f")[1].c_str());
  Omega0    = 2 * Pi * f;
  double  k = Omega0 / C0;

  // FunctionSpace //
  cout << "Functionspace... " << endl << flush;
  size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(dom, order);

  // Formulations for wave //
  cout << "FEM terms... " << endl << flush;
  Wave    waveAir(air,    fs, k);
  Wave waveSrcVol(srcVol, fs, k);
  Wave    waveRod(rod,    fs, k, nuRRod, epsRRod, sVol);

  Wave wavePmlXYZ(pmlXYZ, fs, k, Material::XYZ::Nu, Material::XYZ::Eps, sVol);
  Wave  wavePmlXY(pmlXY,  fs, k,  Material::XY::Nu,  Material::XY::Eps, sVol);
  Wave  wavePmlYZ(pmlYZ,  fs, k,  Material::YZ::Nu,  Material::YZ::Eps, sVol);
  Wave  wavePmlXZ(pmlXZ,  fs, k,  Material::XZ::Nu,  Material::XZ::Eps, sVol);
  Wave   wavePmlX(pmlX,   fs, k,   Material::X::Nu,   Material::X::Eps, sVol);
  Wave   wavePmlY(pmlY,   fs, k,   Material::Y::Nu,   Material::Y::Eps, sVol);
  Wave   wavePmlZ(pmlZ,   fs, k,   Material::Z::Nu,   Material::Z::Eps, sVol);

  // System //
  System<Complex> system;
  system.addFormulation(waveAir);
  system.addFormulation(waveSrcVol);
  system.addFormulation(waveRod);
  system.addFormulation(wavePmlXYZ);
  system.addFormulation(wavePmlXY);
  system.addFormulation(wavePmlYZ);
  system.addFormulation(wavePmlXZ);
  system.addFormulation(wavePmlX);
  system.addFormulation(wavePmlY);
  system.addFormulation(wavePmlZ);

  // Dirichlet //
  cout << "Dirichlet conditions... " << endl << flush;
  SystemHelper<Complex>::dirichlet(system, fs, srcLine, fSrc);

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
    system.getSolution(feSol, fs, dom);
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
  ret(2) = Complex(0, 1) * Omega0 * Mu0 * cos(Pi * xyz(2) / srcL);

  return ret;
}

fullVector<Complex> sVol(fullVector<double>& xyz){
  fullVector<Complex> ret(3);
  ret(0) = 0; ret(1) = 0; ret(2) = 0;

  return ret;
}
