#include "SmallFem.h"
#include "System.h"
#include "SystemHelper.h"
#include "Interpolator.h"
#include "FormulationHelper.h"
#include "FormulationImpedance.h"
#include "FormulationSilverMuller.h"
#include "FormulationSteadyWave.h"

#include <iostream>

using namespace std;

static const int    scal = 0;
static const int    vect = 1;
static       double k;

Complex fSourceScal(fullVector<double>& xyz){
  const double  Pi = M_PI;
  const double  y  = xyz(1);
  const double  z  = xyz(2);

  const Complex E0 = Complex(1, 0);
  const double  a  = 2;
  const double  b  = 1;
  const int     m  = 1;
  const int     n  = 1;

  const double ky = m * Pi / a;
  const double kz = n * Pi / b;

  return E0 * Complex(sin(ky * y) * sin(kz * z), 0);
}

Complex fZeroScal(fullVector<double>& xyz){
  return Complex(0, 0);
}

fullVector<Complex> fSourceVect(fullVector<double>& xyz){
  const Complex I  = Complex(0, 1);
  const double  Pi = M_PI;
  const double  x  = xyz(0);
  const double  y  = xyz(1);
  const double  z  = xyz(2);

  const Complex E0 = Complex(1, 0);
  const double  a  = 1;
  const double  b  = 1;
  const int     m  = 1;
  const int     n  = 1;

  const Complex ky = Complex(m * Pi / a, 0);
  const Complex kz = Complex(n * Pi / b, 0);
  const Complex kx = sqrt(Complex(k * k, 0) - (ky * ky) - (kz * kz));
  /*
  // TEM 2D
  fullVector<Complex> tmp(3);
  tmp(0) = Complex(0, 0);
  tmp(1) = E0;
  tmp(2) = Complex(0, 0);
  */
  /*
  // TMm 2D
  fullVector<Complex> tmp(3);
  tmp(0) = E0 * I * ky / k * sin(ky * y);
  tmp(1) = E0 *     kx / k * cos(ky * y);
  tmp(2) = Complex(0, 0);
  */
  /*
  // TEm0 3D
  fullVector<Complex> tmp(3);
  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(0, 0);
  tmp(2) = E0 * sin(ky * y);
  */
  /*
  // TEmn 3D
  fullVector<Complex> tmp(3);
  tmp(0) = Complex(0, 0);
  tmp(1) = -E0 * cos(ky * y) * sin(kz * z);
  tmp(2) = +E0 * sin(ky * y) * cos(kz * z);
  */
  // TMmn 3D
  fullVector<Complex> tmp(3);
  tmp(0) = E0                                  * sin(ky * y) * sin(kz * z);
  tmp(1) = E0 * (-I * kx * ky) / (k*k - kx*kx) * cos(ky * y) * sin(kz * z);
  tmp(2) = E0 * (-I * kx * kz) / (k*k - kx*kx) * sin(ky * y) * cos(kz * z);

  return tmp;
}

fullVector<Complex> fZeroVect(fullVector<double>& xyz){
  fullVector<Complex> tmp(3);

  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(0, 0);
  tmp(2) = Complex(0, 0);

  return tmp;
}

void compute(const Options& option){
  // Get Type //
  int type;
  if(option.getValue("-type")[1].compare("scalar") == 0){
    cout << "Scalar Waveguide" << endl << flush;
    type = scal;
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    cout << "Vetorial Waveguide" << endl << flush;
    type = vect;
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Get Parameters //
  const size_t nDom  = atoi(option.getValue("-n")[1].c_str());
  k                  = atof(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  const double sigma = atof(option.getValue("-sigma")[1].c_str());

  cout << "Wavenumber: " << k     << endl
       << "Order:      " << order << endl
       << "# Domain:   " << nDom  << endl << flush;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement   source(msh.getFromPhysical(1 * nDom + 1));
  GroupOfElement     zero(msh.getFromPhysical(2 * nDom + 2));
  GroupOfElement infinity(msh.getFromPhysical(2 * nDom + 1));
  GroupOfElement   volume(msh);

  vector<GroupOfElement*> perVolume(nDom);
  for(size_t i = 0; i < nDom; i++){
    perVolume[i] = new GroupOfElement(msh.getFromPhysical(i + 1));
    volume.add(*perVolume[i]);
  }

  // Full Domain //
  vector<const GroupOfElement*> domain(4);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &zero;
  domain[3] = &infinity;

  // Function Space //
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Steady Wave Formulation //
  const double    Z0 = 119.9169832 * M_PI;
  const Complex epsr(1, 1 / k * Z0 * sigma);
  const Complex  mur(1, 0);

  FormulationSteadyWave<Complex>  wave(volume,   *fs, k);
  FormulationImpedance       impedance(infinity, *fs, k, epsr, mur);

  // Solve //
  System<Complex> system;
  system.addFormulation(wave);
  system.addFormulation(impedance);

  // Constraint
  if(fs->isScalar()){
    SystemHelper<Complex>::dirichlet(system, *fs, zero  , fZeroScal);
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceScal);
  }
  else{
    SystemHelper<Complex>::dirichlet(system, *fs, zero  , fZeroVect);
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceVect);
  }

  // Assemble
  system.assemble();
  cout << "Assembled: " << system.getSize() << endl << flush;

  // Sove
  system.solve();
  cout << "Solved!" << endl << flush;

  // Draw Solution //
  try{
    option.getValue("-nopos");
  }

  catch(...){
    cout << "Writing solution..." << endl << flush;

    stringstream stream;
    try{
      vector<string> name = option.getValue("-name");
      stream << name[1];
    }
    catch(...){
      stream << "waveguide";
    }

    try{
      // Get Visu Mesh //
      vector<string> visuStr = option.getValue("-interp");
      Mesh           visuMsh(visuStr[1]);

      // Get Solution //
      map<Dof, Complex> sol;

      FormulationHelper::initDofMap(*fs, volume, sol);
      system.getSolution(sol, 0);

      // Interoplate //
      for(size_t i = 0; i < nDom; i++){
        // GroupOfElement to interoplate on
        GroupOfElement visuGoe(visuMsh.getFromPhysical(i + 1));

        // Interpolation
        stringstream name;
        name << stream.str() << i << ".dat";

        map<const MVertex*, vector<Complex> > map;
        Interpolator<Complex>::interpolate(*perVolume[i], visuGoe, *fs,sol,map);
        Interpolator<Complex>::write(name.str(), map);
      }
    }

    catch(...){
      FEMSolution<Complex> feSol;
      system.getSolution(feSol, *fs, volume);
      feSol.write(stream.str());
    }
  }

  // Clean //
  for(size_t i = 0; i < nDom; i++)
    delete perVolume[i];

  delete fs;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-n,-sigma,-type,-interp,-name,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
