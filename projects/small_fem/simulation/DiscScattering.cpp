#include <iostream>

#include "SmallFem.h"

#include "System.h"
#include "SystemHelper.h"
#include "Interpolator.h"

#include "FormulationHelper.h"
#include "FormulationSommerfeld.h"
#include "FormulationSteadyWave.h"

using namespace std;

static const int scal = 0;
static const int vect = 1;

static double k; // Need to be more sexy !
static double theta = 0;

Complex fSourceScal(fullVector<double>& xyz){
  double p = xyz(0) * cos(theta) + xyz(1) * sin(theta);

  return Complex (1, 0) * Complex(cos(k * p), sin(k * p));
}

fullVector<Complex> fSourceVect(fullVector<double>& xyz){
  double p = xyz(0) * cos(theta) + xyz(1) * sin(theta);

  fullVector<Complex> tmp(3);

  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(1, 0);// * Complex(cos(k * p), sin(k * p));
  tmp(2) = Complex(0, 0);

  return tmp;
}

void compute(const Options& option){
  // Get Type //
  int type;
  if(option.getValue("-type")[1].compare("scalar") == 0){
    cout << "Scalar Disc Scattering" << endl << flush;
    type = scal;
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    cout << "Vectorial Disc Scattering" << endl << flush;
    type = vect;
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Get Parameters //
  const size_t nDom  = atoi(option.getValue("-n")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  k                  = atof(option.getValue("-k")[1].c_str());

  cout << "Wavenumber: " << k     << endl
       << "Order:      " << order << endl
       << "# Domain:   " << nDom  << endl << flush;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume(msh);
  GroupOfElement source(msh);
  GroupOfElement infinity(msh);

  vector<GroupOfElement*> perVolume(nDom);

  // Source
  source.add(msh.getFromPhysical(1000));

  // Infinity
  infinity.add(msh.getFromPhysical(2000 + nDom - 1));

  // Volume
  for(size_t i = 0; i < nDom; i++){
    perVolume[i] = new GroupOfElement(msh.getFromPhysical(100 + i));
    volume.add(*perVolume[i]);
  }

  // Full Domain //
  vector<const GroupOfElement*> domain(3);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &infinity;

  // Dirichlet Border //
  vector<const GroupOfElement*> dirichlet(1);
  dirichlet[0] = &source;

  // Function Space //
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Steady Wave Formulation //
  FormulationSteadyWave<Complex> wave(volume, *fs, k);
  FormulationSommerfeld  sommerfeld(infinity, *fs, k);

  // Solve //
  System<Complex> system;
  system.addFormulation(wave);
  system.addFormulation(sommerfeld);

  // Constraint
  if(fs->isScalar())
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceScal);
  else
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceVect);

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
    stream << "disc";

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
        GroupOfElement visuGoe(visuMsh.getFromPhysical(100 + i));

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
  SmallFem::Keywords("-msh,-o,-k,-type,-n,-interp,-nopos");

  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
