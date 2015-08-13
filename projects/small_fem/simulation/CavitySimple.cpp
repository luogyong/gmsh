#include "SmallFem.h"
#include "Timer.h"
#include "System.h"
#include "SystemHelper.h"
#include "Interpolator.h"
#include "FormulationHelper.h"
#include "FormulationSilverMuller.h"
#include "FormulationSteadyWave.h"
#include "FormulationSteadySlow.h"

#include <iostream>

using namespace std;

static const int   scal = 0;
static const int   vect = 1;

static const double  Pi = M_PI;

static const double  E0 = 1;
static const double  a  = 0.25;
static const double  b  = 0.25;
static const int     m  = 1;
static const int     n  = 1;

static const double  ky = m * Pi / a;
static const double  kz = n * Pi / b;
static       double  k;

////////////////////
// Sources Fileds //
////////////////////
double             fSourceScal(fullVector<double>& xyz);
double               fZeroScal(fullVector<double>& xyz);
fullVector<double> fSourceVect(fullVector<double>& xyz);
fullVector<double>   fZeroVect(fullVector<double>& xyz);

double fSourceScal(fullVector<double>& xyz){
  const double y = xyz(1);
  const double z = xyz(2);

  return E0 * sin(ky * y) * sin(kz * z);
}

double fZeroScal(fullVector<double>& xyz){
  return 0;
}

fullVector<double> fSourceVect(fullVector<double>& xyz){
  const double y  = xyz(1);
  const double z  = xyz(2);

  fullVector<double> tmp(3);

  // TEmn 3D
  tmp(0) = 0;
  tmp(1) = -E0 * cos(ky * y) * sin(kz * z);
  tmp(2) = +E0 * sin(ky * y) * cos(kz * z);

  return tmp;
}

fullVector<double> fZeroVect(fullVector<double>& xyz){
  fullVector<double> tmp(3);

  tmp(0) = 0;
  tmp(1) = 0;
  tmp(2) = 0;

  return tmp;
}

/////////////
// Compute //
/////////////
void compute(const Options& option){
  // Timers
  Timer tAssembly;
  Timer tSolve;

  // Get Type //
  int type;
  if(option.getValue("-type")[1].compare("scalar") == 0){
    cout << "Scalar Cavity" << endl << flush;
    type = scal;
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    cout << "Vetorial Cavity" << endl << flush;
    type = vect;
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  k                  = atof(option.getValue("-k")[1].c_str());

  cout << "Wavenumber: " << k     << endl
       << "Order:      " << order << endl << flush;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement source(msh.getFromPhysical(5));
  GroupOfElement volume(msh.getFromPhysical(7));
  GroupOfElement   zero(msh);

  zero.add(msh.getFromPhysical(4)); // 'Infinity'
  zero.add(msh.getFromPhysical(6)); // 'Walls'


  // Full Domain //
  vector<const GroupOfElement*> domain(3);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &zero;

  // Function Space //
  tAssembly.start();
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Steady Wave Formulation //
  Formulation<double>* wave;

  try{
    option.getValue("-slow");
    cout << "Slow Formulation" << endl << flush;
    wave = new FormulationSteadySlow<double>(volume, *fs, k);
  }

  catch(...){
    cout << "Fast Formulation" << endl << flush;
    wave = new FormulationSteadyWave<double>(volume, *fs, k);
  }

  // Solve //
  System<double> system;
  system.addFormulation(*wave);

  // Constraint
  if(fs->isScalar()){
    SystemHelper<double>::dirichlet(system, *fs, zero  , fZeroScal);
    SystemHelper<double>::dirichlet(system, *fs, source, fSourceScal);
  }
  else{
    SystemHelper<double>::dirichlet(system, *fs, zero  , fZeroVect);
    SystemHelper<double>::dirichlet(system, *fs, source, fSourceVect);
  }

  // Assemble
  system.assemble();
  tAssembly.stop();
  cout << "Assembled: " << system.getSize()
       << "! (" << tAssembly.time() << " " << tAssembly.unit() << ")"
       << endl << flush;

  // Sove
  tSolve.start();
  system.solve();
  tSolve.stop();
  cout << "Solved!"
       << " (" << tSolve.time() << " " << tSolve.unit() << ")"
       << endl << flush;

  // Draw Solution //
  try{
    option.getValue("-nopos");
  }

  catch(...){
    cout << "Writing solution..." << endl << flush;

    FEMSolution<double> feSol;
    system.getSolution(feSol, *fs, volume);
    feSol.write("cavity");
  }

  // Clean //
  delete wave;
  delete fs;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-type,-slow,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
