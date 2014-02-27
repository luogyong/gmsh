#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationSteadyWave.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

static const size_t scal = 0;
static const size_t vect = 1;

fullVector<double> fSourceVec(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 1;
  res(2) = 0;

  return res;
}

fullVector<double> fWallVec(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 0;
  res(2) = 0;

  return res;
}

double fSourceScal(fullVector<double>& xyz){
  return 1;
}

double fWallScal(fullVector<double>& xyz){
  return 0;
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume = msh.getFromPhysical(7);
  GroupOfElement source = msh.getFromPhysical(5);
  GroupOfElement wall   = msh.getFromPhysical(6);

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(source);
  domain.add(wall);

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());

  // Get Type //
  size_t type;

  if(option.getValue("-type")[1].compare("scalar") == 0){
    type = scal;
    cout << "Scalar ";
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    type = vect;
    cout << "Vectorial ";
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Function Space //
  assemble.start();
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Formulation & System //
  FormulationSteadyWave<double> wave(volume, *fs, k);
  System<double> sys;
  sys.addFormulation(wave);

  // Dirichlet //
  if(type == scal){
    SystemHelper<double>::dirichlet(sys, *fs, source, fSourceScal);
    SystemHelper<double>::dirichlet(sys, *fs, wall,   fWallScal);
  }

  else{
    SystemHelper<double>::dirichlet(sys, *fs, source, fSourceVec);
    SystemHelper<double>::dirichlet(sys, *fs, wall,   fWallVec);
  }

  cout << "Steady Wave (Order: " << order
       << " --- Wavenumber: "    << k
       << "): " << sys.getSize() << endl;

  // Assemble and solve //
  sys.assemble();
  assemble.stop();
  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  solve.start();
  sys.solve();
  solve.stop();
  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<double> feSol;
    sys.getSolution(feSol, *fs, volume);
    feSol.write("swave");
  }

  // Clean //
  delete fs;

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-nopos,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
