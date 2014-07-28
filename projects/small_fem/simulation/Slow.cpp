#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationSteadySlow.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

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
  vector<const GroupOfElement*> domain(3);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &wall;

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());

  // Function Space //
  assemble.start();
  FunctionSpaceVector fs(domain, order);

  // Formulation & System //
  FormulationSteadySlow wave(volume, fs, k);
  System<double> sys;
  sys.addFormulation(wave);

  // Dirichlet //
  SystemHelper<double>::dirichlet(sys, fs, source, fSourceVec);
  SystemHelper<double>::dirichlet(sys, fs, wall,   fWallVec);

  cout << "Slow Vectorial "
       << "Steady Wave (Order: " << order
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
    sys.getSolution(feSol, fs, volume);
    feSol.write("slow");
  }

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
