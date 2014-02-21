#include <iostream>
#include <complex>
#include <cmath>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationSteadyWave.h"
#include "FormulationPML.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

complex<double> fSourceScal(fullVector<double>& xyz){
  return complex<double>(1, 0);
  //return complex<double>(fabs(xyz(1)), 0);
  //return complex<double>(exp(-((xyz(1) * 4.2) * (xyz(1) * 4.2) +
  //                             (xyz(2) * 4.2) * (xyz(2) * 4.2))), 0);
}

complex<double> fInfinityScal(fullVector<double>& xyz){
  return complex<double>(0, 0);
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume     = msh.getFromPhysical(7);
  GroupOfElement outerSpace = msh.getFromPhysical(8);
  GroupOfElement source     = msh.getFromPhysical(5);
  GroupOfElement infinity   = msh.getFromPhysical(4);

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(outerSpace);
  domain.add(source);

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());

  // Formulation //
  assemble.start();
  FunctionSpaceScalar fs(domain, order);

  FormulationSteadyWave<complex<double> > wave(volume, fs, k);
  FormulationPML pml(outerSpace, fs, k);

  // System //
  System<complex<double> > sys;
  sys.addFormulation(wave);
  sys.addFormulation(pml);

  SystemHelper<complex<double> >::dirichlet(sys, fs, source,   fSourceScal);
  SystemHelper<complex<double> >::dirichlet(sys, fs, infinity, fInfinityScal);

  cout << "Free Space (Order: "  << order
       << " --- Wavenumber: "    << k
       << "): " << sys.getSize() << endl;

  // Assemble //
  sys.assemble();
  assemble.stop();

  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  // Solve //
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
    FEMSolution<complex<double> > feSol;
    sys.getSolution(feSol);
    feSol.write("free");
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
