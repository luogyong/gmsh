#include <iostream>
#include <cmath>
#include <complex>

#include "Mesh.h"
#include "System.h"

#include "FormulationSommerfeld.h"
#include "FormulationSteadyWave.h"
#include "FormulationFieldLagrange.h"
#include "FormulationLagrangeField.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

double fSourceReal(fullVector<double>& xyz){
  //return 1;
  //return fabs(xyz(1);
  return exp(-((xyz(1) * 4.2) * (xyz(1) * 4.2) +
               (xyz(2) * 4.2) * (xyz(2) * 4.2)));
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume     = msh.getFromPhysical(7);
  GroupOfElement source     = msh.getFromPhysical(5);
  GroupOfElement freeSpace  = msh.getFromPhysical(6);

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(source);
  domain.add(freeSpace);

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());

  // Formulation //
  assemble.start();
  FunctionSpaceScalar fs(domain, order);

  FormulationSteadyWave<complex<double> > wave(volume, fs, k);
  FormulationSommerfeld                   sommerfeld(freeSpace, fs, k);

  // Lagrange //
  FunctionSpaceScalar lagrange(domain, order);

  FormulationFieldLagrange fieldLagrange(source, fs, lagrange, fSourceReal);
  FormulationLagrangeField lagrangeField(source, lagrange, fs);

  // System //
  System<complex<double> > sys;
  sys.addFormulation(wave);
  sys.addFormulation(sommerfeld);
  sys.addFormulation(fieldLagrange);
  sys.addFormulation(lagrangeField);

  cout << "Free Space Lagrange contrainted (Order: "  << order
       << " --- Wavenumber: "                         << k
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
    FEMSolution<complex<double> > feSolVol;
    FEMSolution<complex<double> > feSolSrc;
    FEMSolution<complex<double> > feSolMul;

    sys.getSolution(feSolVol, fs,       volume);
    sys.getSolution(feSolSrc, fs,       source);
    sys.getSolution(feSolMul, lagrange, source);

    feSolVol.write("lagrVol");
    feSolSrc.write("lagrSrc");
    feSolMul.write("lagrMul");
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
