#include <complex>
#include <iostream>

#include "SmallFem.h"

#include "Timer.h"

#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "BasisGenerator.h"
#include "TriLagrangeBasis.h"
#include "LineNodeBasis.h"
#include "LineEdgeBasis.h"
#include "LineNedelecBasis.h"
#include "TriNodeBasis.h"
#include "QuadNedelecBasis.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationNeumann.h"
#include "FormulationEMDA.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"

#include "FormulationSteadyWave.h"
#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;
/*
complex<double> fDirichlet0(fullVector<double>& xyz){
  return complex<double>(0, 0);
}

complex<double> fDirichlet1(fullVector<double>& xyz){
  return complex<double>(1, 0);
}
*/

fullVector<double> fDirichlet0(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 0;
  res(2) = 0;

  return res;
}

fullVector<double> fDirichlet1(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 1;
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
  GroupOfElement    volume = msh.getFromPhysical(7);
  GroupOfElement boundary0 = msh.getFromPhysical(6);
  GroupOfElement boundary1 = msh.getFromPhysical(5);

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(boundary0);
  domain.add(boundary1);

  // Get Order //
  size_t order = atoi(option.getValue("-o")[1].c_str());

  // Function Space //
  assemble.start();
  FunctionSpaceVector fs(domain, order);

  // Compute //
  FormulationSteadyWave<double> wave(volume, fs, 5);

  System<double> sys;
  sys.addFormulation(wave);

  SystemHelper<double>::dirichlet(sys, fs, boundary0, fDirichlet0);
  SystemHelper<double>::dirichlet(sys, fs, boundary1, fDirichlet1);

  cout << "Test -- Order " << order
       << ": " << sys.getSize()
       << endl << flush;

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
    sys.getSolution(feSol);
    feSol.write("test");
  }

  timer.stop();
  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;

}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("-msh,-o,-k,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
