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

#include "FormulationSilverMuller.h"
#include "FormulationEMDA.h"
#include "FormulationProjection.h"

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

#include <complex>
#include <iostream>
#include <cmath>

using namespace std;

void compute(const Options& option){
  Mesh         mesh(option.getValue("-msh")[1]);
  int          myProc;
  stringstream name;
  ofstream     file;

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  name << "Mesh_proc" << myProc;

  file.open(name.str().c_str(), ifstream::out);
  file << "Mesh for process " << myProc << endl
       << "################ "           << endl
       << mesh.toString()               << endl << flush;
  file.close();
}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("-msh");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
