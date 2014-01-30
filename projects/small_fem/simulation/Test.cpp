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
#include "FormulationSteadyWaveScalar.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;


int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("-msh,-o,-k");
  SmallFem::Initialize(argc, argv);


  SmallFem::Finalize();
}
