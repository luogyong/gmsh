#include <complex>
#include <iostream>
#include <cmath>

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

#include "FormulationSommerfeld.h"
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

using namespace std;

Complex f(fullVector<double>& xyz){
  return Complex(+1, -1) * (sin(10 * xyz(0)) +
                            sin(10 * xyz(1)) +
                            sin(10 * xyz(2)));
}
/*
fullVector<Complex> f(fullVector<double>& xyz){
  fullVector<Complex> res(3);

  res(0) =  Complex(+1, -1) * sin(10 * xyz(0));
  res(1) =  Complex(+1, -1) * sin(10 * xyz(1));
  res(2) =  Complex(+1, -1) * sin(10 * xyz(2));

  return res;
}
*/
void compute(const Options& option){
  // Get FEM Orders //
  const size_t nOrder = option.getValue("-o").size() - 1;
  vector<int>   order(nOrder);

  for(size_t i = 0; i < nOrder; i++)
    order[i] = atoi(option.getValue("-o")[i + 1].c_str());

  // Get FEM Meshes //
  const size_t  nMesh = option.getValue("-msh").size() - 1;
  vector<string> mesh(nMesh);

  for(size_t i = 0; i < nMesh; i++)
    mesh[i] = option.getValue("-msh")[i + 1];

  // Iterate on Meshes //
  for(size_t i = 0; i < nMesh; i++){
    cout << " ** Mesh: " << mesh[i] << endl << flush;
    Mesh           msh(mesh[i]);
    GroupOfElement domain = msh.getFromPhysical(7);

    // Iterate on Orders
    for(size_t j = 0; j < nOrder; j++){
      cout << "  -- Order " << order[j] << ": " << flush;

      // Projection
      FunctionSpaceScalar fSpace(domain, order[j]);
      FormulationProjection<Complex> projection(domain, fSpace, f);
      System<Complex> sysProj;

      sysProj.addFormulation(projection);

      // Assemble and Solve //
      sysProj.assemble();
      sysProj.solve();

      // Post-processing //
      FEMSolution<Complex> feSol;
      stringstream stream;
      stream << "projection_Mesh" << domain.getNumber() << "_Order" << order[j];

      sysProj.getSolution(feSol, fSpace, domain);
      feSol.write(stream.str());
    }
  }
}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("-msh,-o");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
