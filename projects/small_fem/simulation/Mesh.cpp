#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "ElementSolution.h"
#include "Mesh.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  // Get Mesh //
  cout << "## Reading Mesh" << endl << flush;

  Mesh   msh(option.getValue("-msh")[1]);
  size_t physical = atoi(option.getValue("-phys")[1].c_str());

  GroupOfElement           domain  = msh.getFromPhysical(physical);
  vector<const MElement*>  element = domain.getAll();
  size_t                  nElement = element.size();

  // Full Mesh //
  cout << "## Full Mesh data" << endl << flush
       << msh.toString()      << endl << flush;

  // Orientations //
  vector<double> orientation(nElement);
  for(size_t i = 0; i < nElement; i++)
    orientation[i] = ReferenceSpaceManager::getOrientation(*element[i]);

  // Mesh //
  size_t nVertex;
  cout << "## Mesh" << endl << flush;

  for(size_t i = 0; i < nElement; i++){
    nVertex = element[i]->getNumPrimaryVertices();

    cout << "  -- " << "Element " << element[i]->getNum()
         << ": "    << "[";

    for(size_t j = 0; j < nVertex - 1; j++)
      cout << element[i]->getVertex(j)->getNum() << ", ";

    cout << element[i]->getVertex(nVertex - 1)->getNum() << "]"
         << ": #" << orientation[i]
         << endl;
  }

  // Printing //
  stringstream name;
  name << "analyze" << physical;

  cout << "## Writing Results" << endl << flush;

  ElementSolution sol;
  sol.addValues(0, 0, domain, orientation);
  sol.write(name.str());

  // Done //
  cout << "## Done" << endl << flush;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-phys");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
