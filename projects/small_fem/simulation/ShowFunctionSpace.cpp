#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "FEMSolution.h"
#include "fullMatrix.h"
#include "DofManager.h"
#include "Mesh.h"

#include "Exception.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement domain =
    msh.getFromPhysical(atoi(option.getValue("-phys")[1].c_str()));

  // Get FunctionSpace //
  const size_t   order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpace* fSpace;

  if(option.getValue("-type")[1].compare("lagrange") == 0)
    fSpace = new FunctionSpaceScalar(domain, order, "lagrange");

  else if(option.getValue("-type")[1].compare("scalar") == 0)
    fSpace = new FunctionSpaceScalar(domain, order);

  else if(option.getValue("-type")[1].compare("vector") == 0)
    fSpace = new FunctionSpaceVector(domain, order);

  else
    throw Exception("Unknown FunctionSpace type: %s",
                    option.getValue("-type")[1].c_str());

  // Get Dofs //
  set<Dof> dof;
  fSpace->getKeys(domain, dof);

  // Enumerate Dofs //
  DofManager<double> dofM;
  dofM.addToDofManager(dof);
  dofM.generateGlobalIdSpace();

  // FunctionSpace is solution of FEM problem                         //
  // For every global basis, every dofs are equal to zero excepte one //

  // Dof - Coef Map
  set<Dof>::iterator dEnd = dof.end();
  set<Dof>::iterator  dIt = dof.begin();
  map<Dof, double>  coef;

  for(; dIt != dEnd; dIt++)
    coef.insert(pair<Dof, double>(*dIt, 0));

  // FEM Solutions
  FEMSolution<double>        sol;
  map<Dof, double>::iterator cEnd = coef.end();
  map<Dof, double>::iterator  cIt = coef.begin();

  for(size_t i = 0; cIt != cEnd; cIt++, i++){
    cIt->second = 1;
    sol.addCoefficients(i, i, domain, *fSpace, coef);
    cIt->second = 0;
  }

  // Write //
  sol.write("function_space");

  // Clean //
  delete fSpace;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-type,-phys");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
