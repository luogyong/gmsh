#include "SmallFem.h"
#include "Mesh.h"
#include "System.h"
#include "SystemPETSc.h"
#include "SystemHelper.h"
#include "FormulationPoisson.h"

#include <iostream>
#include <sstream>

using namespace std;

double fDirichlet0(fullVector<double>& xyz){
  return 1;
}

double fDirichlet1(fullVector<double>& xyz){
  return 0;
}

double fSource(fullVector<double>& xyz){
  return 0;
}

void fMaterial(fullVector<double>& xyz, fullMatrix<double>& tensor){
  tensor.scale(0);

  if(xyz(0) > 0){
    tensor(0, 0) = 1;
    tensor(1, 1) = 1;
    tensor(2, 2) = 1;
  }
  if(xyz(0) <= 0){
    tensor(0, 0) = 2;
    tensor(1, 1) = 2;
    tensor(2, 2) = 2;
  }
}

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement    volume = msh.getFromPhysical(7);
  GroupOfElement boundary0 = msh.getFromPhysical(6);
  GroupOfElement boundary1 = msh.getFromPhysical(5);

  // Full Domain //
  vector<const GroupOfElement*> domain(3);
  domain[0] = &volume;
  domain[1] = &boundary0;
  domain[2] = &boundary1;

  // Get Order //
  size_t order = atoi(option.getValue("-o")[1].c_str());

  // Function Space //
  FunctionSpaceScalar fs(domain, order);

  // Compute //
  FormulationPoisson poisson(volume, fs, fSource, fMaterial);

  System<double> sysPoisson;
  sysPoisson.addFormulation(poisson);

  SystemHelper<double>::dirichlet(sysPoisson, fs, boundary0, fDirichlet0);
  SystemHelper<double>::dirichlet(sysPoisson, fs, boundary1, fDirichlet1);

  cout << "Assembling..." << endl;
  sysPoisson.assemble();

  cout << "Poisson -- Order " << order
       << ": " << sysPoisson.getSize()
       << endl << flush;

  cout << "Solving..." << endl;
  sysPoisson.solve();

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<double> feSol;
    sysPoisson.getSolution(feSol, fs, volume);
    feSol.write("poisson");
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
