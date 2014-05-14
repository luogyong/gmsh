#include <iostream>
#include <complex>

#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include "SmallFem.h"

using namespace std;

static const size_t scal = 0;
static const size_t vect = 1;

void nu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(xyz(0) + 10, 0);
  tensor(1, 1) = Complex(xyz(1) + 10, 0);
  tensor(2, 2) = Complex(xyz(2) + 10, 0);
}

void eps(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(xyz(0) + 10, 0);
  tensor(1, 1) = Complex(xyz(1) + 10, 0);
  tensor(2, 2) = Complex(xyz(2) + 10, 0);
}

/*
void nu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(1, 0);
  tensor(1, 1) = Complex(2, 0);
  tensor(2, 2) = Complex(3, 0);
}

void eps(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(1, 0);
  tensor(1, 1) = Complex(2, 0);
  tensor(2, 2) = Complex(3, 0);
}
*/
fullVector<Complex> fVect(fullVector<double>& xyz){
  fullVector<Complex> f(3);

  f(0) = Complex(0, 0);
  f(1) = Complex(0, 0);
  f(2) = Complex(0, 0);

  return f;
}

Complex fScal(fullVector<double>& xyz){
  return Complex(0, 0);
}

void compute(const Options& option){
  // Get Domain //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(border);

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  const size_t nWave = atoi(option.getValue("-n")[1].c_str());

  // Get Type //
  size_t type;

  if(option.getValue("-type")[1].compare("scalar") == 0){
    type = scal;
    cout << "Scalar ";
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    type = vect;
    cout << "Vectorial ";
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Function Space //
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Formulations & System //
  FormulationStiffness<Complex> stiff(volume, *fs, *fs, nu);
  FormulationMass<Complex>       mass(volume, *fs, *fs, eps);

  SystemEigen sys;
  sys.addFormulation(stiff);
  sys.addFormulationB(mass);

  // Dirichlet //
  if(type == scal)
    SystemHelper<Complex>::dirichlet(sys, *fs, border, fScal);
  else
    SystemHelper<Complex>::dirichlet(sys, *fs, border, fVect);

  // Assemble and Solve //
  cout << "Eigenvalues problem: " << sys.getSize() << endl
       << "Assembling..."         << endl          << flush;
  sys.assemble();

  cout << "Solving..." << endl << flush;
  sys.setNumberOfEigenValues(nWave);
  sys.solve();

  // Display //
  fullVector<Complex> eigenValue;
  const size_t nEigenValue = sys.getNComputedSolution();
  sys.getEigenValues(eigenValue);

  cout << "Number of found Eigenvalues: " << nEigenValue
       << endl
       << endl
       << "Number\tEigen Value" << endl;

  for(size_t i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1  << "\t"
         << eigenValue(i) << endl;

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<Complex> feSol;
    sys.getSolution(feSol, *fs, volume);
    feSol.write("eigen_mode");
  }

  // Clean //
  delete fs;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-nopos,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
