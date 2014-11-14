#include <iostream>
#include <complex>

#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include "MPIOStream.h"
#include "SmallFem.h"

using namespace std;

static const size_t scal = 0;
static const size_t vect = 1;

/*
void nu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(xyz(0) + 10, xyz(0) - 10);
  tensor(1, 1) = Complex(xyz(1) + 10, xyz(1) - 10);
  tensor(2, 2) = Complex(xyz(2) + 10, xyz(2) - 10);
}

void eps(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Complex(xyz(0) + 10, xyz(0) - 10);
  tensor(1, 1) = Complex(xyz(1) + 10, xyz(1) - 10);
  tensor(2, 2) = Complex(xyz(2) + 10, xyz(2) - 10);
}

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
  // MPI //
  MPIOStream cout(0, std::cout);
  int        myProc;
  int         nProc;

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  // Get Domain //
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement* volume;
  GroupOfElement* border;

  if(nProc == 1){
    volume = new GroupOfElement(msh.getFromPhysical(7));
    border = new GroupOfElement(msh.getFromPhysical(5));
  }

  else{
    volume = new GroupOfElement(msh.getFromPhysical(7, myProc + 1));
    border = new GroupOfElement(msh.getFromPhysical(5, myProc + 1));
  }

  // Full Domain //
  vector<const GroupOfElement*> domain(2);
  domain[0] = volume;
  domain[1] = border;

  // Get Order //
  const size_t order = atoi(option.getValue("-o")[1].c_str());

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
  FormulationStiffness<Complex> stiff(*volume, *fs, *fs);//, nu);
  FormulationMass<Complex>       mass(*volume, *fs, *fs);//, eps);

  SystemEigen sys;
  sys.addFormulation(stiff);
  sys.addFormulationB(mass);

  // Dirichlet //
  if(type == scal)
    SystemHelper<Complex>::dirichlet(sys, *fs, *border, fScal);
  else
    SystemHelper<Complex>::dirichlet(sys, *fs, *border, fVect);

  // Assemble and Solve //
  cout << "Eigenvalues problem" << endl << flush;

  sys.assemble();
  cout << "Assembled: " << sys.getSize() << endl << flush;

  // Set number of eigenvalue (if any, else default)
  try{
    const size_t nWave = atoi(option.getValue("-n")[1].c_str());
    sys.setNumberOfEigenValues(nWave);
  }

  catch(...){
  }

  // Set shift (if any, else default)
  try{
    const double shift = atof(option.getValue("-shift")[1].c_str());

    sys.setWhichEigenpairs("target_magnitude");
    sys.setTarget(Complex(shift, 0));
  }

  catch(...){
  }

  // Solve
  sys.solve();
  cout << "Solved" << endl << flush;

  // Display //
  if(myProc == 0){
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

  }
  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<Complex> feSol;
    stringstream         name;
    sys.getSolution(feSol, *fs, *volume);

    name << "eigen_mode_proc" << myProc;
    feSol.write(name.str());
  }

  // Clean //
  delete fs;
  delete volume;
  delete border;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-shift,-nopos,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
