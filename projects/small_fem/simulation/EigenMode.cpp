#include "SmallFem.h"
#include "MPIOStream.h"
#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"
#include "FormulationSteadyWave.h"

#include "TextSolution.h"

#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

static const size_t scal = 0;
static const size_t vect = 1;
static const double Pi   = atan(1.0) * 4;
static const double c0   = 299792458;
static const double mu0  = 4 * Pi * 1e-7;
static const double eps0 = 1.0 / (mu0 * c0 * c0);
static const double nu0  = 1.0 / mu0;

void nu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = nu0 * 1e3; // Scaling for millimeters
  tensor(1, 1) = nu0 * 1e3; // Scaling for millimeters
  tensor(2, 2) = nu0 * 1e3; // Scaling for millimeters
}

void eps(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = eps0 / 1e3; // Scaling for millimeters
  tensor(1, 1) = eps0 / 1e3; // Scaling for millimeters
  tensor(2, 2) = eps0 / 1e3; // Scaling for millimeters
}

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
  int        nProcs;

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  // Get Domain //
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement* volume;
  GroupOfElement* border;

  if(nProcs == 1){
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
  FormulationStiffness<Complex> stiff(*volume, *fs, *fs, nu);
  FormulationMass<Complex>       mass(*volume, *fs, *fs, eps);
  //FormulationSteadyWave<Complex> wave(*volume, *fs, 4.94967);

  SystemEigen sys;
  sys.addFormulation(stiff);
  sys.addFormulationB(mass);
  //sys.addFormulation(wave);

  // Dirichlet //
  if(type == scal)
    SystemHelper<Complex>::dirichlet(sys, *fs, *border, fScal);
  else
    SystemHelper<Complex>::dirichlet(sys, *fs, *border, fVect);

  // Assemble and Solve //
  cout << "Eigenvalues problem" << endl << flush;

  sys.assemble();
  cout << "Assembled: " << sys.getSize() << endl << flush;

  // Just dump matrices ? //
  try{
    string name = option.getValue("-dump")[1];
    string dummy;

    dummy.clear();

    sys.writeMatrix(name, dummy);
  }

  catch(...){
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

    // Set tolerance (if any, else default)
    try{
      sys.setTolerance(atof(option.getValue("-tol")[1].c_str()));
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
      // Fields
      FEMSolution<Complex> feSol;
      sys.getSolution(feSol, *fs, *volume);

      feSol.setSaveMesh(false);
      feSol.setBinaryFormat(true);
      if(nProcs != 1)
        feSol.setParition(myProc + 1);

      feSol.write("eigenModes");

      // Eigenvalues
      const size_t        nEigenValue = sys.getNComputedSolution();
      TextSolution        eigenSol;
      fullVector<Complex> eigenValue;

      sys.getEigenValues(eigenValue);

      for(size_t i = 0; i < nEigenValue; i++){
        stringstream stream;
        stream << "Eigen value " << i + 1 << ": " << eigenValue(i) << " ";

        eigenSol.addValues(i * 2 + 0, stream.str().append("(real part)"));
        eigenSol.addValues(i * 2 + 1, stream.str().append("(imaginary part)"));
      }

      eigenSol.write("eigenModes");
    }
  }

  // Clean //
  delete fs;
  delete volume;
  delete border;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-shift,-nopos,-type,-tol,-dump");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
