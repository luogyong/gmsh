#include <iostream>
#include <cstdio>

#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include "HarocheHelper2D.h"
#include "MPIOStream.h"
#include "SmallFem.h"

using namespace std;

void dump(string filename, fullVector<Complex>& eig){
  FILE* file = fopen(filename.c_str(), "w");
  double pi  = 4 * atan(1);

  fprintf(file, "Eig\tReal(Omega^2)\tImag(Omega^2)\tf[MHz]\tt[ms]\n");
  for(int i = 0; i < eig.size(); i++)
    fprintf(file, "%d:\t%.16e\t%.16e\t%.16e\t%.16e\n",
            i,
            eig(i).real(), eig(i).imag(),
            sqrt(eig(i)).real()  / (2 * pi) * 1e-9,
            1.0 / (sqrt(eig(i)).imag() / (2 * pi)) * 1e3);

  fclose(file);
}

fullVector<Complex> fZero(fullVector<double>& xyz){
  fullVector<Complex> f(3);

  f(0) = Complex(0, 0);
  f(1) = Complex(0, 0);
  f(2) = Complex(0, 0);

  return f;
}

void compute(const Options& option){
  // MPI std::cout //
  MPIOStream cout(0, std::cout);
  int        myProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);

  // Get Domains //
  cout << "Reading domain... " << flush;
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement Air    = msh.getFromPhysical(4000, myProc + 1);

  GroupOfElement PMLx   = msh.getFromPhysical(1000, myProc + 1);
  GroupOfElement PMLxy  = msh.getFromPhysical(2000, myProc + 1);
  GroupOfElement PMLy   = msh.getFromPhysical(3000, myProc + 1);

  GroupOfElement Mirror = msh.getFromPhysical(103, myProc + 1);
  GroupOfElement LineOY = msh.getFromPhysical(101, myProc + 1);
  GroupOfElement LineOX = msh.getFromPhysical(102, myProc + 1);

  /*
  GroupOfElement Air    = msh.getFromPhysical(4000);

  GroupOfElement PMLx   = msh.getFromPhysical(1000);
  GroupOfElement PMLxy  = msh.getFromPhysical(2000);
  GroupOfElement PMLy   = msh.getFromPhysical(3000);

  GroupOfElement Mirror = msh.getFromPhysical(103);
  GroupOfElement LineOY = msh.getFromPhysical(101);
  GroupOfElement LineOX = msh.getFromPhysical(102);
  */
  // Full Domain
  vector<const GroupOfElement*> All_domains(7);
  All_domains[0] = &Air;

  All_domains[1] = &PMLx;
  All_domains[2] = &PMLxy;
  All_domains[3] = &PMLy;

  All_domains[4] = &Mirror;
  All_domains[5] = &LineOY;
  All_domains[6] = &LineOX;

  // Full Surface
  cout << "Done!" << endl << flush;

  // FunctionSpace //
  cout << "FunctionSpace... " << flush;
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  cout << "Done!" << endl << flush;

  // Formulation //
  cout << "Formulations... " << flush;

  Formulation<Complex>* stifAir = new FormulationStiffness<Complex>
                                       (Air,    fs, fs, Material::Air::Nu);
  Formulation<Complex>* stifXY  = new FormulationStiffness<Complex>
                                       (PMLxy,  fs, fs, Material::XY::Nu);
  Formulation<Complex>* stifX   = new FormulationStiffness<Complex>
                                       (PMLx,   fs, fs, Material::X::Nu);
  Formulation<Complex>* stifY   = new FormulationStiffness<Complex>
                                       (PMLy,   fs, fs, Material::Y::Nu);

  Formulation<Complex>* massAir = new FormulationMass<Complex>
                                       (Air,    fs, fs, Material::Air::Epsilon);
  Formulation<Complex>* massXY  = new FormulationMass<Complex>
                                       (PMLxy,  fs, fs, Material::XY::Epsilon);
  Formulation<Complex>* massX   = new FormulationMass<Complex>
                                       (PMLx,   fs, fs, Material::X::Epsilon);
  Formulation<Complex>* massY   = new FormulationMass<Complex>
                                       (PMLy,   fs, fs, Material::Y::Epsilon);
  cout << "Done!" << endl << flush;

  // System //
  cout << "System... " << flush;
  SystemEigen sys;

  sys.addFormulation(*stifAir);
  sys.addFormulation(*stifXY);
  sys.addFormulation(*stifX);
  sys.addFormulation(*stifY);

  sys.addFormulationB(*massAir);
  sys.addFormulationB(*massXY);
  sys.addFormulationB(*massX);
  sys.addFormulationB(*massY);
  cout << "Done!" << endl << flush;

  // Dirichlet //
  // Mirror
  cout << "Dirichlet... " << flush;
  SystemHelper<Complex>::dirichlet(sys, fs, Mirror, fZero);

  // Symmetry
  try{
    int sym = atoi(option.getValue("-sym")[1].c_str());

    if(sym == 0)
      SystemHelper<Complex>::dirichlet(sys, fs, LineOY, fZero);

    else
      SystemHelper<Complex>::dirichlet(sys, fs, LineOX, fZero);
  }

  catch(...){
    // If no symmetry given, use YZ
    cout << "No symmetry given: defaulting to YZ" << endl;
    SystemHelper<Complex>::dirichlet(sys, fs, LineOY, fZero);
  }
  cout << "Done!" << endl << flush;

  // Assemble //
  cout << "True assembling... " << endl << flush;
  sys.assemble();
  cout << "Done!" << endl << flush;

  // Free formulations //
  cout << "Clearing..." << endl << flush;
  delete stifAir;
  delete stifXY;
  delete stifX;
  delete stifY;

  delete massAir;
  delete massXY;
  delete massX;
  delete massY;
  cout << "Done!" << endl << flush;

  // Solve //
  cout << "Solving: " << sys.getSize() << endl << flush;

  // Set number of eigenvalue (if any, else default)
  try{
    sys.setNumberOfEigenValues(atoi(option.getValue("-n")[1].c_str()));
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

  // Set maximun iteration number (if any, else default)
  try{
    sys.setMaxIteration(atoi(option.getValue("-maxit")[1].c_str()));
  }

  catch(...){
  }

  // Do what you have to do !
  sys.solve();
  cout << "Done!" << endl << flush;

  // Post-Pro //
  if(myProc == 0){
    // Display
    cout << "Post Pro" << endl << flush;
    const size_t nEigenValue = sys.getNComputedSolution();
    fullVector<Complex> eigenValue;

    sys.getEigenValues(eigenValue);

    cout << "Number of found Eigenvalues: " << nEigenValue
         << endl
         << endl
         << "Number\tEigen Value" << endl;

    for(size_t i = 0; i < nEigenValue; i++)
      cout << "#" << i + 1  << "\t"
           << std::scientific
           << eigenValue(i) << endl;

    // Dump on disk
    dump("harocheValues.txt", eigenValue);
  }

  // Draw
  try{
    option.getValue("-nopos");
  }

  catch(...){
    FEMSolution<Complex> feSol;
    stringstream         name;
    sys.getSolution(feSol, fs, All_domains);

    name << "harocheModes" << "_proc" << myProc;
    feSol.write(name.str());
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-shift,-sym,-tol,-maxit,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
