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

typedef FormulationStiffness<Complex> FStif;
typedef FormulationMass<Complex>      FMass;

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
  int        nProcs;

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  // Get PML Data //
  cout << "Reading PML... " << flush;
  PML::read(option.getValue("-pml")[1]);
  cout << "Done!" << endl << flush;

  // Get Domains //
  cout << "Reading domain... " << flush;
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement* Air;
  GroupOfElement* PMLx;
  GroupOfElement* PMLxy;
  GroupOfElement* PMLy;
  GroupOfElement* Mirror;
  GroupOfElement* LineOY;
  GroupOfElement* OutPML;

  if(nProcs != 1){
    Air    = new GroupOfElement(msh.getFromPhysical(4000, myProc + 1));
    PMLx   = new GroupOfElement(msh.getFromPhysical(1000, myProc + 1));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(2000, myProc + 1));
    PMLy   = new GroupOfElement(msh.getFromPhysical(3000, myProc + 1));
    Mirror = new GroupOfElement(msh.getFromPhysical( 103, myProc + 1));
    LineOY = new GroupOfElement(msh.getFromPhysical( 101, myProc + 1));
    OutPML = new GroupOfElement(msh.getFromPhysical( 104, myProc + 1));
  }

  else{
    Air    = new GroupOfElement(msh.getFromPhysical(4000));
    PMLx   = new GroupOfElement(msh.getFromPhysical(1000));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(2000));
    PMLy   = new GroupOfElement(msh.getFromPhysical(3000));
    Mirror = new GroupOfElement(msh.getFromPhysical( 103));
    LineOY = new GroupOfElement(msh.getFromPhysical( 101));
    OutPML = new GroupOfElement(msh.getFromPhysical( 104));
  }

  // Full Domain
  vector<const GroupOfElement*> All_domains(7);
  All_domains[0] = Air;
  All_domains[1] = PMLx;
  All_domains[2] = PMLxy;
  All_domains[3] = PMLy;
  All_domains[4] = Mirror;
  All_domains[5] = LineOY;
  All_domains[6] = OutPML;

  // True Domain
  vector<const GroupOfElement*> True_domains(3);
  True_domains[0] = Air;
  True_domains[1] = Mirror;
  True_domains[2] = LineOY;

  cout << "Done!" << endl << flush;

  // FunctionSpace //
  cout << "FunctionSpace... " << flush;
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  cout << "Done!" << endl << flush;

  // Formulation //
  cout << "Formulations... " << flush;

  Formulation<Complex>* stifAir;
  Formulation<Complex>* stifXY;
  Formulation<Complex>* stifX;
  Formulation<Complex>* stifY;

  Formulation<Complex>* massAir;
  Formulation<Complex>* massXY;
  Formulation<Complex>* massX;
  Formulation<Complex>* massY;
  /*
  stifAir = new FStif(*Air,   fs, fs, Material::Air::Nu);
  stifXY  = new FStif(*PMLxy, fs, fs,  Material::XY::Nu);
  stifX   = new FStif(*PMLx,  fs, fs,   Material::X::Nu);
  stifY   = new FStif(*PMLy,  fs, fs,   Material::Y::Nu);

  massAir = new FMass(*Air,   fs, fs, Material::Air::Epsilon);
  massXY  = new FMass(*PMLxy, fs, fs,  Material::XY::Epsilon);
  massX   = new FMass(*PMLx,  fs, fs,   Material::X::Epsilon);
  massY   = new FMass(*PMLy,  fs, fs,   Material::Y::Epsilon);
  */

  stifAir = new FStif(*Air,   fs, fs,Material::Air::OverMuEps);
  stifXY  = new FStif(*PMLxy, fs, fs, Material::XY::OverMuEps);
  stifX   = new FStif(*PMLx,  fs, fs,  Material::X::OverMuEps);
  stifY   = new FStif(*PMLy,  fs, fs,  Material::Y::OverMuEps);

  massAir = new FMass(*Air,   fs, fs);
  massXY  = new FMass(*PMLxy, fs, fs);
  massX   = new FMass(*PMLx,  fs, fs);
  massY   = new FMass(*PMLy,  fs, fs);

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
  cout << "Dirichlet... " << flush;
  SystemHelper<Complex>::dirichlet(sys, fs, *Mirror, fZero);
  SystemHelper<Complex>::dirichlet(sys, fs, *LineOY, fZero);
  SystemHelper<Complex>::dirichlet(sys, fs, *OutPML, fZero);
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
  sys.setProblem("pos_gen_non_hermitian");
  //sys.setProblem("gen_non_hermitian");
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
    sys.getSolution(feSol, fs, All_domains);

    //feSol.setSaveMesh(false);
    feSol.setBinaryFormat(true);
    if(nProcs != 1)
      feSol.setParition(myProc + 1);

    feSol.write("harocheModes");
  }

  // Game over! //
  delete Air;
  delete PMLx;
  delete PMLxy;
  delete PMLy;
  delete Mirror;
  delete LineOY;
  delete OutPML;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-pml,-o,-n,-shift,-tol,-maxit,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
