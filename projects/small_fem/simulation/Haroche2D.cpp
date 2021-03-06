#include "SmallFem.h"
#include "MPIOStream.h"
#include "HarocheHelper2D.h"
#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include <iostream>
#include <cstdio>

using namespace std;

typedef FormulationStiffness<Complex> FStif;
typedef FormulationMass<Complex>      FMass;

void dump(string filename, fullVector<Complex>& eig){
  const double Pi = 4 * atan(1);
  FILE*      file = fopen(filename.c_str(), "w");

  fprintf(file, "Eig\tReal(Omega^2)\tImag(Omega^2)\tf[MHz]\tt[ms]\n");
  for(int i = 0; i < eig.size(); i++)
    fprintf(file, "%d:\t%.16e\t%.16e\t%.16e\t%.16e\n",
            i,
            eig(i).real(), eig(i).imag(),
            sqrt(eig(i)).real()  / (2 * Pi) * 1e-9,
            1.0 / (sqrt(eig(i)).imag() / (2 * Pi)) * 1e3);

  fclose(file);
}

double getPeakMemory(void){
  // Stream //
  ifstream stream("/proc/self/status", ifstream::in);
  char        tmp[1048576];
  double   vmPeak;

  // Is open ? //
  if(!stream.is_open())
    throw Exception("Haroche2D: cannot open /proc/self/status for VmPeak");

  // Look for "VmPeak:" //
  stream >> tmp;
  while(strncmp(tmp, "VmPeak:", 1048576) != 0){
    stream.getline(tmp, 1048576);
    stream >> tmp;
  }

  // Read VmPeak //
  stream >> vmPeak;

  // Close & Return //
  stream.close();
  return vmPeak / 1024;
}

fullVector<Complex> fZero(fullVector<double>& xyz){
  fullVector<Complex> f(3);

  f(0) = Complex(0, 0);
  f(1) = Complex(0, 0);
  f(2) = Complex(0, 0);

  return f;
}

void fImpedance(fullVector<double>& xyz, fullMatrix<Complex>& a){
  const double       Pi = 4 * atan(1);      // Pi
  const double       Z0 = 119.9169832 * Pi; // Air impedance
  const double    sigma = 5.96e7;           // Copper conductivity (S/m)
  const double        k = PML::getK();      // Harcohe wavenumber
  const Complex       I = Complex(0, 1);    // Imaginary number
  const Complex     imp = -I/k * std::sqrt(Complex(1, 0) + I/k * Z0 * sigma);

  a(0, 0) = imp; a(0, 1) = 0;   a(0, 2) = 0;
  a(1, 0) = 0;   a(1, 1) = imp; a(1, 2) = 0;
  a(2, 0) = 0;   a(2, 1) = 0;   a(2, 2) = imp;
}

void compute(const Options& option){
  // MPI std::cout //
  MPIOStream cout(0, std::cout);
  int        myProc;
  int        nProcs;

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  // Get PML Data //
  cout << "Reading PML... " << endl << flush;
  PML::read(option.getValue("-pml")[1]);

  // Get Domains //
  cout << "Reading domain... " << endl << flush;
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement* Air;
  GroupOfElement* PMLx;
  GroupOfElement* PMLxy;
  GroupOfElement* PMLy;
  GroupOfElement* Mirror;
  GroupOfElement* Frame;
  GroupOfElement* LineOY;
  GroupOfElement* OutPML;

  if(nProcs != 1){
    Air    = new GroupOfElement(msh.getFromPhysical(4000, myProc + 1));
    PMLx   = new GroupOfElement(msh.getFromPhysical(1000, myProc + 1));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(2000, myProc + 1));
    PMLy   = new GroupOfElement(msh.getFromPhysical(3000, myProc + 1));
    Mirror = new GroupOfElement(msh.getFromPhysical( 103, myProc + 1));
    Frame  = new GroupOfElement(msh.getFromPhysical( 105, myProc + 1));
    LineOY = new GroupOfElement(msh.getFromPhysical( 101, myProc + 1));
    OutPML = new GroupOfElement(msh.getFromPhysical( 104, myProc + 1));
  }

  else{
    Air    = new GroupOfElement(msh.getFromPhysical(4000));
    PMLx   = new GroupOfElement(msh.getFromPhysical(1000));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(2000));
    PMLy   = new GroupOfElement(msh.getFromPhysical(3000));
    Mirror = new GroupOfElement(msh.getFromPhysical( 103));
    Frame  = new GroupOfElement(msh.getFromPhysical( 105));
    LineOY = new GroupOfElement(msh.getFromPhysical( 101));
    OutPML = new GroupOfElement(msh.getFromPhysical( 104));
  }

  // Full Domain
  vector<const GroupOfElement*> All_domains(8);
  All_domains[0] = Air;
  All_domains[1] = PMLx;
  All_domains[2] = PMLxy;
  All_domains[3] = PMLy;
  All_domains[4] = Mirror;
  All_domains[5] = Frame;
  All_domains[6] = LineOY;
  All_domains[7] = OutPML;

  // FunctionSpace //
  cout << "FunctionSpace... " << endl << flush;
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  // Formulations //
  // Volume
  cout << "Formulations... " << endl << flush;
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

  // Impedance boundary condition
  Formulation<Complex>* impedance = new FMass(*Frame, fs, fs, fImpedance);

  // System //
  cout << "System... " << endl << flush;
  SystemEigen sys;

  sys.addFormulation(*stifAir);
  sys.addFormulation(*stifXY);
  sys.addFormulation(*stifX);
  sys.addFormulation(*stifY);

  sys.addFormulationB(*massAir);
  sys.addFormulationB(*massXY);
  sys.addFormulationB(*massX);
  sys.addFormulationB(*massY);

  sys.addFormulationB(*impedance);

  // Dirichlet //
  cout << "Dirichlet... " << endl << flush;
  SystemHelper<Complex>::dirichlet(sys, fs, *Mirror, fZero);
  //SystemHelper<Complex>::dirichlet(sys, fs, *Frame , fZero);
  SystemHelper<Complex>::dirichlet(sys, fs, *LineOY, fZero);
  SystemHelper<Complex>::dirichlet(sys, fs, *OutPML, fZero);

  // Assemble //
  cout << "True assembling... " << endl << flush;
  sys.assemble();

  // Free formulations //
  delete stifAir;
  delete stifXY;
  delete stifX;
  delete stifY;

  delete massAir;
  delete massXY;
  delete massX;
  delete massY;

  delete impedance;

  // Solve //
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
  //sys.setProblem("pos_gen_non_hermitian");
  sys.setProblem("gen_non_hermitian");

  cout << "Solving: " << sys.getSize() << "..." << endl << flush;
  sys.solve();

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
  delete Frame;
  delete LineOY;
  delete OutPML;

  // Give peak virtual memory //
  cout << "Process " << myProc << " peak VM: " << getPeakMemory() << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-pml,-o,-n,-shift,-tol,-maxit,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
