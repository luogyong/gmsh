#include "SmallFem.h"
#include "MPIOStream.h"
#include "HarocheHelper.h"
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
  GroupOfElement* PMLz;
  GroupOfElement* PMLxyz;
  GroupOfElement* PMLxz;
  GroupOfElement* PMLyz;
  GroupOfElement* Mirror;
  GroupOfElement* SurfYZ;
  GroupOfElement* SurfXZ;
  GroupOfElement* SurfXY;
  GroupOfElement* OutPML;

  if(nProcs != 1){
    Air    = new GroupOfElement(msh.getFromPhysical(138, myProc + 1));
    PMLx   = new GroupOfElement(msh.getFromPhysical(139, myProc + 1));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(140, myProc + 1));
    PMLy   = new GroupOfElement(msh.getFromPhysical(141, myProc + 1));
    PMLz   = new GroupOfElement(msh.getFromPhysical(142, myProc + 1));
    PMLxyz = new GroupOfElement(msh.getFromPhysical(143, myProc + 1));
    PMLxz  = new GroupOfElement(msh.getFromPhysical(144, myProc + 1));
    PMLyz  = new GroupOfElement(msh.getFromPhysical(145, myProc + 1));
    Mirror = new GroupOfElement(msh.getFromPhysical(148, myProc + 1));
    SurfYZ = new GroupOfElement(msh.getFromPhysical(147, myProc + 1));
    SurfXZ = new GroupOfElement(msh.getFromPhysical(146, myProc + 1));
    SurfXY = new GroupOfElement(msh.getFromPhysical(149, myProc + 1));
    OutPML = new GroupOfElement(msh.getFromPhysical(150, myProc + 1));
  }

  else{
    Air    = new GroupOfElement(msh.getFromPhysical(138));
    PMLx   = new GroupOfElement(msh.getFromPhysical(139));
    PMLxy  = new GroupOfElement(msh.getFromPhysical(140));
    PMLy   = new GroupOfElement(msh.getFromPhysical(141));
    PMLz   = new GroupOfElement(msh.getFromPhysical(142));
    PMLxyz = new GroupOfElement(msh.getFromPhysical(143));
    PMLxz  = new GroupOfElement(msh.getFromPhysical(144));
    PMLyz  = new GroupOfElement(msh.getFromPhysical(145));
    Mirror = new GroupOfElement(msh.getFromPhysical(148));
    SurfYZ = new GroupOfElement(msh.getFromPhysical(147));
    SurfXZ = new GroupOfElement(msh.getFromPhysical(146));
    SurfXY = new GroupOfElement(msh.getFromPhysical(149));
    OutPML = new GroupOfElement(msh.getFromPhysical(150));
  }

  // Full Domain
  vector<const GroupOfElement*> All_domains(13);
  All_domains[0]  = Air;
  All_domains[1]  = PMLx;
  All_domains[2]  = PMLxy;
  All_domains[3]  = PMLy;
  All_domains[4]  = PMLz;
  All_domains[5]  = PMLxyz;
  All_domains[6]  = PMLxz;
  All_domains[7]  = PMLyz;
  All_domains[8]  = Mirror;
  All_domains[9]  = SurfYZ;
  All_domains[10] = SurfXZ;
  All_domains[11] = SurfXY;
  All_domains[12] = OutPML;

  // Full Surface
  vector<const GroupOfElement*> All_surfaces(4);
  All_surfaces[0] = SurfYZ;
  All_surfaces[1] = SurfXZ;
  All_surfaces[2] = SurfXY;
  All_surfaces[3] = OutPML;

  // FunctionSpace //
  cout << "FunctionSpace... " << endl << flush;
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  // Formulation //
  cout << "Formulations... " << endl << flush;
  Formulation<Complex>* stifAir;
  Formulation<Complex>* stifXYZ;
  Formulation<Complex>* stifXY;
  Formulation<Complex>* stifYZ;
  Formulation<Complex>* stifXZ;
  Formulation<Complex>* stifX;
  Formulation<Complex>* stifY;
  Formulation<Complex>* stifZ;

  Formulation<Complex>* massAir;
  Formulation<Complex>* massXYZ;
  Formulation<Complex>* massXY;
  Formulation<Complex>* massYZ;
  Formulation<Complex>* massXZ;
  Formulation<Complex>* massX;
  Formulation<Complex>* massY;
  Formulation<Complex>* massZ;
  /*
  stifAir = new FStif(*Air,    fs, fs, Material::Air::Nu);
  stifXYZ = new FStif(*PMLxyz, fs, fs, Material::XYZ::Nu);
  stifXY  = new FStif(*PMLxy,  fs, fs,  Material::XY::Nu);
  stifYZ  = new FStif(*PMLyz,  fs, fs,  Material::YZ::Nu);
  stifXZ  = new FStif(*PMLxz,  fs, fs,  Material::XZ::Nu);
  stifX   = new FStif(*PMLx,   fs, fs,   Material::X::Nu);
  stifY   = new FStif(*PMLy,   fs, fs,   Material::Y::Nu);
  stifZ   = new FStif(*PMLz,   fs, fs,   Material::Z::Nu);

  massAir = new FMass(*Air,    fs, fs, Material::Air::Epsilon);
  massXYZ = new FMass(*PMLxyz, fs, fs, Material::XYZ::Epsilon);
  massXY  = new FMass(*PMLxy,  fs, fs,  Material::XY::Epsilon);
  massYZ  = new FMass(*PMLyz,  fs, fs,  Material::YZ::Epsilon);
  massXZ  = new FMass(*PMLxz,  fs, fs,  Material::XZ::Epsilon);
  massX   = new FMass(*PMLx,   fs, fs,   Material::X::Epsilon);
  massY   = new FMass(*PMLy,   fs, fs,   Material::Y::Epsilon);
  massZ   = new FMass(*PMLz,   fs, fs,   Material::Z::Epsilon);
  */

  stifAir = new FStif(*Air,    fs, fs, Material::Air::OverMuEps);
  stifXYZ = new FStif(*PMLxyz, fs, fs, Material::XYZ::OverMuEps);
  stifXY  = new FStif(*PMLxy,  fs, fs,  Material::XY::OverMuEps);
  stifYZ  = new FStif(*PMLyz,  fs, fs,  Material::YZ::OverMuEps);
  stifXZ  = new FStif(*PMLxz,  fs, fs,  Material::XZ::OverMuEps);
  stifX   = new FStif(*PMLx,   fs, fs,   Material::X::OverMuEps);
  stifY   = new FStif(*PMLy,   fs, fs,   Material::Y::OverMuEps);
  stifZ   = new FStif(*PMLz,   fs, fs,   Material::Z::OverMuEps);

  massAir = new FMass(*Air,    fs, fs);
  massXYZ = new FMass(*PMLxyz, fs, fs);
  massXY  = new FMass(*PMLxy,  fs, fs);
  massYZ  = new FMass(*PMLyz,  fs, fs);
  massXZ  = new FMass(*PMLxz,  fs, fs);
  massX   = new FMass(*PMLx,   fs, fs);
  massY   = new FMass(*PMLy,   fs, fs);
  massZ   = new FMass(*PMLz,   fs, fs);

  // System //
  cout << "System... " << endl << flush;
  SystemEigen sys;

  sys.addFormulation(*stifAir);
  sys.addFormulation(*stifXYZ);
  sys.addFormulation(*stifXY);
  sys.addFormulation(*stifYZ);
  sys.addFormulation(*stifXZ);
  sys.addFormulation(*stifX);
  sys.addFormulation(*stifY);
  sys.addFormulation(*stifZ);

  sys.addFormulationB(*massAir);
  sys.addFormulationB(*massXYZ);
  sys.addFormulationB(*massXY);
  sys.addFormulationB(*massYZ);
  sys.addFormulationB(*massXZ);
  sys.addFormulationB(*massX);
  sys.addFormulationB(*massY);
  sys.addFormulationB(*massZ);

  // Dirichlet //
  // Mirror & PML
  cout << "Dirichlet... " << endl << flush;
  SystemHelper<Complex>::dirichlet(sys, fs, *Mirror, fZero);
  SystemHelper<Complex>::dirichlet(sys, fs, *OutPML, fZero);

  // Symmetry
  try{
    int sym = atoi(option.getValue("-sym")[1].c_str());

    if(sym == 0)
      SystemHelper<Complex>::dirichlet(sys, fs, *SurfYZ, fZero);

    else
      SystemHelper<Complex>::dirichlet(sys, fs, *SurfXZ, fZero);
  }

  catch(...){
    // If no symmetry given, use YZ
    cout << "No symmetry given: defaulting to YZ" << endl;
    SystemHelper<Complex>::dirichlet(sys, fs, *SurfYZ, fZero);
  }

  // Assemble //
  cout << "True assembling... " << endl << flush;
  sys.assemble();

  // Free formulations //
  delete stifAir;
  delete stifXYZ;
  delete stifXY;
  delete stifYZ;
  delete stifXZ;
  delete stifX;
  delete stifY;
  delete stifZ;

  delete massAir;
  delete massXYZ;
  delete massXY;
  delete massYZ;
  delete massXZ;
  delete massX;
  delete massY;
  delete massZ;

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
  sys.setProblem("pos_gen_non_hermitian");
  //sys.setProblem("gen_non_hermitian");

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
    sys.getSolution(feSol, fs, All_surfaces);

    feSol.setSaveMesh(false);
    feSol.setBinaryFormat(true);
    if(nProcs != 1)
      feSol.setParition(myProc + 1);

    feSol.write("harocheModes");
  }

  // Give peak virtual memory //
  double  myVmPeak = getPeakMemory();
  double* alVmPeak = new double[nProcs];

  MPI_Allgather(&myVmPeak,1,MPI_DOUBLE, alVmPeak,1,MPI_DOUBLE, MPI_COMM_WORLD);

  cout << "Peak VM:" << endl << flush;
  for(int i = 0; i < nProcs; i++)
    cout << " ** Process " << i << ": " << alVmPeak[i] << " MB"
         << endl << flush;

  // Game over! //
  delete   Air;
  delete   PMLx;
  delete   PMLxy;
  delete   PMLy;
  delete   PMLz;
  delete   PMLxyz;
  delete   PMLxz;
  delete   PMLyz;
  delete   Mirror;
  delete   SurfYZ;
  delete   SurfXZ;
  delete   SurfXY;
  delete   OutPML;
  delete[] alVmPeak;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-pml,-o,-n,-shift,-sym,-tol,-maxit,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
