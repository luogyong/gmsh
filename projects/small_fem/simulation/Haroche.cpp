#include <iostream>
#include <cstdio>

#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationStiffness.h"
#include "FormulationMass.h"

#include "HarocheHelper.h"
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
  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement Air    = msh.getFromPhysical(138);

  GroupOfElement PMLx   = msh.getFromPhysical(139);
  GroupOfElement PMLxy  = msh.getFromPhysical(140);
  GroupOfElement PMLy   = msh.getFromPhysical(141);
  GroupOfElement PMLz   = msh.getFromPhysical(142);
  GroupOfElement PMLxyz = msh.getFromPhysical(143);
  GroupOfElement PMLxz  = msh.getFromPhysical(144);
  GroupOfElement PMLyz  = msh.getFromPhysical(145);

  GroupOfElement Mirror = msh.getFromPhysical(148);
  GroupOfElement SurfYZ = msh.getFromPhysical(147);
  GroupOfElement SurfXZ = msh.getFromPhysical(146);
  GroupOfElement SurfXY = msh.getFromPhysical(149);

  // Full Domain
  GroupOfElement All_domains(msh);
  All_domains.add(Air);

  All_domains.add(PMLx);
  All_domains.add(PMLxy);
  All_domains.add(PMLy);
  All_domains.add(PMLz);
  All_domains.add(PMLxyz);
  All_domains.add(PMLxz);
  All_domains.add(PMLyz);

  All_domains.add(SurfYZ);
  All_domains.add(SurfXZ);
  All_domains.add(SurfXY);

  // Full Volume
  GroupOfElement All_volumes(msh);
  All_volumes.add(Air);

  All_volumes.add(PMLx);
  All_volumes.add(PMLxy);
  All_volumes.add(PMLy);
  All_volumes.add(PMLz);
  All_volumes.add(PMLxyz);
  All_volumes.add(PMLxz);
  All_volumes.add(PMLyz);

  // Full Surface
  GroupOfElement All_surfaces(msh);
  All_surfaces.add(SurfYZ);
  All_surfaces.add(SurfXZ);
  All_surfaces.add(SurfXY);

  // FunctionSpace //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  cout << "Number of Dofs: " << fs.getAllDofs().size() << endl;

  // Formulation //
  cout << "Assembling" << endl << flush;

  Formulation<Complex>* stifAir = new FormulationStiffness<Complex>
                                       (Air,    fs, fs, Material::Air::Nu);
  Formulation<Complex>* stifXYZ = new FormulationStiffness<Complex>
                                       (PMLxyz, fs, fs, Material::XYZ::Nu);
  Formulation<Complex>* stifXY  = new FormulationStiffness<Complex>
                                       (PMLxy,  fs, fs, Material::XY::Nu);
  Formulation<Complex>* stifYZ  = new FormulationStiffness<Complex>
                                       (PMLyz,  fs, fs, Material::YZ::Nu);
  Formulation<Complex>* stifXZ  = new FormulationStiffness<Complex>
                                       (PMLxz,  fs, fs, Material::XZ::Nu);
  Formulation<Complex>* stifX   = new FormulationStiffness<Complex>
                                       (PMLx,   fs, fs, Material::X::Nu);
  Formulation<Complex>* stifY   = new FormulationStiffness<Complex>
                                       (PMLy,   fs, fs, Material::Y::Nu);
  Formulation<Complex>* stifZ   = new FormulationStiffness<Complex>
                                       (PMLz,   fs, fs, Material::Z::Nu);

  Formulation<Complex>* massAir = new FormulationMass<Complex>
                                       (Air,    fs, fs, Material::Air::Epsilon);
  Formulation<Complex>* massXYZ = new FormulationMass<Complex>
                                       (PMLxyz, fs, fs, Material::XYZ::Epsilon);
  Formulation<Complex>* massXY  = new FormulationMass<Complex>
                                       (PMLxy,  fs, fs, Material::XY::Epsilon);
  Formulation<Complex>* massYZ  = new FormulationMass<Complex>
                                       (PMLyz,  fs, fs, Material::YZ::Epsilon);
  Formulation<Complex>* massXZ  = new FormulationMass<Complex>
                                       (PMLxz,  fs, fs, Material::XZ::Epsilon);
  Formulation<Complex>* massX   = new FormulationMass<Complex>
                                       (PMLx,   fs, fs, Material::X::Epsilon);
  Formulation<Complex>* massY   = new FormulationMass<Complex>
                                       (PMLy,   fs, fs, Material::Y::Epsilon);
  Formulation<Complex>* massZ   = new FormulationMass<Complex>
                                       (PMLz,   fs, fs, Material::Z::Epsilon);
  // System //
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
  // Mirror
  SystemHelper<Complex>::dirichlet(sys, fs, Mirror, fZero);

  // Symmetry
  try{
    int sym = atoi(option.getValue("-sym")[1].c_str());

    if(sym == 0)
      SystemHelper<Complex>::dirichlet(sys, fs, SurfYZ, fZero);

    else
      SystemHelper<Complex>::dirichlet(sys, fs, SurfXZ, fZero);
  }

  catch(...){
    // If no symmetry given, use YZ
    cout << "No symmetry given: defaulting to YZ" << endl;
    SystemHelper<Complex>::dirichlet(sys, fs, SurfYZ, fZero);
  }

  // Assemble //
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
  cout << "Solving: " << sys.getSize() << endl << flush;

  // Set number of eigenvalue (if any, else default)
  try{
    const size_t nWave = atoi(option.getValue("-n")[1].c_str());
    sys.setNumberOfEigenValues(nWave);
  }

  catch(...){
  }

  // Set shift (if any, else default)
  try{
    double shift = atof(option.getValue("-shift")[1].c_str());

    sys.setWhichEigenpairs("target_magnitude");
    sys.setTarget(Complex(shift, 0));
  }

  catch(...){
  }

  // Do what you have to do !
  sys.solve();

  // Post-Pro //
  // Display
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

  // Draw
  try{
    option.getValue("-nopos");
  }

  catch(...){
    FEMSolution<Complex> feSol;
    sys.getSolution(feSol, fs, All_surfaces);
    feSol.write("harocheModes");
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-shift,-sym,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
