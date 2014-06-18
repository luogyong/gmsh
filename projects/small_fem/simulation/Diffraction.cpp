#include <iostream>

#include "DiffractionHelper.h"
#include "SmallFem.h"

#include "GroupOfElement.h"
#include "FunctionSpaceVector.h"
#include "FormulationSteadyWave.h"

#include "System.h"
#include "FEMSolution.h"

using namespace std;

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);

  GroupOfElement PMLxyz   = msh.getFromPhysical(1000);
  GroupOfElement PMLxz    = msh.getFromPhysical(1001);
  GroupOfElement PMLyz    = msh.getFromPhysical(1002);
  GroupOfElement PMLxy    = msh.getFromPhysical(1003);
  GroupOfElement PMLz     = msh.getFromPhysical(1004);
  GroupOfElement PMLy     = msh.getFromPhysical(1005);
  GroupOfElement PMLx     = msh.getFromPhysical(1006);
  GroupOfElement Scat_In  = msh.getFromPhysical(1008);
  GroupOfElement Scat_Out = msh.getFromPhysical(1007);

  // Full Domain
  GroupOfElement All_domains(msh);
  All_domains.add(Scat_In);
  All_domains.add(Scat_Out);
  All_domains.add(PMLxyz);
  All_domains.add(PMLxz);
  All_domains.add(PMLyz);
  All_domains.add(PMLxy);
  All_domains.add(PMLz);
  All_domains.add(PMLy);
  All_domains.add(PMLx);

  // FunctionSpace //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  FunctionSpaceVector fs(All_domains, order);

  // Formulation //
  cout << "Assembling" << endl << flush;
  const double k = (Constant::omega0 / Constant::cel);

  FormulationSteadyWave<Complex> in(Scat_In, fs, k,
                                    Material::In::Nu,
                                    Material::In::Epsilon,
                                    Signal::In::source);

  FormulationSteadyWave<Complex> out(Scat_Out, fs, k,
                                     Material::Out::Nu,
                                     Material::Out::Epsilon,
                                     Signal::Out::source);

  FormulationSteadyWave<Complex> xyz(PMLxyz, fs, k,
                                     Material::XYZ::Nu,
                                     Material::XYZ::Epsilon,
                                     Signal::PML::source);

  FormulationSteadyWave<Complex> xz(PMLxz, fs, k,
                                    Material::XZ::Nu,
                                    Material::XZ::Epsilon,
                                    Signal::PML::source);

  FormulationSteadyWave<Complex> yz(PMLyz, fs, k,
                                    Material::YZ::Nu,
                                    Material::YZ::Epsilon,
                                    Signal::PML::source);

  FormulationSteadyWave<Complex> xy(PMLxy, fs, k,
                                    Material::XY::Nu,
                                    Material::XY::Epsilon,
                                    Signal::PML::source);

  FormulationSteadyWave<Complex> x(PMLx, fs, k,
                                   Material::X::Nu,
                                   Material::X::Epsilon,
                                   Signal::PML::source);

  FormulationSteadyWave<Complex> y(PMLy, fs, k,
                                   Material::Y::Nu,
                                   Material::Y::Epsilon,
                                   Signal::PML::source);

  FormulationSteadyWave<Complex> z(PMLz, fs, k,
                                   Material::Z::Nu,
                                   Material::Z::Epsilon,
                                   Signal::PML::source);
  // System //
  System<Complex> sys;

  sys.addFormulation(in);
  sys.addFormulation(out);
  sys.addFormulation(xyz);
  sys.addFormulation(xz);
  sys.addFormulation(yz);
  sys.addFormulation(xy);
  sys.addFormulation(x);
  sys.addFormulation(y);
  sys.addFormulation(z);

  sys.assemble();

  cout << "Solving: " << sys.getSize() << endl << flush;
  sys.solve();

  // Post-Pro //
  cout << "Drawing" << endl << flush;
  FEMSolution<Complex> feSol;
  sys.getSolution(feSol, fs, All_domains);
  feSol.write("test");
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
