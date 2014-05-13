#include <iostream>
#include <cmath>

#include "SmallFem.h"

#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "FormulationSteadyWave.h"
#include "System.h"
#include "FEMSolution.h"

using namespace std;

// Data //
static const double Pi         = atan(1.0) * 4;
static const double paramaille = 6.0;
static const double nm         = 1000;
static const double lambda0    = 1000000.00000000;
static const double theta0     = 0.0 * Pi/180;
static const double phi0       = 0.0 * Pi/180;
static const double psi0       = 0.0 * Pi/180;
static const double a_lat      = 2000000.000000;
static const double period_x   = a_lat;
static const double period_y   = a_lat;
static const double period_z   = a_lat;
static const double PML_top    = 600000.000000 ;
static const double PML_bot    = 600000.000000;
static const double PML_lat    = 600000.000000;
static const double ro         = 500000.000000;
static const double eps_re_In  = 9.000000000;
static const double eps_im_In  = -1.000000000;
static const double eps_re_Out = 1.000000000;
static const double eps_im_Out = -0.000000000;

// Geomtrical Predicates //
static bool isInScat_In(fullVector<double>& xyz){
  return sqrt(xyz(0) * xyz(0) + xyz(1) * xyz(1) + xyz(2) * xyz(2)) < abs(ro);
}

static bool isInScat_Out(fullVector<double>& xyz){
  return
    (abs(xyz(0)) < period_x / 2) && // Not In PML
    (abs(xyz(1)) < period_y / 2) && //    |
    (abs(xyz(2)) < period_z / 2) && //    _
    !isInScat_In(xyz);              // Not In Scat
}

static bool isInDomain(fullVector<double>& xyz){
  return
    (abs(xyz(0)) < period_x / 2) && // Not In PML
    (abs(xyz(1)) < period_y / 2) && //    |
    (abs(xyz(2)) < period_z / 2);   //    _
}

static bool isInPMLxyz(fullVector<double>& xyz){
  return
    (abs(xyz(0)) > period_x / 2) &&
    (abs(xyz(1)) > period_y / 2) &&
    (abs(xyz(2)) > period_z / 2);
}

static bool isInPMLxz(fullVector<double>& xyz){
  return
    (abs(xyz(0)) > period_x / 2) &&
    (abs(xyz(1)) < period_y / 2) &&
    (abs(xyz(2)) > period_z / 2);
}

static bool isInPMLyz(fullVector<double>& xyz){
  return
    (abs(xyz(0)) < period_x / 2) &&
    (abs(xyz(1)) > period_y / 2) &&
    (abs(xyz(2)) > period_z / 2);
}

static bool isInPMLxy(fullVector<double>& xyz){
  return
    (abs(xyz(0)) > period_x / 2) &&
    (abs(xyz(1)) > period_y / 2) &&
    (abs(xyz(2)) < period_z / 2);
}

static bool isInPMLz(fullVector<double>& xyz){
  return
    (abs(xyz(0)) < period_x / 2) &&
    (abs(xyz(1)) < period_y / 2) &&
    (abs(xyz(2)) > period_z / 2);
}

static bool isInPMLy(fullVector<double>& xyz){
  return
    (abs(xyz(0)) < period_x / 2) &&
    (abs(xyz(1)) > period_y / 2) &&
    (abs(xyz(2)) < period_z / 2);
}

static bool isInPMLx(fullVector<double>& xyz){
  return
    (abs(xyz(0)) > period_x / 2) &&
    (abs(xyz(1)) < period_y / 2) &&
    (abs(xyz(2)) < period_z / 2);
}

static bool isInPML(fullVector<double>& xyz){
  return !isInDomain(xyz);
}

// Constants //
static const double mu0      =  4 * Pi * 100.0 * nm;
static const double epsilon0 =  8.854187817E-3 * nm;
static const double cel      =  1.0 / sqrt(epsilon0 * mu0);
static const double Freq     =  cel / lambda0;
static const double omega0   =  2.0 * Pi * Freq;
static const double k0       =  2.0 * Pi / lambda0;
static const double Ae       =  1.0;
static const double Ah       =  Ae * sqrt(epsilon0 / mu0);
static const double alpha0   =  k0 * sin(theta0) * cos(phi0);
static const double beta0    =  k0 * sin(theta0) * sin(phi0);
static const double gamma0   =  k0 * cos(theta0);
static const double Ex0      =  Ae * cos(psi0) * cos(theta0) * cos(phi0) -
                                Ae * sin(psi0) * sin(phi0);
static const double Ey0      =  Ae * cos(psi0) * cos(theta0) * sin(phi0) +
                                Ae * sin(psi0) * cos(phi0);
static const double Ez0      = -Ae * cos(psi0) * sin(theta0);
static const double Hx0      = -1 / (omega0 * mu0)
                                  * (beta0  * Ez0 - gamma0 * Ey0);
static const double Hy0      = -1 / (omega0 * mu0)
                                  * (gamma0 * Ex0 - alpha0 * Ez0);
static const double Hz0      = -1 / (omega0 * mu0)
                                  * (alpha0 * Ey0 - beta0  * Ex0);
static const double Pinc     =  0.5 * Ae * Ae
                                    * sqrt(epsilon0 / mu0) * cos(theta0);

// Functions //
static Complex Prop(fullVector<double>& xyz){
  return Ae * Complex(cos(alpha0 * xyz(0) + beta0 * xyz(1) + gamma0 * xyz(2)) ,
                      sin(alpha0 * xyz(0) + beta0 * xyz(1) + gamma0 * xyz(2)));
}

static fullVector<Complex> Einc(fullVector<double>& xyz){
  fullVector<Complex> ret(3);
  fullVector<Complex> zero(3);

  zero(0) = Complex(0, 0);
  zero(1) = Complex(0, 0);
  zero(2) = Complex(0, 0);

  if(isInPML(xyz))
    return zero;

  else{
    ret(0) = Ex0 * Prop(xyz);
    ret(1) = Ey0 * Prop(xyz);
    ret(2) = Ez0 * Prop(xyz);

    return ret;
  }
}

// PML //
static const double a_pml = 1.;
static const double b_pml = 1.;

static Complex sx(fullVector<double>& xyz){
  Complex ret = Complex(0, 0);

  if(isInScat_In(xyz))
    ret = Complex(1, 0);

  else if(isInScat_Out(xyz))
    ret = Complex(1, 0);

  else if(isInPMLxyz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLxz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLyz(xyz))
    ret = Complex(1, 0);

  else if(isInPMLxy(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLx(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLy(xyz))
    ret = Complex(1, 0);

  else if(isInPMLz(xyz))
    ret = Complex(1, 0);

  return ret;
}

static Complex sy(fullVector<double>& xyz){
  Complex ret = Complex(0, 0);

  if(isInScat_In(xyz))
    ret = Complex(1, 0);

  else if(isInScat_Out(xyz))
    ret = Complex(1, 0);

  else if(isInPMLxyz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLxz(xyz))
    ret = Complex(1, 0);

  else if(isInPMLyz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLxy(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLx(xyz))
    ret = Complex(1, 0);

  else if(isInPMLy(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLz(xyz))
    ret = Complex(1, 0);

  return ret;
}

static Complex sz(fullVector<double>& xyz){
  Complex ret = Complex(0, 0);

  if(isInScat_In(xyz))
    ret = Complex(1, 0);

  else if(isInScat_Out(xyz))
    ret = Complex(1, 0);

  else if(isInPMLxyz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLxz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLyz(xyz))
    ret = Complex(a_pml, -b_pml);

  else if(isInPMLxy(xyz))
    ret = Complex(1, 0);

  else if(isInPMLx(xyz))
    ret = Complex(1, 0);

  else if(isInPMLy(xyz))
    ret = Complex(1, 0);

  else if(isInPMLz(xyz))
    ret = Complex(a_pml, -b_pml);

  return ret;
}

static Complex Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

static Complex Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

static Complex Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// Materials //
static const Complex epsilon_In  = Complex(eps_re_In , eps_im_In);
static const Complex epsilon_Out = Complex(eps_re_Out, eps_im_Out);

static void epsilon(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  if(isInScat_In(xyz)){
    tensor(0, 0) = epsilon_In;
    tensor(1, 1) = epsilon_In;
    tensor(2, 2) = epsilon_In;
  }

  else if(isInScat_Out(xyz)){
    tensor(0, 0) = epsilon_Out;
    tensor(1, 1) = epsilon_Out;
    tensor(2, 2) = epsilon_Out;
  }

  else if(isInPML(xyz)){
    tensor(0, 0) = epsilon_Out * Lxx(xyz);
    tensor(1, 1) = epsilon_Out * Lzz(xyz);
    tensor(2, 2) = epsilon_Out * Lyy(xyz);
  }
}

static void epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  if(isInScat_In(xyz)){
    tensor(0, 0) = epsilon_Out;
    tensor(1, 1) = epsilon_Out;
    tensor(2, 2) = epsilon_Out;
  }

  else if(isInScat_Out(xyz)){
    tensor(0, 0) = epsilon_Out;
    tensor(1, 1) = epsilon_Out;
    tensor(2, 2) = epsilon_Out;
  }

  else if(isInPML(xyz)){
    tensor(0, 0) = epsilon_Out * Lxx(xyz);
    tensor(1, 1) = epsilon_Out * Lzz(xyz);
    tensor(2, 2) = epsilon_Out * Lyy(xyz);
  }
}

static void mu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  if(isInScat_In(xyz)){
    tensor(0, 0) = Complex(1, 0);
    tensor(1, 1) = Complex(1, 0);
    tensor(2, 2) = Complex(1, 0);
  }

  else if(isInScat_Out(xyz)){
    tensor(0, 0) = Complex(1, 0);
    tensor(1, 1) = Complex(1, 0);
    tensor(2, 2) = Complex(1, 0);
  }

  else if(isInPML(xyz)){
    tensor(0, 0) = Lxx(xyz);
    tensor(1, 1) = Lzz(xyz);
    tensor(2, 2) = Lyy(xyz);
  }
}

static void nu(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  if(isInScat_In(xyz)){
    tensor(0, 0) = Complex(1, 0);
    tensor(1, 1) = Complex(1, 0);
    tensor(2, 2) = Complex(1, 0);
  }

  else if(isInScat_Out(xyz)){
    tensor(0, 0) = Complex(1, 0);
    tensor(1, 1) = Complex(1, 0);
    tensor(2, 2) = Complex(1, 0);
  }

  else if(isInPML(xyz)){
    tensor(0, 0) = Complex(1, 0) / Lxx(xyz);
    tensor(1, 1) = Complex(1, 0) / Lzz(xyz);
    tensor(2, 2) = Complex(1, 0) / Lyy(xyz);
  }
}

// Source //
static fullVector<Complex> source(fullVector<double>& xyz){
  double kSquare = (omega0 / cel) * (omega0 / cel);

  fullVector<Complex> e = Einc(xyz);
  fullMatrix<Complex> eps(3, 3);
  fullMatrix<Complex> eps1(3, 3);

  epsilon(xyz, eps);
  epsilon1(xyz, eps1);

  fullVector<Complex> ret(3);

  ret(0) = kSquare * (eps(0, 0) - eps1(0, 0)) * e(0);
  ret(1) = kSquare * (eps(1, 1) - eps1(1, 1)) * e(1);
  ret(2) = kSquare * (eps(2, 2) - eps1(2, 2)) * e(2);

  return ret;
}

// FEM //
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

  // Full Scat Domain
  GroupOfElement Domain(msh);
  Domain.add(Scat_In);
  Domain.add(Scat_Out);

  // Full PML
  GroupOfElement PMLs(msh);
  PMLs.add(PMLxyz);
  PMLs.add(PMLxz);
  PMLs.add(PMLyz);
  PMLs.add(PMLxy);
  PMLs.add(PMLz);
  PMLs.add(PMLy);
  PMLs.add(PMLx);

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

  // Formulation //
  FunctionSpaceVector fs(All_domains, 2);
  FormulationSteadyWave<Complex>
    wave(All_domains, fs, (omega0 / cel), nu, epsilon, source);

  // System //
  System<Complex> sys;
  sys.addFormulation(wave);

  cout << "Assembling" << endl << flush;
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
  SmallFem::Keywords("-msh,-o,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
