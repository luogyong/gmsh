#include "SmallFem.h"
#include "Timer.h"
#include "System.h"
#include "SystemPETSc.h"
#include "SystemHelper.h"
#include "Interpolator.h"
#include "FormulationHelper.h"
#include "FormulationSilverMuller.h"
#include "FormulationSteadyWave.h"
#include "FormulationSteadySlow.h"

#include <iostream>

using namespace std;

static const int    scal = 0;
static const int    vect = 1;

static const Complex I  = Complex(0, 1);
static const double  Pi = M_PI;

static const Complex E0 = Complex(1, 0);
static const double  a  = 0.25;
static const double  b  = 0.25;
static const int     m  = 1;
static const int     n  = 1;

static const double  ky = m * Pi / a;
static const double  kz = n * Pi / b;
static       double  k;
static       Complex kx;

static       bool    isTE;
static       bool    is2D;

////////////////////
// Sources Fileds //
////////////////////
Complex             fSourceScal(fullVector<double>& xyz);
Complex                fAnaScal(fullVector<double>& xyz);
Complex               fZeroScal(fullVector<double>& xyz);
fullVector<Complex> fSourceVect(fullVector<double>& xyz);
fullVector<Complex>    fAnaVect(fullVector<double>& xyz);
fullVector<Complex>   fZeroVect(fullVector<double>& xyz);

Complex fSourceScal(fullVector<double>& xyz){
  const double y = xyz(1);
  const double z = xyz(2);

  if(is2D)
    return E0 * sin(ky * y);
  else
    return E0 * sin(ky * y) * sin(kz * z);
}

Complex fAnaScal(fullVector<double>& xyz){
  throw Exception("Not implemented");
}

Complex fZeroScal(fullVector<double>& xyz){
  return Complex(0, 0);
}

fullVector<Complex> fSourceVect(fullVector<double>& xyz){
  const double y  = xyz(1);
  const double z  = xyz(2);

  fullVector<Complex> tmp(3);

  if(is2D){
    if(isTE){
      throw Exception("Not implemented");
    }
    else{
      // TMm 2D
      tmp(0) = -E0 * I * ky / k * sin(ky * y);
      tmp(1) = +E0 *     kx / k * cos(ky * y);
      tmp(2) = Complex(0, 0);
    }
  }

  else{
    if(isTE){
      // TEmn 3D
      tmp(0) = Complex(0, 0);
      tmp(1) = -E0 * cos(ky * y) * sin(kz * z);
      tmp(2) = +E0 * sin(ky * y) * cos(kz * z);
    }

    else{
      // TMmn 3D
      tmp(0) = E0                                 * sin(ky * y) * sin(kz * z);
      tmp(1) = E0 * (I * kx * ky) / (k*k - kx*kx) * cos(ky * y) * sin(kz * z);
      tmp(2) = E0 * (I * kx * kz) / (k*k - kx*kx) * sin(ky * y) * cos(kz * z);
    }
  }

  return tmp;
}

fullVector<Complex> fAnaVect(fullVector<double>& xyz){
  const double x  = xyz(0);
  const double y  = xyz(1);
  const double z  = xyz(2);

  fullVector<Complex> tmp(3);

  if(is2D){
    if(isTE){
      throw Exception("Not implemented");
    }
    else{
      // TMm 2D
      tmp(0) = -E0 * I * ky / k * sin(ky * y) * exp(I * kx * x);
      tmp(1) = +E0 *     kx / k * cos(ky * y) * exp(I * kx * x);
      tmp(2) = Complex(0, 0);
    }
  }

  else{
    if(isTE){
      // TEmn 3D
      tmp(0) = Complex(0, 0);
      tmp(1) = -E0 * cos(ky * y) * sin(kz * z) * exp(I * kx * x);
      tmp(2) = +E0 * sin(ky * y) * cos(kz * z) * exp(I * kx * x);
    }

    else{
      // TMmn 3D
      tmp(0) = E0                           * sin(ky*y) * sin(kz*z) * exp(I*kx*x);
      tmp(1) = E0*(I*kx*ky) / (k*k - kx*kx) * cos(ky*y) * sin(kz*z) * exp(I*kx*x);
      tmp(2) = E0*(I*kx*kz) / (k*k - kx*kx) * sin(ky*y) * cos(kz*z) * exp(I*kx*x);
    }
  }

  return tmp;
}

fullVector<Complex> fZeroVect(fullVector<double>& xyz){
  fullVector<Complex> tmp(3);

  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(0, 0);
  tmp(2) = Complex(0, 0);

  return tmp;
}

/////////////
// Helpers //
/////////////
void getKx(void);

void ana(Complex (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<Complex>& eval);

void ana(fullVector<Complex> (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<Complex>& eval);

double modulusSquare(Complex a);

double l2Norm(const fullMatrix<Complex>& val);
double l2Norm(const fullMatrix<Complex>& valA, const fullMatrix<Complex>& valB);
double l2Error(const fullMatrix<Complex>& fem, const fullMatrix<Complex>& ana);

/////////////
// Compute //
/////////////
void compute(const Options& option){
  // Timers
  Timer tAssembly;
  Timer tSolve;

  // Get Type //
  int type;
  if(option.getValue("-type")[1].compare("scalar") == 0){
    cout << "Scalar Waveguide" << endl << flush;
    type = scal;
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    cout << "Vetorial Waveguide" << endl << flush;
    type = vect;
  }

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  k                  = atof(option.getValue("-k")[1].c_str());

  // Get Mode //
  string mode = option.getValue("-mode")[1];
  if(mode.compare("te") == 0)
    isTE = true;
  else if(mode.compare("tm") == 0)
    isTE = false;
  else
    throw Exception("Unknown mode %s", mode.c_str());

  // Is 2D ? //
  try{
    option.getValue("-2d");
    is2D = true;
  }
  catch(...){
    is2D = false;
  }

  // Compute kx //
  getKx();

  // Compute kInfinity for mode matching in silver-muller //
  Complex kInf;
  if(isTE)
    kInf = kx;
  else
    kInf = (k * k) / kx;

  cout << "Wavenumber: " << k            << endl
       << "Mode:       " << mode.c_str() << endl
       << "Order:      " << order        << endl << flush;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement infinity(msh.getFromPhysical(4));
  GroupOfElement   source(msh.getFromPhysical(5));
  GroupOfElement     zero(msh.getFromPhysical(6));
  GroupOfElement   volume(msh.getFromPhysical(7));

  // Full Domain //
  vector<const GroupOfElement*> domain(4);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &zero;
  domain[3] = &infinity;

  // Function Space //
  tAssembly.start();
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // Steady Wave Formulation //
  Formulation<Complex>*        wave;
  FormulationSilverMuller radiation(infinity, *fs, kInf);

  try{
    option.getValue("-slow");
    cout << "Slow Formulation" << endl << flush;
    wave = new FormulationSteadySlow<Complex>(volume, *fs, k);
  }

  catch(...){
    cout << "Fast Formulation" << endl << flush;
    wave = new FormulationSteadyWave<Complex>(volume, *fs, k);
  }

  // Solve //
  System<Complex> system;
  system.addFormulation(*wave);
  system.addFormulation(radiation);

  // Constraint
  if(fs->isScalar()){
    SystemHelper<Complex>::dirichlet(system, *fs, zero  , fZeroScal);
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceScal);
  }
  else{
    SystemHelper<Complex>::dirichlet(system, *fs, zero  , fZeroVect);
    SystemHelper<Complex>::dirichlet(system, *fs, source, fSourceVect);
  }

  // Assemble
  system.assemble();
  tAssembly.stop();
  cout << "Assembled: " << system.getSize()
       << "! (" << tAssembly.time() << " " << tAssembly.unit() << ")"
       << endl << flush;

  // Sove
  tSolve.start();
  system.solve();
  tSolve.stop();
  cout << "Solved!"
       << " (" << tSolve.time() << " " << tSolve.unit() << ")"
       << endl << flush;

  // Draw Solution //
  try{
    option.getValue("-nopos");
  }

  catch(...){
    cout << "Writing solution..." << endl << flush;

    FEMSolution<Complex> feSol;
    system.getSolution(feSol, *fs, volume);
    feSol.write("waveguide");
  }

  // L2 Error //
  try{
    // Visu mesh
    Mesh visu(option.getValue("-l2")[1]);

    // Info
    cout << "L2 Error..." << endl;

    // Get points
    fullMatrix<double> point;
    GroupOfElement     visuGoe = visu.getFromPhysical(7);
    visuGoe.getAllVertexCoordinate(point);

    // Analytic value
    fullMatrix<Complex> anaValue;
    if(type == scal)
      ana(fAnaScal, point, anaValue);
    else
      ana(fAnaVect, point, anaValue);

    // FEM value
    set<Dof> dof;
    fs->getKeys(volume, dof);

    set<Dof>::iterator    end = dof.end();
    set<Dof>::iterator     it = dof.begin();
    map<Dof, Complex>     sol;

    for(; it != end; it++)
      sol.insert(pair<Dof, Complex>(*it, 0));

    system.getSolution(sol, 0);

    vector<bool>        isValid;
    fullMatrix<Complex> femValue;
    Interpolator<Complex>::
      interpolate(volume, *fs, sol, point, femValue, isValid);

    // L2 Error
    cout << "L2 Error: " << scientific << l2Error(femValue, anaValue) << endl;
  }

  catch(...){
  }

  // Clean //
  delete wave;
  delete fs;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-type,-2d,-mode,-slow,-nopos,-l2");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}

////////////////////////////
// Helpers Implementation //
////////////////////////////
void getKx(void){
  kx = sqrt(Complex(k * k, 0) - (ky * ky) - (kz * kz));
}

void ana(Complex (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<Complex>& eval){
  // Alloc eval for Scalar Values //
  const size_t nPoint = point.size1();
  eval.resize(nPoint, 1);

  // Loop on point and evaluate f //
  fullVector<double> xyz(3);
  for(size_t i = 0; i < nPoint; i++){
    xyz(0) = point(i, 0);
    xyz(1) = point(i, 1);
    xyz(2) = point(i, 2);

    eval(i, 0) = f(xyz);
  }
}

void ana(fullVector<Complex> (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<Complex>& eval){
  // Alloc eval for Vectorial Values //
  const size_t nPoint = point.size1();
  eval.resize(nPoint, 3);

  // Loop on point and evaluate f //
  fullVector<double> xyz(3);
  fullVector<Complex> tmp(3);
  for(size_t i = 0; i < nPoint; i++){
    xyz(0) = point(i, 0);
    xyz(1) = point(i, 1);
    xyz(2) = point(i, 2);

    tmp = f(xyz);

    eval(i, 0) = tmp(0);
    eval(i, 1) = tmp(1);
    eval(i, 2) = tmp(2);
  }
}

double modulusSquare(Complex a){
  return (a.real() * a.real()) + (a.imag() * a.imag());
}

double l2Norm(const fullMatrix<Complex>& val){
  const size_t nPoint = val.size1();
  const size_t    dim = val.size2();

  double norm = 0;
  double modSquare = 0;

  for(size_t i = 0; i < nPoint; i++){
    modSquare = 0;

    for(size_t j = 0; j < dim; j++)
      modSquare += modulusSquare(val(i, j));

    norm += modSquare;
  }

  return sqrt(norm);
}

double l2Norm(const fullMatrix<Complex>& valA, const fullMatrix<Complex>& valB){
  const size_t nPoint = valA.size1();
  const size_t    dim = valA.size2();

  double norm = 0;
  double modSquare = 0;

  for(size_t i = 0; i < nPoint; i++){
    modSquare = 0;

    for(size_t j = 0; j < dim; j++)
      modSquare += modulusSquare(valA(i, j) - valB(i, j));

    norm += modSquare;
  }

  return sqrt(norm);
}

double l2Error(const fullMatrix<Complex>& fem, const fullMatrix<Complex>& ana){
  double anaNorm = l2Norm(ana);
  double femNorm = l2Norm(fem, ana);

  return femNorm / anaNorm;
}
